//! kmer counter that uses a blocked bloom filter with hashing of minimizers to determine blocks
//! is tuned for a specific labtop (for now) and only supports up to 31-mers (and not optimal for
//! k<31)

mod input;
mod bloom;
mod unit_tests_one_day;
pub mod utils;
pub mod decyclers;
pub mod super_bitvec;
pub mod minimizers;

use input::{Hell};
use minimizers::{selected_mins_x_pos, MinimizerMode};
use decyclers::{Decycler};
use bloom::{BloomFilter, FrozenBloomFilter};
use utils::{hash_u128_to_u64, roll_u128_kmer, sum_vec_bool, xorshift_u64};
use packed_seq::{ChunkIt, PackedNSeqVec, PackedSeq, PackedSeqVec, PaddedIt, Seq, SeqVec, u32x8};
use std::env; //for backtrace
use rayon::prelude::*;
use clap::Parser;
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering};
use needletail::parse_fastx_file;
use rand::Rng;
use std::time::{Duration, Instant};
use libc::{RUSAGE_SELF, getrusage, rusage};
use std::io::{self, Write};
use seq_hash::{KmerHasher, NtHasher};

const INTRA_RECORD_CHUNK_KMERS: usize = 1 << 20;
const BITS_PER_GIB: usize = 1 << 33;
const HALF_FILL_LOG: f64 = std::f64::consts::LN_2;
const CARDINALITY_SKETCH_SIZE: usize = 1 << 13;
const HASH_SPACE_SIZE: f64 = (u32::MAX as f64) + 1.0;

#[derive(Clone, Copy, Debug, Default)]
struct PhaseStats {
    wall_seconds: f64,
    cpu_seconds: f64,
}

impl PhaseStats {
    fn cpu_usage_percent(self) -> f64 {
        if self.wall_seconds == 0.0 {
            0.0
        } else {
            100.0 * self.cpu_seconds / self.wall_seconds
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct PhaseTimer {
    wall_start: Instant,
    cpu_start_seconds: f64,
}

impl PhaseTimer {
    fn start() -> Self {
        Self {
            wall_start: Instant::now(),
            cpu_start_seconds: current_process_cpu_seconds(),
        }
    }

    fn finish(self) -> PhaseStats {
        PhaseStats {
            wall_seconds: self.wall_start.elapsed().as_secs_f64(),
            cpu_seconds: current_process_cpu_seconds() - self.cpu_start_seconds,
        }
    }
}

struct BottomSketchEstimator {
    values: Vec<u32>,
    sample_size: usize,
    bound: u32,
}

impl BottomSketchEstimator {
    fn new(sample_size: usize) -> Self {
        assert!(sample_size > 0);
        Self {
            values: Vec::with_capacity(sample_size * 2),
            sample_size,
            bound: u32::MAX,
        }
    }

    fn observe_hash(&mut self, hash: u32) {
        if hash == u32::MAX || hash >= self.bound {
            return;
        }
        self.values.push(hash);
        if self.values.len() >= self.sample_size * 2 {
            self.compact();
        }
    }

    fn compact(&mut self) {
        if self.values.is_empty() {
            self.bound = u32::MAX;
            return;
        }
        self.values.sort_unstable();
        self.values.dedup();
        if self.values.len() > self.sample_size {
            self.values.truncate(self.sample_size);
        }
        self.bound = if self.values.len() == self.sample_size {
            self.values[self.sample_size - 1]
        } else {
            u32::MAX
        };
    }

    fn merge(&mut self, other: Self) {
        self.values.extend(other.values);
        self.compact();
    }

    fn estimate_distinct_hashes(&mut self) -> f64 {
        self.compact();
        if self.values.len() < self.sample_size {
            self.values.len() as f64
        } else {
            let threshold = self.values[self.sample_size - 1] as f64 + 1.0;
            ((self.sample_size - 1) as f64) * HASH_SPACE_SIZE / threshold
        }
    }
}

fn current_process_cpu_seconds() -> f64 {
    unsafe {
        let mut usage: rusage = std::mem::zeroed();
        let status = getrusage(RUSAGE_SELF, &mut usage);
        assert_eq!(status, 0, "getrusage failed");
        timeval_to_seconds(usage.ru_utime) + timeval_to_seconds(usage.ru_stime)
    }
}

fn timeval_to_seconds(value: libc::timeval) -> f64 {
    value.tv_sec as f64 + value.tv_usec as f64 / 1_000_000.0
}

fn auto_size_bits_from_ram_gb(ram_gb: usize) -> usize {
    ram_gb
        .checked_mul(BITS_PER_GIB)
        .expect("requested RAM budget is too large")
}

fn max_superkmer_kmers(k: u16, m: u16) -> usize {
    k as usize - m as usize + 1
}

fn superkmer_smer_count(k: u16, m: u16, s: u16) -> usize {
    max_superkmer_kmers(k, m) + k as usize - s as usize
}

fn max_superkmer_address_capacity(k: u16, m: u16, s: u16, n_hashes: usize) -> usize {
    superkmer_smer_count(k, m, s) * n_hashes
}

fn optimal_hash_count(block_size: usize, k: u16, m: u16, s: u16) -> usize {
    let smer_count = superkmer_smer_count(k, m, s) as f64;
    ((block_size as f64 * HALF_FILL_LOG) / smer_count)
        .round()
        .max(1.0) as usize
}

fn optimal_smer_length(block_size: usize, k: u16, m: u16, n_hashes: usize) -> u16 {
    let max_superkmer = max_superkmer_kmers(k, m) as f64;
    let raw_s = k as f64 + max_superkmer - (block_size as f64 * HALF_FILL_LOG) / n_hashes as f64;
    let max_supported_s = usize::min(k as usize, 61) as f64;
    raw_s
        .round()
        .clamp(1.0, max_supported_s) as u16
}

fn default_smer_length(k: u16) -> u16 {
    k.saturating_sub(10).max(1)
}

fn infer_block_size_bits(size: usize, input_cardinality: usize, k: u16, m: u16) -> usize {
    let raw_block_size =
        (size as f64 * max_superkmer_kmers(k, m) as f64 / input_cardinality as f64).max(1.0);
    let max_block_size = usize::max(1, size / 1024) as f64;
    let block_size_exponent = raw_block_size
        .clamp(1.0, max_block_size)
        .log2()
        .round()
        .clamp(0.0, max_block_size.log2()) as u32;
    1usize << block_size_exponent
}

fn correct_hash_space_cardinality(distinct_hashes: f64) -> f64 {
    let saturated = (distinct_hashes / HASH_SPACE_SIZE).clamp(0.0, 1.0 - f64::EPSILON);
    -HASH_SPACE_SIZE * (-saturated).ln_1p()
}

fn add_valid_hashes_to_sketch(
    estimator: &mut BottomSketchEstimator,
    hashes: PaddedIt<impl ChunkIt<u32x8>>,
) {
    let lane_len = hashes.it.len();
    let valid_hash_count = lane_len * 8 - hashes.padding;
    let mut seen = 0usize;
    for lane_hashes in hashes.it {
        let remaining = valid_hash_count.saturating_sub(seen);
        if remaining == 0 {
            break;
        }
        let take = remaining.min(8);
        for &hash in &lane_hashes.as_array_ref()[..take] {
            estimator.observe_hash(hash);
        }
        seen += 8;
    }
}

fn estimate_input_cardinality(filename: &str, k: usize, chunk_size: usize) -> (usize, PhaseStats) {
    let timer = PhaseTimer::start();
    let reader = parse_fastx_file(filename).expect("valid path/file");
    let chunked_lines = Hell {
        fxreader: reader,
        chunk_size,
    };
    let mut estimator = BottomSketchEstimator::new(CARDINALITY_SKETCH_SIZE);
    for chunk in chunked_lines {
        let chunk_estimator = chunk
            .into_par_iter()
            .filter(|line| line.len() >= k)
            .flat_map_iter(|line| {
                split_sequence_for_parallelism(line, k as u16, INTRA_RECORD_CHUNK_KMERS).into_iter()
            })
            .fold(
                || (BottomSketchEstimator::new(CARDINALITY_SKETCH_SIZE), <NtHasher>::new(k)),
                |(mut local_estimator, hasher), line| {
                    let sequence = PackedNSeqVec::from_ascii(&line);
                    let hashes = hasher.hash_valid_kmers_simd(sequence.as_slice(), 1);
                    add_valid_hashes_to_sketch(&mut local_estimator, hashes);
                    (local_estimator, hasher)
                },
            )
            .map(|(local_estimator, _)| local_estimator)
            .reduce(
                || BottomSketchEstimator::new(CARDINALITY_SKETCH_SIZE),
                |mut left, right| {
                    left.merge(right);
                    left
                },
            );
        estimator.merge(chunk_estimator);
    }

    let distinct_hashes = estimator.estimate_distinct_hashes();
    let estimated_cardinality = correct_hash_space_cardinality(distinct_hashes);
    (estimated_cardinality.round() as usize, timer.finish())
}

fn flush_stdout_and_exit(code: i32) -> ! {
    io::stdout().flush().expect("stdout flush failed");
    std::process::exit(code);
}


///taking care of all the needed command line arguments
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    ///path to the file to fill the bloom filter
    input: Option<String>,

    ///path to the indexed file to fill the bloom filter
    #[arg(long)]
    indexed_file: Option<String>,

    ///path to a file containing the sequences to be queried
    #[arg(long, default_value_t = String::from(""))]
    query_file: String,

    ///length of the kmers
    #[arg(short, long, default_value_t = 31)]
    k: u16,

    ///quality versus performance parameter, positive integer, usually between 1 and 3, a higher
    ///value will lead to less false positive but slower execution than higher values
    #[arg(short, long, default_value_t = 4)]
    b: u16,

    ///max amount of RAM (integer in GB) to be used, must be at least 1
    #[arg(long, default_value_t = 1)]
    ram: usize,

    ///number of threads to use (default 1)
    #[arg(short, long, default_value_t = String::from("1"))]
    threads: String,

    ///number of hashes for the bloom filter; if omitted while s is set, it is chosen
    ///automatically from block-size occupancy
    #[arg(short, long)]
    n_hashes: Option<usize>,

    ///length of the s-mers, needs to be inferior or equal to k; if omitted while n_hashes is set,
    ///it is chosen automatically from block-size occupancy
    #[arg(short, long)]
    s: Option<u16>,

    ///length of the minimizers for grouping the kmers, has to be inferior to k and inferior to 32
    ///if omitted it is chosen automatically from the bloom geometry and b
    #[arg(short, long)]
    m: Option<u16>,

    ///size (in bits) of the bloom filter, expressed as a power of 2, overrides the ram
    ///parameter when provided
    #[arg(short, long)]
    size: Option<usize>,

    ///size (in bits) of each bloom block, expressed as a power of 2; if omitted while m is fixed,
    ///it can be inferred from the RAM budget and input cardinality
    #[arg(short, long)]
    block_size: Option<usize>,

    ///input cardinality in k-mers, used to infer block size when m is fixed and block size is not
    #[arg(long)]
    input_cardinality: Option<usize>,

    ///number of reads to be distributed in a row to each thread
    #[arg(long, default_value_t = 100)]
    sequential_fallback: usize,

    //to enable counting outputs, that take time outside of the actual algorithm, but allow to have
    //more insights on the bloom and the choice of settings
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    counting: bool,

    ///to disable all code referring to the bloom filter
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    no_bloom: bool,

    ///to disable all code after the parsing part, for bench purposes
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    only_parse: bool,

    ///to change the standard output to just sequence of numbers to be read by a benchmark programm
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    auto_bench: bool,

    ///enables using simd_minimizer double decycling minimizers instead of random
    ///is slower but slightly better in terms of false positives
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    no_simd_minimizer: bool,

    /// enables the open-closed minimizer order from the open-closed minimizer paper
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    open_closed_minimizer: bool,

    /// t-mer length used by open-closed minimizers; implies --open-closed-minimizer
    #[arg(long)]
    open_closed_t: Option<u16>,

}

pub fn main() {
    //for debug
    unsafe {
        env::set_var("RUST_BACKTRACE", "full");
    }
    let whole_run_start = Instant::now();
    //checking the arguments do make some sense
    let args = Args::parse();
    assert!(args.ram >= 1);

    //defining all variables constants that are based on the argument input
    let k: u16 = args.k;
    let requested_s = args.s.filter(|&s| s > 0);
    if let Some(s) = requested_s {
        assert!(args.k >= s); //cant have the s-mers be longer than the kmers
    }
    if let Some(n_hashes) = args.n_hashes {
        assert!(n_hashes > 0);
    }
    if let Some(input_cardinality) = args.input_cardinality {
        assert!(input_cardinality > 0);
    }
    if let Some(t) = args.open_closed_t {
        assert!(t > 0);
    }

    //number of threads allowed
    unsafe {
        env::set_var("RAYON_NUM_THREADS", args.threads.clone());
    }

    let filename = resolve_input_path(args.input.as_deref(), args.indexed_file.as_deref());
    let estimated_input_cardinality = if args.input_cardinality.is_none() {
        Some(estimate_input_cardinality(
            &filename,
            k as usize,
            args.sequential_fallback,
        ))
    } else {
        None
    };

    let m: u16;
    let n_hashes: usize;
    let s: u16;
    let mut size: usize;
    let mut block_size: usize;
    let mut nb_blocks: usize;
    let sequential_fallback: usize;
    let only_parse: bool;
    let no_bloom: bool;
    let no_simd_minimizer: bool;
    let minimizer_mode: MinimizerMode;
    let counting: bool;
    let auto_bench: bool;

    size = match args.size {
        Some(size_exponent) => 1 << size_exponent,
        None => auto_size_bits_from_ram_gb(args.ram),
    };
    match (args.m, args.block_size) {
        (Some(fixed_m), Some(block_size_exponent)) => {
            m = fixed_m;
            block_size = 1 << block_size_exponent;
        }
        (Some(fixed_m), None) => {
            m = fixed_m;
            let input_cardinality = args
                .input_cardinality
                .or_else(|| estimated_input_cardinality.map(|(estimate, _)| estimate))
                .expect("input cardinality estimation failed");
            block_size = infer_block_size_bits(size, input_cardinality, k, m);
        }
        (None, Some(block_size_exponent)) => {
            block_size = 1 << block_size_exponent;
            nb_blocks = size / block_size;
            m = (nb_blocks as f64).log2() as u16 / 2 + args.b;
        }
        (None, None) => {
            block_size = 1 << 13;
            nb_blocks = size / block_size;
            m = (nb_blocks as f64).log2() as u16 / 2 + args.b;
        }
    }
    nb_blocks = size / block_size;
    assert!(m < 32);
    assert!(m > 0);
    assert!(m <= k);
    if let Some(t) = args.open_closed_t {
        assert!(t <= m);
    }
    match (requested_s, args.n_hashes) {
        (Some(resolved_s), Some(resolved_hashes)) => {
            s = resolved_s;
            n_hashes = resolved_hashes;
        }
        (Some(resolved_s), None) => {
            s = resolved_s;
            n_hashes = optimal_hash_count(block_size, k, m, s);
        }
        (None, Some(resolved_hashes)) => {
            n_hashes = resolved_hashes;
            s = optimal_smer_length(block_size, k, m, n_hashes);
        }
        (None, None) => {
            s = default_smer_length(k);
            n_hashes = optimal_hash_count(block_size, k, m, s);
        }
    }
    assert!(s > 0);
    assert!(s < 62);
    sequential_fallback = args.sequential_fallback;
    only_parse = args.only_parse;
    no_bloom = args.no_bloom || args.only_parse;
    no_simd_minimizer = args.no_simd_minimizer;
    minimizer_mode = if let Some(t) = args.open_closed_t {
        MinimizerMode::OpenClosed { t }
    } else if args.open_closed_minimizer {
        MinimizerMode::OpenClosed { t: 5 }
    } else if no_simd_minimizer {
        MinimizerMode::Decycling
    } else {
        MinimizerMode::Simd
    };
    counting = args.counting;
    auto_bench = args.auto_bench;


    if no_bloom {
        size = 1024;
        block_size = 1;
        nb_blocks = 1024;
    }

    if !auto_bench {
        if let Some((estimated_input_cardinality, estimation_stats)) = estimated_input_cardinality {
            println!("estimated input cardinality : {estimated_input_cardinality}");
            println!("cardinality estimation time (s) : {}", estimation_stats.wall_seconds);
            println!("cardinality estimation cpu time (s) : {}", estimation_stats.cpu_seconds);
            println!("cardinality estimation cpu usage (%) : {}", estimation_stats.cpu_usage_percent());
        }
        println!("Parameters : ");
        println!("k : {k}, m: {m}, s : {s}, n_hashes : {n_hashes}");
        println!("bf size : {size}, block size {block_size}, nb_blocks {nb_blocks}");
        match minimizer_mode {
            MinimizerMode::Simd => println!("minimizer mode : simd"),
            MinimizerMode::Decycling => println!("minimizer mode : decycling"),
            MinimizerMode::OpenClosed { t } => {
                println!("minimizer mode : open-closed (t = {t})");
            }
        }
    }

    //we create the needed data structures to store everything
    let bloom = BloomFilter::new(size, n_hashes, k as usize, block_size, nb_blocks);

    //calculating decycling sets as well as its time
    let duration_overhead_decycling: Duration;
    let mut decycler_set: Decycler;
    let debut = Instant::now();
    if matches!(minimizer_mode, MinimizerMode::Decycling) {
        decycler_set = Decycler::new(m);
        decycler_set.compute_blocks();
    } else {
        decycler_set = Decycler::new(1);
    }
    duration_overhead_decycling = debut.elapsed();

    //anti optims variable
    let kmer_sum = AtomicU64::new(0);
    let inserted_kmer_count = AtomicU64::new(0);
    let queried_kmer_count = AtomicU64::new(0);

    //used to check for false negative rate at the end
    let false_neg_list: Mutex<Vec<PackedSeqVec>> = Mutex::new(Vec::new());

    let indexing_timer = PhaseTimer::start();
    {
        let reader = parse_fastx_file(&filename).expect("valid path/file");

        let chunked_lines = Hell {
            fxreader : reader,
            chunk_size : sequential_fallback,
        };

        chunked_lines.par_bridge().for_each(|chunk| {
            if only_parse {
                let block_lines_counter: usize = chunk.iter().map(Vec::len).sum();
                kmer_sum.fetch_add(block_lines_counter as u64, Ordering::Relaxed);
                return;
            }

            let mut sequence_chunks: Vec<Vec<u8>> = Vec::with_capacity(chunk.len());
            for line in chunk {
                if line.len() < k as usize {
                    continue;
                }

                if counting {
                    let dice_roll = rand::rng().random_range(0..5000);
                    if dice_roll == 0 {
                        let mut false_negs = false_neg_list.lock().unwrap();
                        false_negs.push(PackedSeqVec::from_ascii(&line));
                    }
                }

                sequence_chunks.extend(split_sequence_for_parallelism(
                    line,
                    k,
                    INTRA_RECORD_CHUNK_KMERS,
                ));
            }

            sequence_chunks.into_par_iter().for_each_init(
                || vec![0; max_superkmer_address_capacity(k, m, s, n_hashes)],
                |all_addresses, line| {
                    let sequence = PackedSeqVec::from_ascii(&line);
                    inserted_kmer_count.fetch_add(
                        (sequence.len() + 1 - k as usize) as u64,
                        Ordering::Relaxed,
                    );
                    let local_kmer_sum = handle_sequence(
                        &bloom,
                        sequence,
                        k,
                        m,
                        nb_blocks,
                        no_bloom,
                        all_addresses,
                        &decycler_set,
                        minimizer_mode,
                        s,
                    );
                    if no_bloom {
                        kmer_sum.fetch_add(local_kmer_sum, Ordering::Relaxed);
                    }
                },
            );
        })
    }
    let indexing_stats = indexing_timer.finish();
    let frozen_bloom = bloom.into_frozen();

    let query_timer = PhaseTimer::start();
    if args.query_file != "" {
        //this means that we do have to query
        let query_counter = AtomicU64::new(0);
        let positive_query_counter = AtomicU64::new(0);

        let reader = parse_fastx_file(&args.query_file).expect("valid path/file");

        let chunked_lines = Hell {
            fxreader : reader,
            chunk_size : sequential_fallback,
        };

        chunked_lines.par_bridge().for_each(|chunk| {
            let sequence_chunks: Vec<Vec<u8>> = chunk
                .into_iter()
                .filter(|line| line.len() >= k as usize)
                .flat_map(|line| split_sequence_for_parallelism(line, k, INTRA_RECORD_CHUNK_KMERS))
                .collect();

            sequence_chunks.into_par_iter().for_each(|line| {
                let sequence = PackedSeqVec::from_ascii(&line);
                let presence_vec = if s <= 31 {
                    frozen_bloom.check_sequence(sequence, k, m, s, &decycler_set, minimizer_mode)
                } else {
                    frozen_bloom.check_sequence_u128(sequence, k, m, s, &decycler_set, minimizer_mode)
                };
                let local_count = presence_vec.len() as u64;
                let local_pos_count = sum_vec_bool(&presence_vec) as u64;
                query_counter.fetch_add(local_count, Ordering::Relaxed);
                positive_query_counter.fetch_add(local_pos_count, Ordering::Relaxed);
                queried_kmer_count.fetch_add(local_count, Ordering::Relaxed);
            });
        });

        if !auto_bench {
            let q_count = query_counter.load(Ordering::Relaxed);
            let q_count_pos = positive_query_counter.load(Ordering::Relaxed);
            println!("Number of kmer queried : {q_count}");
            println!("Number of positives : {q_count_pos}");
        }
    }
    let query_stats = query_timer.finish();

    let (false_negative_rate, false_positive_rate) = if counting {
        let false_negs = false_neg_list.lock().unwrap().to_vec();
        frozen_bloom.count_false_bloom(false_negs, k, m, s, &decycler_set, minimizer_mode)
    } else {
        (0.0, 0.0)
    };
    let whole_run_duration = whole_run_start.elapsed();
    if !auto_bench {
        let inserted_kmers = inserted_kmer_count.load(Ordering::Relaxed);
        let queried_kmers = queried_kmer_count.load(Ordering::Relaxed);
        let whole_run_ns = whole_run_duration.as_nanos() as f64;
        println!("inserted kmers : {inserted_kmers}");
        if inserted_kmers > 0 {
            println!("ns per inserted kmer : {}", whole_run_ns / inserted_kmers as f64);
        } else {
            println!("ns per inserted kmer : N/A");
        }
        if queried_kmers > 0 {
            println!("ns per queried kmer : {}", whole_run_ns / queried_kmers as f64);
        } else {
            println!("ns per queried kmer : N/A");
        }
        println!("total indexing time (s) : {}", indexing_stats.wall_seconds);
        println!("index cpu time (s) : {}", indexing_stats.cpu_seconds);
        println!("index cpu usage (%) : {}", indexing_stats.cpu_usage_percent());
        println!("total query time (s) : {}", query_stats.wall_seconds);
        println!("query cpu time (s) : {}", query_stats.cpu_seconds);
        println!("query cpu usage (%) : {}", query_stats.cpu_usage_percent());
    }
    
    //to prevent optims
    //printing only a line for the benchmark evaluating programm if option --auto-bench if on
    if auto_bench {
        write_auto_bench_stdout(
            no_bloom, 
            frozen_bloom,
            nb_blocks,
            block_size,
            false_negative_rate,
            false_positive_rate,
            duration_overhead_decycling,
            )
    }
    else if counting {
        if !no_bloom {
            let (n_z_bloom, max_bloom, median_bloom, average_bloom, fill_counter) = frozen_bloom.count_it_all();
            let n_z_bloom_rate: f64 = n_z_bloom as f64/nb_blocks as f64;
            let max_bloom_rate: f64 = max_bloom as f64/block_size as f64;
            let median_bloom_rate: f64 = median_bloom as f64/block_size as f64;
            let average_bloom_rate: f64 = average_bloom as f64/block_size as f64;
            let overfilled_rate: f64 = fill_counter as f64/n_z_bloom as f64;

            println!("Non zero bf amount : {n_z_bloom}");
            println!("Non zero bloom filter block rates : {n_z_bloom_rate}");
            println!("Max bloom fill rate : {max_bloom_rate}");
            println!("Median fill rate : {median_bloom_rate}");
            println!("Average fill rate : {average_bloom_rate}");
            println!("Overfilled rate : {overfilled_rate}");
        }
    }
    flush_stdout_and_exit(0);
}

fn resolve_input_path(input: Option<&str>, indexed_file: Option<&str>) -> String {
    match (input, indexed_file) {
        (Some(path), None) => path.to_owned(),
        (None, Some(path)) => path.to_owned(),
        (None, None) => panic!("an input file path is required"),
        (Some(_), Some(_)) => panic!("use either the positional input or --indexed-file, not both"),
    }
}

fn split_sequence_for_parallelism(
    sequence: Vec<u8>,
    k: u16,
    target_chunk_kmers: usize,
) -> Vec<Vec<u8>> {
    if sequence.len() < k as usize {
        return vec![sequence];
    }

    let total_kmers = sequence.len() + 1 - k as usize;
    if total_kmers <= target_chunk_kmers {
        return vec![sequence];
    }

    let mut chunks = Vec::with_capacity((total_kmers + target_chunk_kmers - 1) / target_chunk_kmers);
    let mut start_kmer = 0;
    while start_kmer < total_kmers {
        let end_kmer = usize::min(start_kmer + target_chunk_kmers, total_kmers);
        let chunk_start = start_kmer;
        let chunk_end = end_kmer + k as usize - 1;
        chunks.push(sequence[chunk_start..chunk_end].to_vec());
        start_kmer = end_kmer;
    }

    chunks
}

fn handle_sequence(
    bloom: &BloomFilter,
    original_sequence: PackedSeqVec,
    k: u16,
    m: u16,
    nb_blocks: usize,
    no_bloom: bool,
    all_addresses: &mut Vec<usize>,
    decycler_set: &Decycler,
    minimizer_mode: MinimizerMode,
    l: u16,
    ) -> u64 {
    if original_sequence.len() < k as usize {
        return 0;
    }
    let address_mask: usize = bloom.block_size-1;
    let mut local_kmer_sum: u64 = 0;
    let (super_kmers_positions, minimizer_values, sequence): (Vec<u32>, Vec<u64>, PackedSeqVec) =
        selected_mins_x_pos(original_sequence, k, m, decycler_set, minimizer_mode);

    //quick check that we don't have abherrent results
    assert!(super_kmers_positions.len()==minimizer_values.len(), 
        "Superkmers and minimizers have different length.");

    let mut kmer_number: usize = 0;
    //compute all hashes at once to g faster than computing them 1 by 1
    for i in 0..super_kmers_positions.len()-1 {
        //using minimizer hashing for now to be sure its not a source of problems, will see if
        //removing it doesn't break anything later
        let hashed_minimizer: u64 = xorshift_u64(minimizer_values[i])&(nb_blocks as u64-1);
        if no_bloom {
            //prevent optims
            local_kmer_sum = local_kmer_sum.wrapping_add(hashed_minimizer);
        } else {
            if l <= 31 {
                kmer_number = 
                    handle_super_kmer(super_kmers_positions[i], super_kmers_positions[i+1], &sequence, 
                    bloom, k, hashed_minimizer, kmer_number,
                    all_addresses, address_mask, l);
            } else {
                kmer_number = 
                    handle_super_kmer_u128(super_kmers_positions[i], super_kmers_positions[i+1], &sequence, 
                    bloom, k, hashed_minimizer, kmer_number,
                    all_addresses, address_mask, l);
            }
        }
    }
    //not forgetting the last element of the list
    let hashed_minimizer: u64 = xorshift_u64(minimizer_values[minimizer_values.len()-1])&(nb_blocks as u64-1);
    if no_bloom {
        //prevent optims
        local_kmer_sum = local_kmer_sum.wrapping_add(hashed_minimizer);
    } else {
        if l <= 31 {
            let _ = 
                handle_super_kmer(super_kmers_positions[super_kmers_positions.len()-1], 
                (sequence.len()+1-k as usize) as u32,
                &sequence, 
                bloom, k, hashed_minimizer,
                kmer_number, all_addresses, address_mask, l);
        } else {
            let _ = 
                handle_super_kmer_u128(super_kmers_positions[super_kmers_positions.len()-1], 
                (sequence.len()+1-k as usize) as u32,
                &sequence,
                bloom, k, hashed_minimizer,
                kmer_number, all_addresses, address_mask, l);
        }
    }

    //is here only to prevent optimisations in case no bloom filters
    local_kmer_sum
}

fn handle_super_kmer(start_pos: u32, end_pos: u32, sequence: &PackedSeqVec,
    bloom: &BloomFilter, 
    k: u16, hashed_minimizer: u64, 
    mut kmer_number: usize, all_addresses: &mut Vec<usize>, 
    address_mask: usize, l: u16) -> usize {
    let smer_count = end_pos as usize + (k - l) as usize - start_pos as usize;
    let required_addresses = smer_count * bloom.n_hashes;
    if all_addresses.len() < required_addresses {
        all_addresses.resize(required_addresses, 0);
    }
    let mut last_relevant_index: usize = 0;
    for j in (start_pos as usize)..(end_pos as usize) + (k-l) as usize{
        let smer: PackedSeq = sequence.slice(j..j+l as usize);
        let mut hash: u64 = xorshift_u64(smer.as_u64());


        for _i in 0..bloom.n_hashes {
            let address = hash as usize & address_mask;
            all_addresses[last_relevant_index] = address;
            last_relevant_index += 1;
            hash = xorshift_u64(hash);
        }

        kmer_number+=1;

    }
    let relevant_addresses = &mut all_addresses[..last_relevant_index];
    let blocknum: usize = (hashed_minimizer as usize)&1023;
    let subblocknum: usize = ((hashed_minimizer as usize)>>10)&((bloom.nb_blocks>>10)-1);
    let mut block = bloom.filter[blocknum].lock().unwrap();
    for address in relevant_addresses {
            if !block.get(subblocknum, *address) {
                block.set(subblocknum, *address, true);
            }
    }
    drop(block);
    kmer_number
}

fn handle_super_kmer_u128(start_pos: u32, end_pos: u32, sequence: &PackedSeqVec,
    bloom: &BloomFilter, 
    k: u16, hashed_minimizer: u64, 
    mut kmer_number: usize, _all_addresses: &mut Vec<usize>, 
    address_mask: usize, l: u16) -> usize {
    let smer_count = end_pos as usize + (k - l) as usize - start_pos as usize;
    let blocknum: usize = (hashed_minimizer as usize)&1023;
    let subblocknum: usize = ((hashed_minimizer as usize)>>10)&((bloom.nb_blocks>>10)-1);
    let mut block = bloom.filter[blocknum].lock().unwrap();

    let mut current_smer = sequence.slice(start_pos as usize..start_pos as usize + l as usize).as_u128();
    for offset in 0..smer_count {
        let mut hash = hash_u128_to_u64(current_smer);
        for _ in 0..bloom.n_hashes {
            let address = hash as usize & address_mask;
            if !block.get(subblocknum, address) {
                block.set(subblocknum, address, true);
            }
            hash = xorshift_u64(hash);
        }

        kmer_number += 1;
        if offset + 1 < smer_count {
            let next_base = sequence.as_slice().get(start_pos as usize + offset + l as usize);
            current_smer = roll_u128_kmer(current_smer, next_base, l);
        }
    }
    kmer_number
}


fn write_auto_bench_stdout(
    no_bloom : bool, 
    bloom: FrozenBloomFilter,
    nb_blocks: usize,
    block_size: usize,
    false_positive_rate: f64,
    false_negative_rate: f64,
    duration_overhead_decycling: Duration,
    ) {
    let mut print_string = String::new();
    //writes every number looked for by the benchmark programm in a single line
    //also does all the counting
    if !no_bloom {
        let (n_z_bloom, max_bloom, median_bloom, average_bloom, fill_counter) = bloom.count_it_all();
        let n_z_bloom_rate: f64 = n_z_bloom as f64/nb_blocks as f64;
        let max_bloom_rate: f64 = max_bloom as f64/block_size as f64;
        let median_bloom_rate: f64 = median_bloom as f64/block_size as f64;
        let average_bloom_rate: f64 = average_bloom as f64/block_size as f64;
        let overfilled_rate: f64 = fill_counter as f64/n_z_bloom as f64;

        print_string += 
            &format!("{n_z_bloom_rate}|{max_bloom_rate}|{average_bloom_rate}|{median_bloom_rate}|{overfilled_rate}");
    } else {
        print_string += &format!("0|0|0|0");
    }


    //false negatives and false potitives rates
    print_string += &format!("|{:.3}|{:.3}", false_positive_rate, false_negative_rate);

    //duration of overhead decycling set calculation
    print_string += &format!("|{}", duration_overhead_decycling.as_secs());

    println!("{print_string}");
}

#[cfg(test)]
mod tests {
    use super::{
        PhaseStats, auto_size_bits_from_ram_gb, handle_sequence, handle_super_kmer,
        handle_super_kmer_u128, infer_block_size_bits, max_superkmer_address_capacity,
        optimal_hash_count, optimal_smer_length, default_smer_length, resolve_input_path,
        split_sequence_for_parallelism, BloomFilter, Decycler,
        estimate_input_cardinality, MinimizerMode,
    };
    use packed_seq::{PackedSeqVec, SeqVec};
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_path(name: &str) -> PathBuf {
        let mut path = std::env::temp_dir();
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        path.push(format!("bloomybloom_main_{name}_{unique}_{}", std::process::id()));
        path
    }

    fn build_bloom(k: usize) -> BloomFilter {
        BloomFilter::new(1 << 20, 3, k, 1 << 10, 1 << 10)
    }

    fn build_decycler(m: u16) -> Decycler {
        let mut decycler = Decycler::new(m);
        decycler.compute_blocks();
        decycler
    }

    #[test]
    fn handle_sequence_returns_zero_for_short_inputs() {
        let bloom = build_bloom(5);
        let decycler = build_decycler(3);
        let sequence = PackedSeqVec::from_ascii(b"ACGT");
        let mut addresses = vec![0; 64];

        assert_eq!(
            handle_sequence(
                &bloom,
                sequence,
                5,
                3,
                1 << 10,
                false,
                &mut addresses,
                &decycler,
                MinimizerMode::Decycling,
                3,
            ),
            0
        );
    }

    #[test]
    fn handle_sequence_no_bloom_returns_non_zero_hash_sum() {
        let bloom = build_bloom(5);
        let decycler = build_decycler(3);
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGT");
        let mut addresses = vec![0; 64];

        let sum = handle_sequence(
            &bloom,
            sequence,
            5,
            3,
            1 << 10,
            true,
            &mut addresses,
            &decycler,
            MinimizerMode::Decycling,
            3,
        );
        assert_ne!(sum, 0);
    }

    #[test]
    fn handle_sequence_populates_bloom_for_u64_queries() {
        let bloom = build_bloom(5);
        let decycler = build_decycler(3);
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGT");
        let mut addresses = vec![0; 64];

        let _ = handle_sequence(
            &bloom,
            sequence.clone(),
            5,
            3,
            1 << 10,
            false,
            &mut addresses,
            &decycler,
            MinimizerMode::Decycling,
            3,
        );
        let results = bloom.check_sequence(sequence, 5, 3, 3, &decycler, MinimizerMode::Decycling);

        assert_eq!(results.len(), 4);
        assert!(results.iter().all(|present| *present));
    }

    #[test]
    fn handle_sequence_populates_bloom_for_u128_queries() {
        let bloom = build_bloom(40);
        let decycler = build_decycler(3);
        let sequence = PackedSeqVec::from_ascii(
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        );
        let mut addresses = vec![0; 256];

        let _ = handle_sequence(
            &bloom,
            sequence.clone(),
            40,
            3,
            1 << 10,
            false,
            &mut addresses,
            &decycler,
            MinimizerMode::Decycling,
            32,
        );
        let results =
            bloom.check_sequence_u128(sequence, 40, 3, 32, &decycler, MinimizerMode::Decycling);

        assert_eq!(results.len(), 9);
        assert!(results.iter().all(|present| *present));
    }

    #[test]
    fn handle_super_kmer_returns_updated_counter() {
        let bloom = build_bloom(5);
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGT");
        let mut addresses = vec![0; 64];

        let count = handle_super_kmer(0, 2, &sequence, &bloom, 5, 0, 7, &mut addresses, (1 << 10) - 1, 3);
        assert_eq!(count, 11);
    }

    #[test]
    fn handle_super_kmer_u128_returns_updated_counter() {
        let bloom = build_bloom(40);
        let sequence = PackedSeqVec::from_ascii(
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        );
        let mut addresses = vec![0; 256];

        let count = handle_super_kmer_u128(0, 2, &sequence, &bloom, 40, 0, 5, &mut addresses, (1 << 10) - 1, 32);
        assert_eq!(count, 15);
    }

    #[test]
    fn resolve_input_path_returns_positional_input() {
        let path = resolve_input_path(Some("reads.fa"), None);
        assert_eq!(path, "reads.fa");
    }

    #[test]
    fn resolve_input_path_returns_indexed_option_path() {
        let path = resolve_input_path(None, Some("indexed.fa"));
        assert_eq!(path, "indexed.fa");
    }

    #[test]
    #[should_panic(expected = "an input file path is required")]
    fn resolve_input_path_requires_one_input_source() {
        let _ = resolve_input_path(None, None);
    }

    #[test]
    #[should_panic(expected = "use either the positional input or --indexed-file, not both")]
    fn resolve_input_path_rejects_both_input_sources() {
        let _ = resolve_input_path(Some("reads.fa"), Some("indexed.fa"));
    }

    #[test]
    fn split_sequence_for_parallelism_keeps_short_sequences_intact() {
        let chunks = split_sequence_for_parallelism(b"ACGTAC".to_vec(), 7, 4);
        assert_eq!(chunks, vec![b"ACGTAC".to_vec()]);
    }

    #[test]
    fn split_sequence_for_parallelism_keeps_small_kmer_sets_intact() {
        let chunks = split_sequence_for_parallelism(b"ACGTAC".to_vec(), 4, 3);
        assert_eq!(chunks, vec![b"ACGTAC".to_vec()]);
    }

    #[test]
    fn split_sequence_for_parallelism_splits_into_overlapping_kmer_ranges() {
        let chunks = split_sequence_for_parallelism(b"ACGTACGTAC".to_vec(), 4, 3);
        assert_eq!(
            chunks,
            vec![
                b"ACGTAC".to_vec(),
                b"TACGTA".to_vec(),
                b"GTAC".to_vec(),
            ]
        );
    }

    #[test]
    fn phase_stats_cpu_usage_percent_handles_zero_wall_time() {
        let stats = PhaseStats {
            wall_seconds: 0.0,
            cpu_seconds: 12.0,
        };

        assert_eq!(stats.cpu_usage_percent(), 0.0);
    }

    #[test]
    fn phase_stats_cpu_usage_percent_scales_cpu_time() {
        let stats = PhaseStats {
            wall_seconds: 2.0,
            cpu_seconds: 6.0,
        };

        assert_eq!(stats.cpu_usage_percent(), 300.0);
    }

    #[test]
    fn auto_size_bits_from_ram_gb_maps_one_gib_correctly() {
        assert_eq!(auto_size_bits_from_ram_gb(1), 1 << 33);
    }

    #[test]
    fn auto_size_bits_from_ram_gb_maps_thirty_two_gib_correctly() {
        assert_eq!(auto_size_bits_from_ram_gb(32), 1 << 38);
    }

    #[test]
    fn optimal_hash_count_matches_half_fill_formula() {
        assert_eq!(optimal_hash_count(1 << 13, 31, 13, 28), 258);
    }

    #[test]
    fn optimal_smer_length_inverts_hash_count_formula() {
        assert_eq!(optimal_smer_length(1 << 13, 31, 13, 258), 28);
    }

    #[test]
    fn default_smer_length_is_k_minus_ten() {
        assert_eq!(default_smer_length(31), 21);
    }

    #[test]
    fn default_smer_length_clamps_to_one() {
        assert_eq!(default_smer_length(5), 1);
    }

    #[test]
    fn infer_block_size_bits_matches_formula_for_exact_power_of_two() {
        assert_eq!(infer_block_size_bits(1 << 38, 11 << 25, 31, 21), 1 << 13);
    }

    #[test]
    fn estimate_input_cardinality_tracks_canonical_kmers() {
        let path = temp_path("cardinality.fasta");
        fs::write(&path, b">seq1\nACGTACGT\n").unwrap();

        let (estimate, _) = estimate_input_cardinality(&path.to_string_lossy(), 3, 2);
        assert!(estimate >= 1);
        assert!(estimate <= 4);

        fs::remove_file(path).unwrap();
    }

    #[test]
    fn max_superkmer_address_capacity_matches_smer_count_times_hashes() {
        assert_eq!(max_superkmer_address_capacity(31, 24, 21, 315), 5670);
    }
}
