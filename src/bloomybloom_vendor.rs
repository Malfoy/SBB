use std::sync::Arc;

#[derive(Clone)]
enum FrozenShardStorage {
    Owned(Vec<u64>),
    Shared {
        words: Arc<Vec<u64>>,
        start: usize,
        len: usize,
    },
}

#[derive(Clone)]
pub(crate) struct FrozenBlockShard {
    storage: FrozenShardStorage,
    subblock_count: usize,
    words_per_subblock: usize,
    block_size: usize,
}

impl FrozenBlockShard {
    fn from_owned(
        words: Vec<u64>,
        subblock_count: usize,
        words_per_subblock: usize,
        block_size: usize,
    ) -> Self {
        Self {
            storage: FrozenShardStorage::Owned(words),
            subblock_count,
            words_per_subblock,
            block_size,
        }
    }

    fn from_shared(
        words: Arc<Vec<u64>>,
        start: usize,
        len: usize,
        subblock_count: usize,
        words_per_subblock: usize,
        block_size: usize,
    ) -> Self {
        Self {
            storage: FrozenShardStorage::Shared { words, start, len },
            subblock_count,
            words_per_subblock,
            block_size,
        }
    }

    pub(crate) fn as_slice(&self) -> &[u64] {
        match &self.storage {
            FrozenShardStorage::Owned(words) => words.as_slice(),
            FrozenShardStorage::Shared { words, start, len } => &words[*start..*start + *len],
        }
    }

    #[inline(always)]
    fn word_index(&self, subblock: usize, address: usize) -> usize {
        debug_assert!(subblock < self.subblock_count);
        debug_assert!(address < self.block_size);
        subblock * self.words_per_subblock + (address >> 6)
    }

    #[inline(always)]
    pub(crate) fn get(&self, subblock: usize, address: usize) -> bool {
        let word = self.as_slice()[self.word_index(subblock, address)];
        ((word >> (63 - (address & 63))) & 1) == 1
    }
}

#[derive(Clone)]
pub(crate) struct FrozenBloomFilter {
    pub(crate) filter: Vec<FrozenBlockShard>,
    pub(crate) block_size: usize,
    pub(crate) nb_blocks: usize,
    pub(crate) n_hashes: usize,
    pub(crate) block_size_mask: usize,
}

impl FrozenBloomFilter {
    pub(crate) fn word_len(&self) -> usize {
        self.filter.iter().map(|shard| shard.as_slice().len()).sum()
    }
}

mod inner {
    #![allow(dead_code)]

    include!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/vendor/bloomybloom/src/bloom.rs"
    ));

    pub(crate) fn into_wrapped_frozen(filter: FrozenBloomFilter) -> super::FrozenBloomFilter {
        let block_size = filter.block_size;
        let nb_blocks = filter.nb_blocks;
        let n_hashes = filter.n_hashes;
        let filter = filter
            .filter
            .into_iter()
            .map(|shard| {
                let words_per_subblock = shard.words_per_subblock;
                let subblock_count = shard.subblock_count;
                super::FrozenBlockShard::from_owned(
                    shard.words,
                    subblock_count,
                    words_per_subblock,
                    block_size,
                )
            })
            .collect();

        super::FrozenBloomFilter {
            filter,
            block_size,
            nb_blocks,
            n_hashes,
            block_size_mask: block_size - 1,
        }
    }
}

pub(crate) use inner::{BlockShard, BloomFilter as ConcurrentBloomFilter};

pub(crate) fn export_filter_words(filter: &FrozenBloomFilter) -> Vec<u64> {
    let mut words = Vec::with_capacity(filter.word_len());
    for shard in &filter.filter {
        words.extend_from_slice(shard.as_slice());
    }
    words
}

pub(crate) fn import_filter_from_words(
    words: Vec<u64>,
    block_size: usize,
    nb_blocks: usize,
    n_hashes: usize,
) -> FrozenBloomFilter {
    assert!(nb_blocks >= 1024);
    assert_eq!(nb_blocks % 1024, 0);

    let subblocks_per_shard = nb_blocks / 1024;
    let words_per_subblock = block_size.div_ceil(64);
    let words_per_shard = subblocks_per_shard * words_per_subblock;
    assert_eq!(words.len(), 1024 * words_per_shard);

    let words = Arc::new(words);
    let filter = (0..1024)
        .map(|shard_idx| {
            let start = shard_idx * words_per_shard;
            FrozenBlockShard::from_shared(
                Arc::clone(&words),
                start,
                words_per_shard,
                subblocks_per_shard,
                words_per_subblock,
                block_size,
            )
        })
        .collect();

    FrozenBloomFilter {
        filter,
        block_size,
        nb_blocks,
        n_hashes,
        block_size_mask: block_size - 1,
    }
}

pub(crate) fn import_filter_from_shards(
    shard_words: Vec<Vec<u64>>,
    block_size: usize,
    nb_blocks: usize,
    n_hashes: usize,
) -> FrozenBloomFilter {
    assert!(nb_blocks >= 1024);
    assert_eq!(nb_blocks % 1024, 0);
    assert_eq!(shard_words.len(), 1024);

    let subblocks_per_shard = nb_blocks / 1024;
    let words_per_subblock = block_size.div_ceil(64);
    let words_per_shard = subblocks_per_shard * words_per_subblock;

    let filter = shard_words
        .into_iter()
        .map(|words| {
            assert_eq!(words.len(), words_per_shard);
            FrozenBlockShard::from_owned(words, subblocks_per_shard, words_per_subblock, block_size)
        })
        .collect();

    FrozenBloomFilter {
        filter,
        block_size,
        nb_blocks,
        n_hashes,
        block_size_mask: block_size - 1,
    }
}

pub(crate) fn into_frozen(filter: inner::FrozenBloomFilter) -> FrozenBloomFilter {
    inner::into_wrapped_frozen(filter)
}
