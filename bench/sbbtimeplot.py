#!/usr/bin/env python3

from io import StringIO
import csv
import matplotlib.pyplot as plt


TSV_DATA = """impl	hash_num	maker_wall_sec	query_wall_sec	maker_cpu_sec	query_cpu_sec	maker_max_mem_mb	query_max_mem_mb
cpp_classic	1	173.39	195.98	422.63	5951.23	22745.07	16388.11
rust_classic	1	24.49	135.27	559.96	2881.69	16620.60	16900.25
rust_blocked	1	23.54	165.43	550.30	3805.42	16621.29	16900.25
rust_bloomybloom	1	19.42	77.24	389.80	1735.93	16763.74	16900.77
cpp_classic	2	181.06	251.96	589.24	7636.95	22766.62	16388.32
rust_classic	2	33.79	211.68	837.38	4630.51	16621.22	16899.45
rust_blocked	2	30.80	191.50	745.84	4206.64	16621.01	16899.13
rust_bloomybloom	2	23.80	85.54	513.85	2067.20	16762.61	16899.68
cpp_classic	3	176.75	302.49	735.11	9076.66	22773.20	16387.92
rust_classic	3	45.42	224.04	1188.13	5926.42	16621.95	16900.05
rust_blocked	3	36.31	200.53	935.05	4480.28	16621.94	16899.91
rust_bloomybloom	3	25.47	95.82	556.22	2137.83	16764.72	16901.19
cpp_classic	4	181.41	359.47	879.50	10684.44	22718.44	16387.14
rust_classic	4	55.78	270.83	1504.10	7259.93	16620.64	16899.68
rust_blocked	4	37.97	198.50	992.39	4926.95	16620.90	16899.90
rust_bloomybloom	4	31.07	101.34	687.88	2317.74	16761.99	16901.61
cpp_classic	5	202.21	448.29	1133.44	12842.73	22724.57	16388.45
rust_classic	5	70.76	332.00	1935.79	8716.26	16620.11	16901.04
rust_blocked	5	42.85	237.41	1111.97	5661.53	16620.53	16901.06
rust_bloomybloom	5	38.30	112.54	817.17	2520.18	16768.79	16900.46
cpp_classic	6	213.92	478.39	1547.11	14054.64	22718.18	16388.07
rust_classic	6	71.88	375.47	2005.57	9969.87	16621.69	16901.66
rust_blocked	6	46.87	282.99	1176.38	6940.85	16621.72	16900.99
rust_bloomybloom	6	29.33	129.87	649.74	2794.84	16763.11	16901.50
cpp_classic	7	217.36	587.17	1500.53	17055.34	22749.80	16387.89
rust_classic	7	89.07	429.69	2464.23	11311.98	16620.38	16901.38
rust_blocked	7	42.19	238.49	1110.39	6123.50	16620.58	16901.23
rust_bloomybloom	7	29.85	113.76	682.49	2575.80	16760.66	16901.50
cpp_classic	8	234.70	577.05	1668.21	16935.20	22745.05	16388.07
rust_classic	8	109.52	490.40	2997.99	12600.99	16620.98	16901.56
rust_blocked	8	39.02	185.66	1012.35	4408.19	16620.89	16901.34
rust_bloomybloom	8	31.05	126.60	716.93	2887.94	16763.98	16901.12
cpp_classic	9	247.01	647.81	1855.29	19040.67	22741.25	16388.42
rust_classic	9	128.52	517.29	3563.42	14142.36	16621.07	16901.41
rust_blocked	9	44.66	238.93	1177.36	6218.69	16621.20	16901.12
rust_bloomybloom	9	32.16	128.26	746.37	3021.59	16766.77	16901.59
cpp_classic	10	261.76	669.36	1978.57	19894.99	22718.33	16388.43
rust_classic	10	128.22	586.45	3600.07	15623.49	16621.41	16901.87
rust_blocked	10	44.44	243.47	1179.85	6437.21	16621.26	16902.16
rust_bloomybloom	10	32.65	133.77	753.51	3146.36	16761.86	16900.46
"""


def load_rows(tsv_text: str):
    reader = csv.DictReader(StringIO(tsv_text), delimiter="\t")
    rows = []
    for row in reader:
        rows.append(
            {
                "impl": row["impl"],
                "hash_num": int(row["hash_num"]),
                "maker_wall_sec": float(row["maker_wall_sec"]),
                "query_wall_sec": float(row["query_wall_sec"]),
                "maker_cpu_sec": float(row["maker_cpu_sec"]),
                "query_cpu_sec": float(row["query_cpu_sec"]),
            }
        )
    return rows


def ordered_impls(rows):
    seen = []
    for row in rows:
        impl = row["impl"]
        if impl not in seen:
            seen.append(impl)
    return seen


def ordered_hash_nums(rows):
    return sorted({row["hash_num"] for row in rows})


def make_lookup(rows):
    return {(row["impl"], row["hash_num"]): row for row in rows}


def plot_metric(rows, maker_key, query_key, ylabel, title, output_file):
    impls = ordered_impls(rows)
    hash_nums = ordered_hash_nums(rows)
    lookup = make_lookup(rows)

    fig, ax = plt.subplots(figsize=(10, 6))

    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    impl_to_color = {
        impl: color_cycle[i % len(color_cycle)]
        for i, impl in enumerate(impls)
    }

    for impl in impls:
        color = impl_to_color[impl]
        maker_values = [lookup[(impl, h)][maker_key] for h in hash_nums]
        query_values = [lookup[(impl, h)][query_key] for h in hash_nums]

        ax.plot(
            hash_nums,
            maker_values,
            marker="o",
            linestyle="-",
            linewidth=2,
            color=color,
            label=f"{impl} index",
        )
        ax.plot(
            hash_nums,
            query_values,
            marker="s",
            linestyle="--",
            linewidth=2,
            color=color,
            label=f"{impl} query",
        )

    ax.set_title(title)
    ax.set_xlabel("Number of hash functions")
    ax.set_ylabel(ylabel)
    ax.set_xticks(hash_nums)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout()
    fig.savefig(output_file, dpi=200, bbox_inches="tight")


def main():
    rows = load_rows(TSV_DATA)

    plot_metric(
        rows=rows,
        maker_key="maker_cpu_sec",
        query_key="query_cpu_sec",
        ylabel="CPU time (s)",
        title="Index and query CPU time vs number of hash functions",
        output_file="bench_cpu_time.png",
    )

    plot_metric(
        rows=rows,
        maker_key="maker_wall_sec",
        query_key="query_wall_sec",
        ylabel="Wall clock time (s)",
        title="Index and query wall clock time vs number of hash functions",
        output_file="bench_wall_time.png",
    )

    plt.show()


if __name__ == "__main__":
    main()