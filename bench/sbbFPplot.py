#!/usr/bin/env python3

import csv
import io
from collections import defaultdict

import matplotlib.pyplot as plt

TSV_DATA = """tool\tbit_exp\tfp_log10
cpp_classic\t25\t7.573239
rust_classic\t25\t7.572865
rust_blocked\t25\t7.587370
rust_bloomybloom_s31\t25\t7.586905
rust_bloomybloom_s30\t25\t6.520655
rust_bloomybloom_s28\t25\t5.095678
rust_bloomybloom_s24\t25\t3.893873
cpp_classic\t26\t6.795888
rust_classic\t26\t6.795970
rust_blocked\t26\t6.840626
rust_bloomybloom_s31\t26\t6.828267
rust_bloomybloom_s30\t26\t5.088561
rust_bloomybloom_s28\t26\t2.482874
cpp_classic\t27\t5.956872
rust_classic\t27\t5.958324
rust_blocked\t27\t6.089524
rust_bloomybloom_s31\t27\t6.026318
rust_bloomybloom_s30\t27\t3.575880
cpp_classic\t28\t5.084612
rust_classic\t28\t5.090018
rust_blocked\t28\t5.434812
rust_bloomybloom_s31\t28\t5.223940
rust_bloomybloom_s30\t28\t2.082785
cpp_classic\t29\t4.199700
rust_classic\t29\t4.201452
rust_blocked\t29\t4.938505
rust_bloomybloom_s31\t29\t4.455043
rust_bloomybloom_s30\t29\t0.698970
cpp_classic\t30\t3.309843
rust_classic\t30\t3.323458
rust_blocked\t30\t4.559536
rust_bloomybloom_s31\t30\t3.756940
rust_bloomybloom_s30\t30\t0.000000
cpp_classic\t31\t2.431364
rust_classic\t31\t2.423246
rust_blocked\t31\t4.231138
rust_bloomybloom_s31\t31\t3.122871
cpp_classic\t32\t1.342423
rust_classic\t32\t1.397940
rust_blocked\t32\t3.923296
rust_bloomybloom_s31\t32\t2.593286
cpp_classic\t33\t0.845098
rust_classic\t33\t0.602060
rust_blocked\t33\t3.623353
rust_bloomybloom_s31\t33\t2.149219
rust_blocked\t34\t3.313445
rust_bloomybloom_s31\t34\t1.792392
rust_blocked\t35\t3.020361
rust_bloomybloom_s31\t35\t1.431364
"""

def main() -> None:
    reader = csv.DictReader(io.StringIO(TSV_DATA), delimiter="\t")
    series = defaultdict(list)

    for row in reader:
        tool = row["tool"]
        bit_exp = int(row["bit_exp"])
        fp_log10 = float(row["fp_log10"])
        series[tool].append((bit_exp, fp_log10))

    fig, ax = plt.subplots(figsize=(10, 6))

    all_x = set()
    for tool, points in sorted(series.items()):
        points.sort()
        xs = [x for x, _ in points]
        ys = [y for _, y in points]
        all_x.update(xs)
        ax.plot(xs, ys, marker="o", label=tool)

    xticks = sorted(all_x)
    ax.set_xticks(xticks)
    ax.set_xticklabels([rf"$2^{{{x}}}$" for x in xticks])

    ax.set_xlabel("Filter size in bits")
    ax.set_ylabel("False positives (log10)")
    ax.set_title("False positives vs filter size")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()

    plt.savefig("fp_vs_bits.png", dpi=200)
    plt.show()

if __name__ == "__main__":
    main()