"""Таблицы и графики по сохранённым парным прогонам gkls_mco.py."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from gkls_mco import frontier, hv


def series(folder, budget, local):
    return np.genfromtxt(folder / f"solve-{budget}-{local:g}.csv", delimiter=",", skip_header=1, ndmin=2)


def normalized(values, row):
    return (values + 1) / np.array([row["scale1"], row["scale2"]])


def volume(trials, row, count):
    prefix = trials[:count, 3:5]
    return hv(normalized(prefix[np.isfinite(prefix[:, 0])], row)) / row["reference_hv"]


def save_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def summarize(out, rows):
    pairs, summaries = [], []
    for budget in sorted({int(r["budget"]) for r in rows}):
        group = [r for r in rows if r["budget"] == budget]
        one = {int(r["instance"]): r for r in group if r["r_local"] == 0}
        two = {int(r["instance"]): r for r in group if r["r_local"] > 0}
        if one.keys() != two.keys():
            raise ValueError("Incomplete pairs")
        for i, a in one.items():
            b = two[i]
            folder = out / f"instance-{i:03d}"
            n = int(min(a["trials"], b["trials"]))
            common = (volume(series(folder, budget, 1.5), b, n)
                      - volume(series(folder, budget, 0), a, n))
            pairs.append(dict(instance=i, budget=budget, common_trials=n, delta_common_hv=common,
                              delta_relative_hv=b["relative_hv"] - a["relative_hv"],
                              delta_igd_plus=b["igd_plus"] - a["igd_plus"],
                              delta_trials=b["trials"] - a["trials"]))
        sub = [p for p in pairs if p["budget"] == budget]
        for field in ["delta_relative_hv", "delta_common_hv", "delta_igd_plus", "delta_trials"]:
            d = np.array([p[field] for p in sub])
            q1, median, q3 = np.percentile(d, [25, 50, 75])
            summaries.append(dict(budget=budget, metric=field, q1=q1, median=median, q3=q3,
                                  positive=sum(d > 1e-10), zero=sum(abs(d) <= 1e-10),
                                  negative=sum(d < -1e-10)))
    save_csv(out / "pairs.csv", pairs)
    save_csv(out / "summary.csv", summaries)
    return pairs


def figures(out, rows, pairs):
    budget = max(int(r["budget"]) for r in rows)
    group = sorted((p for p in pairs if p["budget"] == budget), key=lambda p: p["delta_relative_hv"])
    selected = [group[k]["instance"] for k in sorted({0, len(group)//2, len(group)-1})]
    fig, axes = plt.subplots(2, len(selected), squeeze=False, figsize=(4.4*len(selected), 7.2))
    for column, i in enumerate(selected):
        row = next(r for r in rows if r["instance"] == i and r["budget"] == budget)
        folder = out / f"instance-{i:03d}"
        ref = np.loadtxt(folder / f"grid-{int(row['grid'])}.csv", delimiter=",", skiprows=1, ndmin=2)
        ax, curve = axes[:, column]
        ax.scatter(ref[:, 0], ref[:, 1], s=4, color="0.7", label="Сетка")
        for loc, color, label in [(0, "#ae4a22", "Один r"), (1.5, "#1e668c", "Два r")]:
            trials = series(folder, budget, loc)
            front = frontier(trials[np.isfinite(trials[:, 3]), 3:5])
            ax.scatter(front[:, 0], front[:, 1], s=13, color=color, label=label)
            counts = np.unique(np.geomspace(1, len(trials), 70).astype(int))
            curve.plot(counts, [volume(trials, row, n) for n in counts], color=color, label=label)
        ax.set(title=f"Экземпляр №{i}", xlabel="Первый критерий", ylabel="Второй критерий")
        curve.set(xlabel="Испытания", ylabel="HV / HV сетки", xscale="log")
        ax.legend()
        curve.axhline(0.95, linewidth=.8, linestyle="--", color="0.4")
    fig.tight_layout()
    fig.savefig(out / "fronts.png", dpi=160)
    plt.close(fig)

    pilots = [i for i in (1, 50, 100) if (out / f"instance-{i:03d}" / "map.csv").exists()]
    if not pilots:
        return
    fig, axes = plt.subplots(1, len(pilots), squeeze=False, figsize=(4.2*len(pilots), 4))
    for ax, i in zip(axes[0], pilots):
        folder = out / f"instance-{i:03d}"
        grid = np.loadtxt(folder / "map.csv", delimiter=",", skiprows=1)
        x = np.unique(grid[:, 0])
        a, b = [grid[:, k].reshape(len(x), len(x)).T for k in (2, 3)]
        ax.contourf(x, x, (a <= 0) & (b <= 0), levels=[-.5,.5,1.5], colors=["#eeeeee", "#cbe2c9"])
        for z, color in [(a, "#536b89"), (b, "#784c7a")]:
            if z.min() < 0 < z.max():
                ax.contour(x, x, z, levels=[0], colors=color, linewidths=1)
        row = next(r for r in rows if r["instance"] == i)
        ref = np.loadtxt(folder / f"grid-{int(row['grid'])}.csv", delimiter=",", skiprows=1, ndmin=2)
        ax.scatter(ref[:, 2], ref[:, 3], s=3, color="0.3", label="Парето по сетке")
        for loc, color, label in [(0, "#ae4a22", "Один r"), (1.5, "#1e668c", "Два r")]:
            trials = series(folder, budget, loc)
            feasible = trials[np.isfinite(trials[:, 3])]
            ordered = feasible[np.lexsort((feasible[:, 4], feasible[:, 3]))]
            if len(ordered):
                chosen = ordered[np.r_[True, ordered[1:, 4] < np.minimum.accumulate(ordered[:-1, 4])]]
                ax.scatter(chosen[:, 1], chosen[:, 2], s=10, color=color, label=label)
        ax.set(title=f"№{i} · допустимо {row['feasible_fraction']:.1%}", xlabel="x₁", ylabel="x₂",
               aspect="equal", xlim=(-1,1), ylim=(-1,1))
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "domains.png", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    out = parser.parse_args().results
    with (out / "metrics.csv").open() as f:
        rows = [{k: float(v) if v else np.nan for k, v in row.items()} for row in csv.DictReader(f)]
    figures(out, rows, summarize(out, rows))


if __name__ == "__main__":
    main()
