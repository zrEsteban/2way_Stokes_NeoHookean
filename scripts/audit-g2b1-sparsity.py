#!/usr/bin/env python3
"""Compare H^T W Z H connectivity with the current cell-based FE pattern."""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import g2a_dual_transfer as dual


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    args = parser.parse_args()
    rows, *_ = dual.build_operator(args.case)
    _, cells, stable = dual.gmsh_volume(args.case / "dealiiSolid/solid.msh")
    cell_pattern = set()
    for cell in cells:
        ids = [stable[node] for node in cell]
        cell_pattern.update((i, j) for i in ids for j in ids)
    missing = set()
    for row in rows:
        ids = list(row)
        missing.update((i, j) for i in ids for j in ids
                       if (i, j) not in cell_pattern)
    print(f"rows_H={len(rows)}")
    print(f"missing_ordered_scalar_pairs={len(missing)}")
    print("sample=" + " ".join(f"({i},{j})" for i, j in sorted(missing)[:10]))
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
