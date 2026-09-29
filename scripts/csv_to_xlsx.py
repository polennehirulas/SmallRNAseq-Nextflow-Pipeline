#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser(description="Convert every CSV under a results directory to XLSX.")
parser.add_argument("--results", required=True)
args = parser.parse_args()

root = Path(args.results)
if not root.exists():
    raise SystemExit(f"Results directory does not exist: {root}")

for csv_path in root.rglob("*.csv"):
    xlsx_path = csv_path.with_suffix(".xlsx")
    try:
        df = pd.read_csv(csv_path)
        df.to_excel(xlsx_path, index=False)
        print(f"Converted: {csv_path} -> {xlsx_path}")
    except Exception as exc:
        print(f"WARNING: could not convert {csv_path}: {exc}")
