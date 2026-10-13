import argparse
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--input", required=True)
p.add_argument("--out", default="top5_receptors.csv")
args = p.parse_args()

df = pd.read_csv(args.input)

receptors = df["DEG_as_receptor"].dropna().drop_duplicates().head(20)
receptors.to_frame(name="gene").to_csv(args.out, index=False)

print(f"Saved {len(receptors)} receptors to {args.out}")