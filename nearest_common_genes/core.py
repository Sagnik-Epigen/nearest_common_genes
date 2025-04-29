import argparse
import pandas as pd
from .utils import bed_to_genes, run_em_algorithm

def main():
    parser = argparse.ArgumentParser(description="Find common nearest genes from two BED files.")
    parser.add_argument("bed1", help="First BED file")
    parser.add_argument("bed2", help="Second BED file")
    parser.add_argument("--assembly", default="hg38", help="Genome assembly (default: hg38)")
    parser.add_argument("--upstream", type=int, default=1000, help="Upstream search distance (bp)")
    parser.add_argument("--downstream", type=int, default=1000, help="Downstream search distance (bp)")
    parser.add_argument("--output", default="common_genes.tab", help="Output .tab file path")

    args = parser.parse_args()

    genes1 = bed_to_genes(args.bed1, args.assembly, args.upstream, args.downstream)
    genes2 = bed_to_genes(args.bed2, args.assembly, args.upstream, args.downstream)

    optimized_genes1 = run_em_algorithm(genes1)
    optimized_genes2 = run_em_algorithm(genes2)

    common = sorted(set(optimized_genes1).intersection(set(optimized_genes2)))

    df = pd.DataFrame({'common_genes': common})
    df.to_csv(args.output, sep="\t", index=False)

    print(f"Saved {len(common)} common genes to {args.output}")
