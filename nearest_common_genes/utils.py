import pybedtools
import gget
import pandas as pd
import numpy as np
from collections import Counter

def bed_to_genes(bed_file, assembly, upstream=1000, downstream=1000):
    bed = pybedtools.BedTool(bed_file)
    gene_names = []

    for feature in bed:
        start = max(0, feature.start - upstream)
        end = feature.end + downstream
        coord = f"{feature.chrom}:{start}-{end}"

        try:
            result = gget.info(coord, assembly=assembly)
            if isinstance(result, pd.DataFrame) and not result.empty:
                genes = result['gene_name'].dropna().tolist()
                gene_names.extend(genes)
        except Exception as e:
            print(f"Warning: Failed to fetch for {coord}: {e}")

    return gene_names

def run_em_algorithm(gene_list, max_iter=100, tol=1e-4):
    """
    EM algorithm to estimate the most likely gene set based on frequency likelihood.
    This is a toy EM that reweights genes based on occurrence frequency until convergence.
    """
    gene_counts = Counter(gene_list)
    genes = list(gene_counts.keys())
    probs = np.array([gene_counts[g] for g in genes], dtype=np.float64)
    probs /= probs.sum()

    for _ in range(max_iter):
        # E-step: expected counts based on current probabilities
        expected = probs * len(gene_list)

        # M-step: update probabilities
        new_probs = expected / expected.sum()

        if np.linalg.norm(new_probs - probs) < tol:
            break
        probs = new_probs

    # Final selection: genes above median probability
    threshold = np.median(probs)
    selected = [genes[i] for i, p in enumerate(probs) if p >= threshold]
    return selected
