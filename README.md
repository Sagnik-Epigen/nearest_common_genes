# nearest_common_genes

A Python package to find the **nearest genes** from two BED files using Ensembl annotations, and return the **common genes** between them. Uses [gget](https://github.com/pachterlab/gget) under the hood.

## 🧪 Installation

```bash
git clone [https://github.com/Sagnik-Epigen/nearest_common_genes.git]
cd nearest_common_genes
pip install .

# Usage
nearest-genes peaks1.bed peaks2.bed --assembly hg38 --upstream 2000 --downstream 2000 --output shared.tab

