<p align="center">
    <a href="https://github.com/maxibor/pydamage"><img src="https://github.com/maxibor/pydamage/raw/master/docs/img/logo.png" alt="PyDamage"></a>
</p>

<p align="center">
    <a href="https://github.com/maxibor/pydamage/releases"><img src="https://img.shields.io/github/v/release/maxibor/pydamage?include_prereleases&label=version"/></a>
    <a href="https://github.com/maxibor/pydamage/actions"><img src="https://github.com/maxibor/pydamage/workflows/pydamage_ci/badge.svg"/></a>
    <a href="https://pydamage.readthedocs.io"><img src="https://readthedocs.org/projects/pydamage/badge/?version=latest"/></a>
    <a href="https://pypi.org/project/pydamage/"><img src="https://img.shields.io/badge/install%20with-pip-blue"/></a>
    <a href="https://anaconda.org/bioconda/pydamage"><img src="https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat"/></a>
</p>

# PyDamage

PyDamage is a fast Python tool to automate ancient DNA (aDNA) damage identification and estimation for assembled contigs in *de novo* metagenomic assembly.

After modeling ancient DNA damage using C to T / G to A transitions, PyDamage uses a likelihood ratio test and GLM accuracy prediction to discriminate between truly ancient and modern contigs originating from sample contamination.

---

## Installation

### With [`uv`](https://docs.astral.sh/uv/) (recommended)

```bash
uv pip install pydamage
# or run directly without manual install:
uv run --with pydamage pydamage --help
```

### With [bioconda](https://anaconda.org/bioconda/pydamage)

```bash
conda install -c bioconda pydamage
```

### With `pip`

```bash
pip install pydamage
```

---

## Development Setup

Clone the repository and synchronize the virtual environment using `uv`:

```bash
git clone https://github.com/maxibor/pydamage.git
cd pydamage

# Install all dependencies and development/profiling tools
uv sync --all-groups

# Run the test suite
uv run pytest
```

---

## Quick Start

```bash
# Run analysis on a BAM/SAM/CRAM file
pydamage --outdir result_directory analyze aligned.bam

# Run multi-threaded analysis with plotting
pydamage --outdir result_directory analyze -p 8 --plot aligned.bam

# Filter results by accuracy and q-value
pydamage filter result_directory/pydamage_results.csv

# Generate bin damage plot
pydamage binplot result_directory/pydamage_results.csv bin.fasta
```

---

## Benchmarks: 2026 Update vs Master

Direct head-to-head empirical measurements comparing the legacy `master` branch and the `2026_update` branch on the standard test dataset (`tests/data/aligned.bam`, 31 contigs, 17,504 reads) using Python 3.12:

### Pure Algorithmic Compute Time (In-Memory Processing)

| Worker Processes | `master` (Legacy) | `2026_update` (New) | Speedup |
| :--- | :--- | :--- | :--- |
| **`p=1` (Single process)** | 1,786.10 ms (17.4 contigs/s) | **179.73 ms (172.5 contigs/s)** | **9.94x faster** |
| **`p=2` (2 workers)** | 1,748.45 ms (17.7 contigs/s) | **114.93 ms (269.7 contigs/s)** | **15.21x faster** |
| **`p=4` (4 workers)** | 1,769.95 ms (17.5 contigs/s) | **68.67 ms (451.4 contigs/s)** | **25.78x faster** |
| **`p=8` (8 workers)** | 2,525.63 ms (12.3 contigs/s) | **57.85 ms (535.8 contigs/s)** | **43.66x faster** |

### End-to-End CLI Execution (Including CLI Startup & I/O)

| CLI Workload / Configuration | `master` (Legacy) | `2026_update` (New) | Speedup |
| :--- | :--- | :--- | :--- |
| **`pydamage analyze -p 1` (Single process)** | 2.83 s | **1.01 s** | **2.81x faster** |
| **`pydamage analyze -p 2` (2 workers)** | 2.86 s | **0.94 s** | **3.05x faster** |
| **`pydamage analyze -p 4` (4 workers)** | 2.85 s | **0.93 s** | **3.08x faster** |
| **`pydamage analyze -p 8` (8 workers)** | 3.73 s | **0.90 s** | **4.16x faster** |
| **`pydamage analyze -p 4 --plot` (Plotting)** | 4.37 s | **1.64 s** | **2.66x faster** |
| **`pydamage analyze -p 1 -r` (BAM Rescaling)** | 3.01 s | **1.17 s** | **2.57x faster** |
| **`pydamage cite` (Citation)** | 1.07 s | **0.81 s** | **1.32x faster** |
| **`pydamage --help` (CLI Help)** | 1.05 s | **0.80 s** | **1.31x faster** |

### Key Improvements
1. **Single-Pass BAM Traversal**: Base coverage, aligned read count, and damage transitions are computed in one traversal (`al_to_damage.get_damage`), eliminating 2 redundant BAM scans per contig.
2. **Persistent Worker BAM Handles**: Worker pools initialize a persistent `pysam.AlignmentFile` handle once upon startup, eliminating thousands of repeated file open/close cycles across contigs.
3. **Zero-Allocation In-Place Mutation**: Replaced temporary dictionary/list allocations in `damage_al` with in-place accumulators, removing >210,000 allocations per standard run.
4. **Vectorized NumPy PMF**: Replaced `np.vectorize` with pure array operations in `damage_model.fit`.
5. **Zero-Overhead Single Process Execution**: Fast-path for `-p 1` executes sequentially without `multiprocessing.Pool` spawn delay.

---

## Documentation

Full documentation is available at [pydamage.readthedocs.io](https://pydamage.readthedocs.io).

---

## Cite

PyDamage is published in PeerJ: [10.7717/peerj.11845](https://doi.org/10.7717/peerj.11845)

```bibtex
@article{borry_pydamage_2021,
    author = {Borry, Maxime and Hübner, Alexander and Rohrlach, Adam B. and Warinner, Christina},
    doi = {10.7717/peerj.11845},
    issn = {2167-8359},
    journal = {PeerJ},
    language = {en},
    month = {July},
    note = {Publisher: PeerJ Inc.},
    pages = {e11845},
    shorttitle = {PyDamage},
    title = {PyDamage: automated ancient damage identification and estimation for contigs in ancient DNA de novo assembly},
    url = {https://peerj.com/articles/11845},
    urldate = {2021-07-27},
    volume = {9},
    year = {2021}
}
```