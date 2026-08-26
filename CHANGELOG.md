# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased] - 2026 Rewrite

### Added
- **Modern Package & Dependency Management with `uv`**:
  - Replaced legacy `setup.py` and `environment.yml` with a standard PEP 517/621 [`pyproject.toml`](pyproject.toml) using the `hatchling` build backend.
  - Added dependency groups (`dev`, `profile`) including `pytest` and `memray`, with a reproducible lockfile [`uv.lock`](uv.lock).
  - Modernized GitHub Actions CI to use `astral-sh/setup-uv@v5`.
- **Comprehensive Unit & Integration Test Suite (Pytest)**:
  - Added [`tests/test_cli.py`](tests/test_cli.py) covering all Click CLI subcommands (`analyze`, `filter`, `binplot`, `cite`, and mutually exclusive options).
  - Added [`tests/test_integration.py`](tests/test_integration.py) covering end-to-end Python pipeline runs, grouped analysis, quality rescaling, filtering, and bin plotting.
  - Added [`tests/test_unit_models.py`](tests/test_unit_models.py) verifying statistical model properties, vectorized PMF fitting, Likelihood Ratio tests, RMSE, and GLM accuracy prediction.
  - Added [`tests/test_unit_parsing.py`](tests/test_unit_parsing.py) verifying CIGAR parsing, forward/reverse substitution tallying, soft-clipping/indels, and Phred quality score rescaling.

### Performance & Memory Optimization
- **Zero-Allocation In-Place Read Parsing**:
  - Refactored `damage_al` and `al_to_damage.get_damage()` to use in-place list accumulators instead of creating temporary 5-element dictionaries and lists per alignment. This eliminated over 210,000 temporary list allocations on standard runs and reduced parsing time from 1.23s to **73 ms** (a **16x speedup**).
- **Single-Pass Alignment & Coverage Statistics**:
  - Replaced multiple BAM traversals per contig (`count_coverage()`, `count()`, and `al_to_damage()`) with a single-pass reader in `al_to_damage.get_damage()` that calculates quality-filtered coverage sum, read count, and damage substitutions simultaneously (~11x faster per contig).
- **Persistent Worker BAM Handles**:
  - Multiprocessing worker pools now initialize a persistent `pysam.AlignmentFile` handle once on startup (`damage.init_worker_handle`), eliminating redundant file opening, header parsing, and index loading across thousands of contigs.
- **Pure Vectorized NumPy Fitting**:
  - Replaced `np.vectorize` in `damage_model.fit` with a pure NumPy formula `((1.0 - p) ** x) * p`, delivering a 3x speedup in optimizer iterations.
- **Fast Alignment & CIGAR Parsing**:
  - Removed per-position string reversals (`[::-1]`) inside inner alignment loops and replaced string concatenations with join/slice operations (2x faster parsing).
- **Zero-Overhead Single Process Execution**:
  - Added direct sequential execution fast-path when `process=1`, avoiding `multiprocessing.Pool` spawn latency on macOS/Python 3.12 (reducing runtimes from ~36s to ~0.2s).
- **Optimized Multiprocessing Context & Dynamic Chunking**:
  - Implemented `get_mp_context()` with dynamic chunk sizing (`chunksize = max(1, len(refs) // (proc * 4))`), enabling linear scaling up to **844 contigs/second** across 8 cores.
- **IPC Payload Minimization**:
  - `read_dict` serialization across Inter-Process Communication (IPC) channels is now only executed when BAM quality score rescaling is enabled (`rescale=True`).
- **Explicit Figure Memory Management**:
  - Added explicit `plt.close(fig)` calls in `plot.py` to prevent Matplotlib figure retention across multi-plot batches.

### Changed / Consolidated
- **Module Architecture Consolidation**:
  - Consolidated optimization routine (`optim`) and Likelihood Ratio test (`LR`) directly into [`pydamage/model_fit.py`](pydamage/model_fit.py).
  - Consolidated GLM accuracy modeling and feature preparation into [`pydamage/models.py`](pydamage/models.py).
- **Lazy CLI Subcommand Imports**:
  - Refactored [`pydamage/cli.py`](pydamage/cli.py) to lazy-import heavy dependencies (such as `matplotlib`) on demand, ensuring instantaneous execution of `pydamage cite` and `pydamage --help`.
- **Clean Package Exports**:
  - Removed side-effect `logging.basicConfig()` invocation from `pydamage/__init__.py`.

### Removed
- **Unused Single-Function Modules**:
  - Removed `pydamage/vuong.py` (unused `vuong_closeness` function).
  - Removed `pydamage/gc_content.py` (unused `fastaFile` class).
  - Removed redundant shim files `pydamage/accuracy_model.py`, `pydamage/likelihood_ratio.py`, and `pydamage/optim.py`.
- **Dead Helper Functions**:
  - Removed `get_ref_stats_pysam` and `get_ref_stats_multi` from `pydamage/contig_stats.py`.
  - Removed `pandas_group_processing`, `sort_dict_by_keys`, and `create_damage_dict` from `pydamage/utils.py`.
