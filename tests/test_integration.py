import os
import pytest
import pysam
import pandas as pd
from pydamage.main import pydamage_analyze
from pydamage.filter import apply_filter, define_threshold, filter_pydamage_results
from pydamage.plot import bin_plot


def test_integration_standard_analysis(tmp_path):
    outdir = str(tmp_path / "out_std")
    bam = "tests/data/aligned.bam"
    
    df = pydamage_analyze(
        bam=bam,
        outdir=outdir,
        wlen=20,
        minlen=0,
        process=1,
        force=True,
        plot=True,
        verbose=False,
    )
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 31
    assert "predicted_accuracy" in df.columns
    assert "null_model_p0" in df.columns
    assert "damage_model_p" in df.columns
    assert "pvalue" in df.columns
    assert "qvalue" in df.columns
    assert "RMSE" in df.columns
    assert "CtoT-0" in df.columns
    assert "GtoA-0" in df.columns
    assert "CtoT-19" in df.columns
    assert "GtoA-19" in df.columns
    assert "CtoT-20" not in df.columns
    
    # Check output files
    assert os.path.exists(os.path.join(outdir, "pydamage_results.csv"))
    assert os.path.exists(os.path.join(outdir, "plots"))
    # Check that plots were generated
    plot_files = os.listdir(os.path.join(outdir, "plots"))
    assert len(plot_files) == 31


def test_integration_grouped_analysis(tmp_path):
    outdir = str(tmp_path / "out_grp")
    bam = "tests/data/aligned.bam"
    
    df = pydamage_analyze(
        bam=bam,
        outdir=outdir,
        wlen=15,
        group=True,
        process=1,
        force=True,
    )
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1
    assert df.index[0] == "reference"
    assert "predicted_accuracy" in df.columns
    assert df.loc["reference", "nb_reads_aligned"] > 0
    assert os.path.exists(os.path.join(outdir, "pydamage_results.csv"))


def test_integration_rescaling(tmp_path):
    outdir = str(tmp_path / "out_rescale")
    bam = "tests/data/aligned.bam"
    
    df = pydamage_analyze(
        bam=bam,
        outdir=outdir,
        wlen=20,
        rescale=True,
        process=1,
        force=True,
    )
    
    rescaled_bam = os.path.join(outdir, "pydamage_rescaled.bam")
    assert os.path.exists(rescaled_bam)
    
    # Verify the output BAM is valid and readable without needing a .bai index
    with pysam.AlignmentFile(rescaled_bam, "rb") as af:
        reads = list(af.fetch(until_eof=True))
        assert len(reads) > 0


def test_integration_filtering_workflow(tmp_path):
    outdir = str(tmp_path / "out_filt")
    csv_file = "tests/data/pydamage_results_sequence.csv"
    
    # Test apply_filter
    filt_df = apply_filter(csv=csv_file, threshold=0.6, outdir=outdir, alpha=0.05)
    assert isinstance(filt_df, pd.DataFrame)
    assert os.path.exists(os.path.join(outdir, "pydamage_filtered_results.csv"))
    
    # Test define_threshold
    raw_df = pd.read_csv(csv_file)
    kneedle_thresh = define_threshold(raw_df)
    assert kneedle_thresh is not None
    assert 0.0 <= kneedle_thresh <= 1.0


def test_integration_binplot():
    csv = "tests/data/pydamage_results_sequence.csv"
    fasta = "tests/data/sequence.fa"
    ct_mean, ct_std, ga_mean, ga_std = bin_plot(csv, fasta, None, write_fig=False)
    
    assert len(ct_mean) > 0
    assert len(ct_std) > 0
    assert len(ga_mean) > 0
    assert len(ga_std) > 0
