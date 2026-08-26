import os
import pytest
from click.testing import CliRunner
from pydamage.cli import cli
from pydamage import __version__


@pytest.fixture
def runner():
    return CliRunner()


def test_cli_help(runner):
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "PyDamage: Damage parameter estimation for ancient DNA" in result.output
    assert "analyze" in result.output
    assert "filter" in result.output
    assert "binplot" in result.output
    assert "cite" in result.output


def test_cli_version(runner):
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


def test_cli_cite(runner):
    result = runner.invoke(cli, ["cite"])
    assert result.exit_code == 0
    assert "@article{borry_pydamage_2021" in result.output


def test_cli_analyze_basic(runner, tmp_path):
    outdir = str(tmp_path / "out_basic")
    result = runner.invoke(
        cli,
        ["-o", outdir, "analyze", "-f", "-w", "15", "-p", "1", "tests/data/aligned.bam"],
    )
    assert result.exit_code == 0
    csv_file = os.path.join(outdir, "pydamage_results.csv")
    assert os.path.exists(csv_file)


def test_cli_analyze_grouped(runner, tmp_path):
    outdir = str(tmp_path / "out_grouped")
    result = runner.invoke(
        cli,
        ["-o", outdir, "analyze", "-f", "-g", "-p", "1", "tests/data/aligned.bam"],
    )
    assert result.exit_code == 0
    csv_file = os.path.join(outdir, "pydamage_results.csv")
    assert os.path.exists(csv_file)


def test_cli_analyze_no_ga(runner, tmp_path):
    outdir = str(tmp_path / "out_noga")
    result = runner.invoke(
        cli,
        ["-o", outdir, "analyze", "-f", "-n", "-p", "1", "tests/data/aligned.bam"],
    )
    assert result.exit_code == 0
    csv_file = os.path.join(outdir, "pydamage_results.csv")
    assert os.path.exists(csv_file)


def test_cli_analyze_mutually_exclusive_options(runner, tmp_path):
    outdir = str(tmp_path / "out_excl")
    result = runner.invoke(
        cli,
        [
            "-o",
            outdir,
            "analyze",
            "-f",
            "-r",
            "-s",
            "0.5",
            "tests/data/aligned.bam",
        ],
    )
    assert result.exit_code != 0
    assert "mutually exclusive" in result.output.lower() or result.exit_code == 2


def test_cli_filter(runner, tmp_path):
    outdir = str(tmp_path / "out_filter")
    csv_file = "tests/data/pydamage_results_sequence.csv"
    result = runner.invoke(
        cli,
        ["-o", outdir, "-t", "0.6", "filter", csv_file],
    )
    assert result.exit_code == 0
    filtered_csv = os.path.join(outdir, "pydamage_filtered_results.csv")
    assert os.path.exists(filtered_csv)


def test_cli_binplot(runner, tmp_path):
    outdir = str(tmp_path / "out_binplot")
    os.makedirs(outdir, exist_ok=True)
    csv_file = "tests/data/pydamage_results_sequence.csv"
    fasta_file = "tests/data/sequence.fa"
    result = runner.invoke(
        cli,
        ["-o", outdir, "binplot", csv_file, fasta_file],
    )
    assert result.exit_code == 0
    png_file = os.path.join(outdir, "sequence.png")
    assert os.path.exists(png_file)
