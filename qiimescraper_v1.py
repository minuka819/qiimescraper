#!/usr/bin/env python3
"""
QIIME Scraper
=============

A lightweight command-line wrapper for running a paired-end QIIME2 amplicon
workflow from raw FASTQ files to taxonomy barplots.

Main steps:
1. Generate QIIME2 paired-end manifest
2. Import FASTQs as a QIIME2 artifact
3. Trim primers with cutadapt
4. Denoise with DADA2
5. Assign taxonomy with a trained classifier
6. Generate taxa barplot

This script is intentionally marker-agnostic. It can be used for 16S, ITS,
COI, or other amplicons as long as appropriate primers and a compatible
QIIME2 classifier are supplied.
"""

import argparse
import os
import shlex
import subprocess
from pathlib import Path
from typing import List, Tuple


# -----------------------------------------------------------------------------
# Command runner
# -----------------------------------------------------------------------------
def run_command(cmd: List[str], log_file: Path | None = None) -> None:
    """Run a shell command safely from a list of arguments."""
    printable = " ".join(shlex.quote(str(x)) for x in cmd)
    print(f"\n[Running]\n{printable}\n")

    if log_file:
        with open(log_file, "w") as lf:
            process = subprocess.run(cmd, stdout=lf, stderr=lf, text=True)
    else:
        process = subprocess.run(cmd, text=True)

    if process.returncode != 0:
        raise subprocess.CalledProcessError(process.returncode, printable)


# -----------------------------------------------------------------------------
# Primer handling
# -----------------------------------------------------------------------------
def parse_comma_list(value: str | None) -> List[str]:
    """Parse comma-separated primer strings from the CLI."""
    if not value:
        return []
    return [item.strip() for item in value.split(",") if item.strip()]


def load_primers_from_tsv(path: Path, primer_set: str) -> Tuple[List[str], List[str]]:
    """
    Load primers from a TSV file.

    Expected columns:
        primer_set    direction    sequence

    direction must be one of:
        forward, reverse

    Example:
        primer_set  direction  sequence
        ITS         forward    CTTGGTCATTTAGAGGAAGTAA
        ITS         reverse    GCTGCGTTCTTCATCGATGC
    """
    forward: List[str] = []
    reverse: List[str] = []

    with open(path, "r") as handle:
        header = handle.readline().strip().split("\t")
        required = {"primer_set", "direction", "sequence"}
        missing = required.difference(header)
        if missing:
            raise ValueError(
                f"Primer file is missing required column(s): {', '.join(sorted(missing))}"
            )

        col = {name: idx for idx, name in enumerate(header)}

        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            row_set = fields[col["primer_set"]]
            direction = fields[col["direction"]].lower()
            sequence = fields[col["sequence"]]

            if row_set != primer_set:
                continue

            if direction == "forward":
                forward.append(sequence)
            elif direction == "reverse":
                reverse.append(sequence)
            else:
                raise ValueError(
                    f"Invalid direction '{direction}' in primer file. Use 'forward' or 'reverse'."
                )

    if not forward or not reverse:
        raise ValueError(
            f"No complete primer set found for '{primer_set}' in {path}. "
            "At least one forward and one reverse primer are required."
        )

    return forward, reverse


def get_primers(args: argparse.Namespace) -> Tuple[List[str], List[str]]:
    """Resolve primers from either CLI arguments or a primer TSV file."""
    if args.primers_file:
        if not args.primer_set:
            raise ValueError("--primer-set is required when using --primers-file")
        return load_primers_from_tsv(Path(args.primers_file), args.primer_set)

    forward = parse_comma_list(args.forward_primers)
    reverse = parse_comma_list(args.reverse_primers)

    if not forward or not reverse:
        raise ValueError(
            "You must provide primers using either:\n"
            "  1) --forward-primers and --reverse-primers\n"
            "  OR\n"
            "  2) --primers-file and --primer-set"
        )

    return forward, reverse


# -----------------------------------------------------------------------------
# Manifest generation
# -----------------------------------------------------------------------------
def infer_sample_id(r1_filename: str, r1_suffix: str) -> str:
    """Infer sample ID by removing the R1 suffix from the filename."""
    if not r1_filename.endswith(r1_suffix):
        raise ValueError(f"R1 file does not end with expected suffix: {r1_filename}")
    return r1_filename[: -len(r1_suffix)]


def generate_manifest(
    fastq_dir: Path,
    manifest_path: Path,
    r1_suffix: str = "_R1_001.fastq.gz",
    r2_suffix: str = "_R2_001.fastq.gz",
    sample_contains: str | None = None,
) -> None:
    """
    Generate a QIIME2 paired-end manifest file.

    Parameters
    ----------
    fastq_dir:
        Directory containing paired-end FASTQ files.
    manifest_path:
        Output manifest TSV path.
    r1_suffix:
        Filename suffix identifying forward reads.
    r2_suffix:
        Filename suffix identifying reverse reads.
    sample_contains:
        Optional substring filter for FASTQ filenames.
    """
    fastq_dir = fastq_dir.resolve()
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    r1_files = sorted(
        f for f in os.listdir(fastq_dir)
        if f.endswith(r1_suffix) and (sample_contains is None or sample_contains in f)
    )

    if not r1_files:
        raise FileNotFoundError(
            f"No R1 files found in {fastq_dir} with suffix '{r1_suffix}'."
        )

    with open(manifest_path, "w") as handle:
        handle.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")

        for r1 in r1_files:
            sample_id = infer_sample_id(r1, r1_suffix)
            r2 = r1.replace(r1_suffix, r2_suffix)

            r1_path = fastq_dir / r1
            r2_path = fastq_dir / r2

            if not r2_path.exists():
                raise FileNotFoundError(f"Missing matching R2 file for {r1}: {r2_path}")

            handle.write(f"{sample_id}\t{r1_path.resolve()}\t{r2_path.resolve()}\n")

    print(f"[Manifest created] {manifest_path}")


# -----------------------------------------------------------------------------
# QIIME2 pipeline steps
# -----------------------------------------------------------------------------
def step_import(manifest: Path, outdir: Path) -> None:
    run_command([
        "qiime", "tools", "import",
        "--type", "SampleData[PairedEndSequencesWithQuality]",
        "--input-path", str(manifest),
        "--output-path", str(outdir / "demux.qza"),
        "--input-format", "PairedEndFastqManifestPhred33V2",
    ])

    run_command([
        "qiime", "demux", "summarize",
        "--i-data", str(outdir / "demux.qza"),
        "--o-visualization", str(outdir / "demux.qzv"),
    ])


def step_cutadapt(
    outdir: Path,
    forward_primers: List[str],
    reverse_primers: List[str],
    cores: int,
    discard_untrimmed: bool,
) -> None:
    cmd = [
        "qiime", "cutadapt", "trim-paired",
        "--i-demultiplexed-sequences", str(outdir / "demux.qza"),
        "--verbose",
        "--p-cores", str(cores),
        "--o-trimmed-sequences", str(outdir / "demux_trimmed.qza"),
    ]

    for primer in forward_primers:
        cmd.extend(["--p-front-f", primer])

    for primer in reverse_primers:
        cmd.extend(["--p-front-r", primer])

    if discard_untrimmed:
        cmd.append("--p-discard-untrimmed")

    run_command(cmd, log_file=outdir / "cutadapt_verbose.log")


def step_dada2(
    outdir: Path,
    cores: int,
    trunc_len_f: int,
    trunc_len_r: int,
    trim_left_f: int,
    trim_left_r: int,
) -> None:
    run_command([
        "qiime", "dada2", "denoise-paired",
        "--i-demultiplexed-seqs", str(outdir / "demux_trimmed.qza"),
        "--p-trunc-len-f", str(trunc_len_f),
        "--p-trunc-len-r", str(trunc_len_r),
        "--p-trim-left-f", str(trim_left_f),
        "--p-trim-left-r", str(trim_left_r),
        "--p-n-threads", str(cores),
        "--o-table", str(outdir / "table.qza"),
        "--o-representative-sequences", str(outdir / "repseqs.qza"),
        "--o-denoising-stats", str(outdir / "dada2stats.qza"),
    ])

    run_command([
        "qiime", "metadata", "tabulate",
        "--m-input-file", str(outdir / "dada2stats.qza"),
        "--o-visualization", str(outdir / "dada2stats.qzv"),
    ])


def step_taxonomy(outdir: Path, classifier: Path, cores: int) -> None:
    if not classifier:
        raise ValueError("--classifier is required for taxonomy classification")

    run_command([
        "qiime", "feature-classifier", "classify-sklearn",
        "--i-classifier", str(classifier),
        "--i-reads", str(outdir / "repseqs.qza"),
        "--p-n-jobs", str(cores),
        "--o-classification", str(outdir / "taxonomy.qza"),
    ])

    run_command([
        "qiime", "metadata", "tabulate",
        "--m-input-file", str(outdir / "taxonomy.qza"),
        "--o-visualization", str(outdir / "taxonomy.qzv"),
    ])


def step_barplot(outdir: Path, metadata: Path) -> None:
    if not metadata:
        raise ValueError("--metadata is required for taxa barplot generation")

    run_command([
        "qiime", "taxa", "barplot",
        "--i-table", str(outdir / "table.qza"),
        "--i-taxonomy", str(outdir / "taxonomy.qza"),
        "--m-metadata-file", str(metadata),
        "--o-visualization", str(outdir / "taxa-barplot.qzv"),
    ])


# -----------------------------------------------------------------------------
# Workflow controller
# -----------------------------------------------------------------------------
def run_full_pipeline(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    manifest = Path(args.manifest) if args.manifest else outdir / "manifest.tsv"

    outdir.mkdir(parents=True, exist_ok=True)

    forward_primers, reverse_primers = get_primers(args)

    print("\n=== QIIME Scraper: full pipeline ===")
    print(f"FASTQ directory: {args.fastq_dir}")
    print(f"Output directory: {outdir}")
    print(f"Manifest: {manifest}")
    print(f"Forward primers: {', '.join(forward_primers)}")
    print(f"Reverse primers: {', '.join(reverse_primers)}")
    print(f"Discard untrimmed: {args.discard_untrimmed}")
    print(f"Cores: {args.cores}\n")

    print("Step 1/6: Generating manifest")
    generate_manifest(
        fastq_dir=Path(args.fastq_dir),
        manifest_path=manifest,
        r1_suffix=args.r1_suffix,
        r2_suffix=args.r2_suffix,
        sample_contains=args.sample_contains,
    )

    print("Step 2/6: Importing FASTQs")
    step_import(manifest, outdir)

    print("Step 3/6: Trimming primers with cutadapt")
    step_cutadapt(outdir, forward_primers, reverse_primers, args.cores, args.discard_untrimmed)

    print("Step 4/6: Denoising with DADA2")
    step_dada2(
        outdir=outdir,
        cores=args.cores,
        trunc_len_f=args.trunc_len_f,
        trunc_len_r=args.trunc_len_r,
        trim_left_f=args.trim_left_f,
        trim_left_r=args.trim_left_r,
    )

    print("Step 5/6: Assigning taxonomy")
    step_taxonomy(outdir, Path(args.classifier), args.cores)

    print("Step 6/6: Creating taxa barplot")
    step_barplot(outdir, Path(args.metadata))

    print("\n=== Pipeline completed successfully ===\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run a configurable paired-end QIIME2 amplicon workflow."
    )

    parser.add_argument(
        "--step",
        choices=["manifest", "import", "cutadapt", "dada2", "taxonomy", "barplot"],
        help="Run one pipeline step only. If omitted with --full, the full pipeline is run.",
    )
    parser.add_argument("--full", action="store_true", help="Run the complete pipeline.")

    parser.add_argument("--fastq-dir", default="raw_fastqs", help="Directory containing FASTQ files.")
    parser.add_argument("--outdir", default="outputs/qiimescraper_run", help="Output directory.")
    parser.add_argument("--manifest", help="Optional manifest output/input path.")
    parser.add_argument("--metadata", help="Metadata TSV file for taxa barplot.")
    parser.add_argument("--classifier", help="QIIME2 trained classifier .qza file.")

    parser.add_argument(
        "--forward-primers",
        help="Comma-separated forward primer sequence(s). Example: ACTG,GTCA",
    )
    parser.add_argument(
        "--reverse-primers",
        help="Comma-separated reverse primer sequence(s). Example: TGAC,CAGT",
    )
    parser.add_argument(
        "--primers-file",
        help="Optional TSV primer file with columns: primer_set, direction, sequence.",
    )
    parser.add_argument(
        "--primer-set",
        help="Primer set name to use from --primers-file.",
    )

    parser.add_argument("--cores", type=int, default=8, help="Number of CPU cores/threads.")
    parser.add_argument("--discard-untrimmed", action="store_true", help="Discard reads without detected primers.")

    parser.add_argument("--r1-suffix", default="_R1_001.fastq.gz", help="Forward read filename suffix.")
    parser.add_argument("--r2-suffix", default="_R2_001.fastq.gz", help="Reverse read filename suffix.")
    parser.add_argument(
        "--sample-contains",
        help="Optional substring filter for sample FASTQ names. Useful for processing one marker at a time.",
    )

    parser.add_argument("--trunc-len-f", type=int, default=0, help="DADA2 forward truncation length.")
    parser.add_argument("--trunc-len-r", type=int, default=0, help="DADA2 reverse truncation length.")
    parser.add_argument("--trim-left-f", type=int, default=0, help="DADA2 trim-left forward bases.")
    parser.add_argument("--trim-left-r", type=int, default=0, help="DADA2 trim-left reverse bases.")

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest = Path(args.manifest) if args.manifest else outdir / "manifest.tsv"

    if args.full:
        run_full_pipeline(args)
        return

    if not args.step:
        parser.error("Choose either --full or --step")

    if args.step == "manifest":
        generate_manifest(
            fastq_dir=Path(args.fastq_dir),
            manifest_path=manifest,
            r1_suffix=args.r1_suffix,
            r2_suffix=args.r2_suffix,
            sample_contains=args.sample_contains,
        )
        return

    if args.step == "import":
        step_import(manifest, outdir)
        return

    if args.step == "cutadapt":
        forward_primers, reverse_primers = get_primers(args)
        step_cutadapt(outdir, forward_primers, reverse_primers, args.cores, args.discard_untrimmed)
        return

    if args.step == "dada2":
        step_dada2(
            outdir=outdir,
            cores=args.cores,
            trunc_len_f=args.trunc_len_f,
            trunc_len_r=args.trunc_len_r,
            trim_left_f=args.trim_left_f,
            trim_left_r=args.trim_left_r,
        )
        return

    if args.step == "taxonomy":
        step_taxonomy(outdir, Path(args.classifier), args.cores)
        return

    if args.step == "barplot":
        step_barplot(outdir, Path(args.metadata))
        return


if __name__ == "__main__":
    main()