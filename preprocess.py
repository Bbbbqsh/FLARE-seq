#!/usr/bin/env python3
"""
Usage: python preprocess.py --sample-info sample_info.tsv --threads 16

Paired-end preprocessing for FLARE-seq and FLARE-SHAPE-MaP.

Raw sequencing data can be downloaded from the Gene Expression Omnibus (GEO) database.

The default sample_info.tsv is located beside this script. Edit its r1, r2 and output_dir
columns to match your data locations before running. Relative paths are resolved against
the sample-info directory, not the working directory. Samples are processed sequentially;
threads controls cutadapt parallelism within each sample. Intermediate FASTQ files are retained.

Python: 3.9.18
cutadapt: 4.5
"""

import argparse
import csv
import subprocess
from pathlib import Path


def preprocess_flare_seq(r1, r2, sample_id, output_dir, threads=16):
    """Run cutadapt, collapse and cutBarcode for one FLARE-seq sample.

    r1 and r2 are raw paired-end FASTQ paths. sample_id is the output filename prefix.
    Returns the two preprocessed FASTQ paths.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    trimmed_r1 = output_dir / f"{sample_id}_R1.cutadapt.fq.gz"
    trimmed_r2 = output_dir / f"{sample_id}_R2.cutadapt.fq.gz"
    collapsed_r1 = output_dir / f"{sample_id}_R1.cutadapt.collapse.fq.gz"
    collapsed_r2 = output_dir / f"{sample_id}_R2.cutadapt.collapse.fq.gz"
    final_r1 = output_dir / f"{sample_id}_R1.cutadapt.collapse.cutBarcodes.fq.gz"
    final_r2 = output_dir / f"{sample_id}_R2.cutadapt.collapse.cutBarcodes.fq.gz"

    # 1. Remove the R1/R2 adapters; retain reads of at least 10 nt (barcodes still included).
    subprocess.run([
        "cutadapt", "--quiet", "-a", "TGGAATTCTCGGGTGCCAAGG", "-A", "GATCGTCGGACTGTAGAACTCTGAAC",
        "-m", "10", "-e", "0.15", "-j", str(threads), "-o", trimmed_r1, "-p", trimmed_r2, r1, r2,
    ], check=True)

    # 2. Collapse duplicate reads before removing barcodes, preserving barcode information for this step.
    subprocess.run([
        "collapseRead", "-f", sample_id, "-c", "-1", trimmed_r1, "-2", trimmed_r2, "-o", output_dir,
    ], check=True)
    # collapseRead uses its own output suffix; restore the cumulative processing suffix.
    (output_dir / f"{sample_id}_R1.collapse.fq.gz").rename(collapsed_r1)
    (output_dir / f"{sample_id}_R2.collapse.fq.gz").rename(collapsed_r2)

    # 3. cutBarcode: remove 4 nt from both ends of each mate, then require at least 10 nt.
    subprocess.run([
        "cutadapt", "--quiet", "-j", str(threads), "-u", "4", "-u", "-4", "-U", "4", "-U", "-4",
        "-m", "10", "-o", final_r1, "-p", final_r2, collapsed_r1, collapsed_r2,
    ], check=True)

    return final_r1, final_r2


def preprocess_flare_shape_map(r1, r2, sample_id, output_dir, threads=16):
    """Run cutadapt, collapse and cutBarcode for one FLARE-SHAPE-MaP sample.

    r1 and r2 are raw paired-end FASTQ paths. sample_id is the output filename prefix.
    Returns the two preprocessed FASTQ paths.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    trimmed_r1 = output_dir / f"{sample_id}_R1.cutadapt.fq.gz"
    trimmed_r2 = output_dir / f"{sample_id}_R2.cutadapt.fq.gz"
    collapsed_r1 = output_dir / f"{sample_id}_R1.cutadapt.collapse.fq.gz"
    collapsed_r2 = output_dir / f"{sample_id}_R2.cutadapt.collapse.fq.gz"
    final_r1 = output_dir / f"{sample_id}_R1.cutadapt.collapse.cutBarcodes.fq.gz"
    final_r2 = output_dir / f"{sample_id}_R2.cutadapt.collapse.cutBarcodes.fq.gz"

    # 1. Remove the R1/R2 adapters; retain reads of at least 10 nt (barcodes still included).
    subprocess.run([
        "cutadapt", "--quiet", "-a", "TGGAATTCTCGGGTGCCAAGG", "-A", "GATCGTCGGACTGTAGAACTCTGAAC",
        "-m", "10", "-e", "0.15", "-j", str(threads), "-o", trimmed_r1, "-p", trimmed_r2, r1, r2,
    ], check=True)

    # 2. Collapse duplicate reads before removing barcodes, preserving barcode information for this step.
    subprocess.run([
        "collapseRead", "-f", sample_id, "-c", "-1", trimmed_r1, "-2", trimmed_r2, "-o", output_dir,
    ], check=True)
    # collapseRead uses its own output suffix; restore the cumulative processing suffix.
    (output_dir / f"{sample_id}_R1.collapse.fq.gz").rename(collapsed_r1)
    (output_dir / f"{sample_id}_R2.collapse.fq.gz").rename(collapsed_r2)

    # 3. cutBarcode: remove 6 nt from both ends of each mate, then require at least 10 nt.
    subprocess.run([
        "cutadapt", "--quiet", "-j", str(threads), "-u", "6", "-u", "-6", "-U", "6", "-U", "-6",
        "-m", "10", "-o", final_r1, "-p", final_r2, collapsed_r1, collapsed_r2,
    ], check=True)

    return final_r1, final_r2


def main():
    """Read the sample table and run the appropriate workflow for each sample."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample-info", type=Path, default=Path(__file__).with_name("sample_info.tsv"))
    parser.add_argument("--threads", type=int, default=16, help="Number of cutadapt threads per sample (default: 16)")
    args = parser.parse_args()
    sample_info = args.sample_info.resolve()
    workflows = {"FLAREseq": preprocess_flare_seq, "FLARE-SHAPE-MaP": preprocess_flare_shape_map}

    with sample_info.open(newline="") as handle:
        for sample in csv.DictReader(handle, delimiter="\t"):
            # Resolve the portable paths relative to the table; absolute paths are also accepted.
            r1 = sample_info.parent / sample["r1"]
            r2 = sample_info.parent / sample["r2"]
            output_dir = sample_info.parent / sample["output_dir"]
            print(f"Processing {sample['assay']}: {sample['sample_id']}", flush=True)
            workflows[sample["assay"]](r1, r2, sample["sample_id"], output_dir, threads=args.threads)


if __name__ == "__main__":
    main()
