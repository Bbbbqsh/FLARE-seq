#!/usr/bin/env python3
"""
Usage: python flareseq_process_data.py --genome-index dm6 indexes/dm6
       --genome-fasta dm6 genome/dm6.fa --chrom-sizes dm6 genome/dm6.chrom.sizes
       --genome-gtf dm6 genome/Drosophila_melanogaster.BDGP6.28.101.chr.gtf
       
Repeat the genome options for hg38 when processing the full sample table. FASTA requires a .fai index.
STAR indexes are built before alignment.

Align preprocessed FLARE-seq paired-end FASTQ files with STAR and generate six bigWig tracks.

Python: 3.9.18
STAR: 2.7.10a
bamToWig: (in-house script)
wigToBigWig: UCSC utility (v4)
samtools: 1.23.1
"""

import argparse
import csv
import subprocess
from pathlib import Path


def build_star_index(genome, genome_fasta, genome_gtf, genome_index, threads=20):
    """Build STAR indexes using the original dm6/hg38 genomeGenerate parameters."""
    genome_index = Path(genome_index)
    genome_index.mkdir(parents=True, exist_ok=True)
    # Original annotations: dm6 BDGP6.28.101.chr.gtf; hg38 GENCODE v39.
    command = [
        "STAR", "--runMode", "genomeGenerate", "--runThreadN", str(threads), "--genomeDir", genome_index,
        "--genomeFastaFiles", genome_fasta, "--sjdbGTFfile", genome_gtf, "--sjdbOverhang", "100",
    ]
    if genome == "dm6":
        command.extend(["--genomeSAindexNbases", "12"])
    # hg38 uses STAR's default genomeSAindexNbases=14; keep each index's log in its own directory.
    command.extend(["--outFileNamePrefix", str(genome_index) + "/"])
    subprocess.run(command, check=True)


def align_flare_seq(r1, r2, genome_index, output_dir, prefix, threads=16):
    """Run the original STAR alignment parameters, then sort and index the BAM."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    star_prefix = output_dir / f"{prefix}.EndToEnd.STAR."

    # Preserve the original EndToEnd settings, including multimapping and splice-junction parameters.
    subprocess.run([
        "STAR", "--runThreadN", str(threads), "--genomeLoad", "NoSharedMemory",
        "--limitBAMsortRAM", "60000000000", "--alignEndsType", "EndToEnd", "--outFilterType", "Normal",
        "--outFilterMultimapScoreRange", "0", "--outFilterMultimapNmax", "20",
        "--outFilterMismatchNmax", "15", "--outFilterMismatchNoverLmax", "0.1",
        "--outFilterScoreMin", "0", "--outFilterScoreMinOverLread", "0", "--outFilterMatchNmin", "18",
        "--outFilterMatchNminOverLread", "0.8", "--alignIntronMin", "5", "--seedSearchStartLmax", "15",
        "--seedSearchStartLmaxOverLread", "1", "--seedSearchLmax", "0",
        "--alignTranscriptsPerReadNmax", "20000", "--alignWindowsPerReadNmax", "20000",
        "--seedMultimapNmax", "20000", "--seedPerReadNmax", "1000", "--seedPerWindowNmax", "100",
        "--seedNoneLociPerWindow", "20", "--outSAMtype", "BAM", "Unsorted", "--outSAMmode", "Full",
        "--outSAMattributes", "All", "--outSAMunmapped", "None", "--outSAMorder", "Paired",
        "--outSAMprimaryFlag", "AllBestScore", "--outSAMreadID", "Standard", "--outReadsUnmapped", "Fastx",
        "--scoreGapNoncan", "-8", "--scoreGapATAC", "-8", "--scoreGapGCAG", "-4",
        "--alignSJoverhangMin", "15", "--alignSJDBoverhangMin", "5",
        "--alignEndsProtrude", "150", "ConcordantPair", "--scoreGenomicLengthLog2scale", "-1",
        "--readFilesCommand", "zcat", "--genomeDir", genome_index, "--readFilesIn", r1, r2,
        "--outFileNamePrefix", star_prefix,
    ], check=True)

    # Retain a sorted, indexed BAM for downstream analysis.
    bam = Path(f"{star_prefix}Aligned.out.bam")
    sorted_bam = Path(f"{star_prefix}Aligned.out.sorted.bam")
    subprocess.run(["samtools", "sort", "-@", str(threads), "-o", sorted_bam, bam], check=True)
    subprocess.run(["samtools", "index", "-@", str(threads), sorted_bam], check=True)
    return bam


def generate_flare_seq_bigwigs(bam, output_dir, prefix, genome_fasta, chrom_sizes):
    """Generate six strand-specific RPM bigWigs using bamToWig and wigToBigWig."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Preserve the original paired-end, barcode and mispriming-filter parameters.
    subprocess.run([
        "bamToWig", "-s", "--rpm", "--primer", "TGGAA", "-L", "-1000", "-P", "-u", "4",
        "--fa", genome_fasta, "--fai", f"{genome_fasta}.fai", "-p", prefix, "-o", output_dir, "--bam", bam,
    ], check=True)
    for strand in ("plus", "minus"):
        for signal in ("coverage", "start", "end"):
            wig = output_dir / f"{prefix}.rpm.{strand}.{signal}.wg"
            subprocess.run(["wigToBigWig", wig, chrom_sizes, wig.with_suffix(".bw")], check=True)


def main():
    """Read the shared sample table and process FLARE-seq samples sequentially."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample-info", type=Path, default=Path(__file__).with_name("sample_info.tsv"))
    parser.add_argument("--threads", type=int, default=16, help="Number of threads per sample (default: 16)")
    parser.add_argument("--genome-index", nargs=2, action="append", required=True, metavar=("GENOME", "PATH"),
                        help="STAR index for each genome in the FLAREseq sample rows; repeat for multiple genomes")
    parser.add_argument("--genome-fasta", nargs=2, action="append", required=True, metavar=("GENOME", "PATH"),
                        help="Reference FASTA with an existing .fai index; repeat for multiple genomes")
    parser.add_argument("--chrom-sizes", nargs=2, action="append", required=True, metavar=("GENOME", "PATH"),
                        help="Two-column chromosome sizes file; repeat for multiple genomes")
    parser.add_argument("--genome-gtf", nargs=2, action="append", required=True, metavar=("GENOME", "PATH"),
                        help="GTF for STAR indexing: dm6 BDGP6.28.101.chr.gtf or hg38 GENCODE v39")
    args = parser.parse_args()
    sample_info = args.sample_info.resolve()
    genome_indexes = {genome: Path(path).resolve() for genome, path in args.genome_index}

    genome_fastas = {genome: Path(path).resolve() for genome, path in args.genome_fasta}
    chrom_sizes = {genome: Path(path).resolve() for genome, path in args.chrom_sizes}

    with sample_info.open(newline="") as handle:
        samples = [row for row in csv.DictReader(handle, delimiter="\t") if row["assay"] == "FLAREseq"]
    genomes = {sample["genome"] for sample in samples}
    genome_gtfs = {genome: Path(path).resolve() for genome, path in args.genome_gtf}
    for genome in sorted(genomes):
        print(f"Building STAR index: {genome}", flush=True)
        build_star_index(genome, genome_fastas[genome], genome_gtfs[genome], genome_indexes[genome])

    for sample in samples:
        output_dir = sample_info.parent / sample["output_dir"]
        sample_id = sample["sample_id"]
        r1 = output_dir / f"{sample_id}_R1.cutadapt.collapse.cutBarcodes.fq.gz"
        r2 = output_dir / f"{sample_id}_R2.cutadapt.collapse.cutBarcodes.fq.gz"
        prefix = f"{sample['genome']}_{sample_id}_{sample['sample_name']}"
        print(f"Processing FLAREseq: {sample_id}", flush=True)
        bam = align_flare_seq(r1, r2, genome_indexes[sample["genome"]], output_dir, prefix, args.threads)
        generate_flare_seq_bigwigs(bam, output_dir, prefix, genome_fastas[sample["genome"]],
                                   chrom_sizes[sample["genome"]])


if __name__ == "__main__":
    main()
