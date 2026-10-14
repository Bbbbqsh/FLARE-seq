#!/usr/bin/env python3
"""
Usage: python flareshapemap_process_data.py --reference-table ncRNA_information.tsv
       --rnaframework-dir tools/RNAFramework --rf-fold tools/RNAFramework-fold-latest/rf-fold

Process FLARE-SHAPE-MaP reads and fold all RNAs using all replicates in multi-experiment mode.
The reference table contains name and seq columns. External tools must be available on PATH.

Python: 3.9.18
RNAFramework: 2.9.6
Bowtie2: 2.5.1
cutadapt: 4.5
samtools: 1.23.1
R: 4.2.1
ViennaRNA: 2.7.0
"""

import argparse
import csv
import math
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from shlex import quote
from textwrap import dedent


def build_index_for_transcriptomic_reference(reference_table, output_dir):
    """Write and sort the target transcriptome, then build its Bowtie2 index."""
    directory = output_dir / "01_index_building"
    directory.mkdir(parents=True, exist_ok=True)
    fasta = directory / "1P.filtered.groupNcrna.fa"

    with reference_table.open(newline="") as table, fasta.open("w") as handle:
        for row in csv.DictReader(table, delimiter="\t"):
            handle.write(f">{row['name']}\n{row['seq']}\n")

    # Preserve the original second-field sort and whole-line tie-breaking under the C locale.
    sorted_fasta = directory / "reference.sorted.fa"
    sort_command = dedent(rf"""
        set -e -o pipefail
        awk 'BEGIN{{RS=">"}} NR>1 {{gsub("\n", "\t"); print ">"$0}}' {quote(str(fasta))} | \
        LC_ALL=C sort -t ' ' -k 2,2 | \
        awk '{{sub("\t", "\n"); gsub("\t", ""); print $0}}' > {quote(str(sorted_fasta))}
    """).lstrip()
    subprocess.run(["bash", "-c", sort_command], check=True)

    index_dir = directory / "bowtie2_index"
    index_dir.mkdir(exist_ok=True)
    index_prefix = index_dir / "reference_sorted"
    subprocess.run(["bowtie2-build", "--threads", "20", sorted_fasta, index_prefix], check=True)
    return sorted_fasta, index_prefix


def map_reads_to_reference_transcriptome(samples, index_prefix, output_dir, rnaframework):
    """Map the paired-end FASTQs produced by preprocess.py to the target transcriptome."""
    directory = output_dir / "02_read_mapping"
    directory.mkdir(parents=True, exist_ok=True)

    for sample in samples:
        sample_id = sample["sample_id"]
        fastqs = [sample["directory"] / f"{sample_id}_{mate}.cutadapt.collapse.cutBarcodes.fq.gz"
                  for mate in ("R1", "R2")]
        print(f"Mapping {sample_id}", flush=True)
        subprocess.run([
            rnaframework / "rf-map", "--output-dir", directory / sample_id,
            "--bowtie", shutil.which("bowtie2"), "--cutadapt", shutil.which("cutadapt"),
            "--samtools", shutil.which("samtools"), "--bowtie-index", index_prefix,
            "--bowtie2", "--processors", "1", "--working-threads", "30", "--overwrite", "--keep-logs",
            "--cutadapt-5quality", "20", "--cutadapt-quality-only", "--clipped", "--bowtie-norc", "--bowtie-all",
            ",".join(str(path) for path in fastqs),
        ], check=True)


def calculate_mutations_and_read_coverage(samples, reference_fasta, output_dir, rnaframework):
    """Count mutations and read coverage and write mutation maps for each sample."""
    directory = output_dir / "03_mutations_and_read_coverage"
    directory.mkdir(parents=True, exist_ok=True)

    for sample in samples:
        sample_id = sample["sample_id"]
        sample_dir = directory / sample_id
        bam = output_dir / "02_read_mapping" / sample_id / f"{sample_id}.bam"
        print(f"Counting mutations: {sample_id}", flush=True)

        # --fast avoids FASTA lookup failures after rf-count transforms reference identifiers.
        subprocess.run([
            rnaframework / "rf-count", "--output-dir", sample_dir, "--samtools", shutil.which("samtools"),
            "--fasta", reference_fasta, "--processors", "10", "--working-threads", "2", "--overwrite", "--fast",
            "--count-mutations", "--mutation-map", bam,
        ], check=True)

        # rf-count can return zero even when a child process fails; preserve the original output checks.
        error_path = sample_dir / "error.out"
        if error_path.is_file() and error_path.stat().st_size > 0:
            raise RuntimeError(f"rf-count reported errors or warnings: {error_path}")
        for suffix in ("rc", "mm"):
            result = sample_dir / f"{sample_id}.{suffix}"
            if not result.is_file() or result.stat().st_size == 0:
                raise RuntimeError(f"rf-count output missing or empty: {result}")


def normalize_reactivities(samples, output_dir, rnaframework):
    """Pair DMSO and NAIN3 by cell line, end type and replicate, then normalize."""
    directory = output_dir / "04_normalized_reactivities"
    directory.mkdir(parents=True, exist_ok=True)
    pairs = defaultdict(dict)

    for sample in samples:
        source = f"{sample['cell_line']}-{sample['end_type']}"
        pairs[(source, sample["rep"])][sample["treatment"]] = sample["sample_id"]

    groups = defaultdict(list)
    for (source, replicate), treatments in sorted(pairs.items()):
        untreated_id, treated_id = treatments["DMSO"], treatments["NAIN3"]

        # Keep the sample-ID links used by the original normalization step.
        for sample_id in (untreated_id, treated_id):
            rc = output_dir / "03_mutations_and_read_coverage" / sample_id / f"{sample_id}.rc"
            link = directory / f"{sample_id}.rc"
            link.unlink(missing_ok=True)
            link.symlink_to(rc)

        group_name = f"{source}_{replicate}"
        norm_dir = directory / group_name
        print(f"Normalizing {group_name}", flush=True)
        subprocess.run([
            rnaframework / "rf-norm", "--untreated", directory / f"{untreated_id}.rc",
            "--treated", directory / f"{treated_id}.rc", "--output-dir", norm_dir,
            "--processors", "40", "--overwrite", "--scoring-method", "3", "--norm-method", "2",
        ], check=True)
        groups[source].append(norm_dir)

    groups["all"] = [path for paths in groups.values() for path in paths]
    return dict(groups)


def combine_reactivities(groups, output_dir, rnaframework):
    """Combine normalized replicates for each cell-line source and for all sources."""
    directory = output_dir / "05_combined_reactivities"
    directory.mkdir(parents=True, exist_ok=True)

    for group_name, norm_dirs in groups.items():
        print(f"Combining reactivities: {group_name}", flush=True)
        subprocess.run([
            rnaframework / "rf-combine", "--processors", "30", "--output-dir", directory / group_name,
            "--overwrite", *norm_dirs,
        ], check=True)


def fold_structure(norm_dirs, output_dir, rf_fold):
    """Fold all RNAs using the separate normalized replicates in default multi-experiment mode."""
    directory = output_dir / "06_folded_structures"
    directory.mkdir(parents=True, exist_ok=True)

    # rf-fold groups XMLs by RNA name across directories; use every RNA and every normalized replicate.
    subprocess.run([
        rf_fold, "--output-dir", directory / "all_multi_experiment", "--R-path", shutil.which("R"),
        "--vienna-rnasubopt", shutil.which("RNAsubopt"), "--vienna-rnafold", shutil.which("RNAfold"),
        "--processors", "40", "--overwrite", "--slope", "2.2", "--intercept", "-0.6", "--img", *norm_dirs,
    ], check=True)


def modify_reactivities_files(groups, output_dir):
    """Convert combined XML reactivities into position, nucleotide and reactivity columns."""
    directory = output_dir / "05_combined_reactivities"

    for group_name in groups:
        shape_dir = directory / f"{group_name}_shape"
        shape_dir.mkdir(parents=True, exist_ok=True)

        for xml_path in sorted((directory / group_name).glob("*.xml")):
            root = ET.parse(xml_path).getroot()
            shape_path = shape_dir / f"{xml_path.stem}.shape"

            for transcript in root.findall("transcript"):
                sequence = re.sub("[^A-Za-z]", "", transcript.find("sequence").text).upper().replace("T", "U")
                values = transcript.find("reactivity").text.split(",")

                # Preserve the original conversion: invalid/nonfinite values become zero; positions start at one.
                with shape_path.open("w") as handle:
                    for position, value in enumerate(values, 1):
                        try:
                            reactivity = float(value.strip())
                        except ValueError:
                            reactivity = 0.0
                        if not math.isfinite(reactivity):
                            reactivity = 0.0
                        handle.write(f"{position}\t{sequence[position - 1]}\t{reactivity}\n")


def main():
    """Run the seven FLARE-SHAPE-MaP processing steps sequentially."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample-info", type=Path, default=Path(__file__).with_name("sample_info.tsv"))
    parser.add_argument("--reference-table", type=Path, required=True, help="Transcript table: name and seq columns")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).with_name("shapemap_results"))
    parser.add_argument("--rnaframework-dir", type=Path, required=True, help="RNAFramework for map/count/norm/combine")
    parser.add_argument("--rf-fold", type=Path, required=True, help="Path to the newer multi-experiment rf-fold")
    args = parser.parse_args()

    sample_info = args.sample_info.resolve()
    output_dir, rnaframework = args.output_dir.resolve(), args.rnaframework_dir.resolve()
    with sample_info.open(newline="") as handle:
        samples = [row for row in csv.DictReader(handle, delimiter="\t") if row["assay"] == "FLARE-SHAPE-MaP"]

    for sample in samples:
        sample["directory"] = sample_info.parent / sample["output_dir"]

    reference_fasta, index_prefix = build_index_for_transcriptomic_reference(args.reference_table.resolve(), output_dir)
    map_reads_to_reference_transcriptome(samples, index_prefix, output_dir, rnaframework)
    calculate_mutations_and_read_coverage(samples, reference_fasta, output_dir, rnaframework)
    groups = normalize_reactivities(samples, output_dir, rnaframework)
    combine_reactivities(groups, output_dir, rnaframework)
    fold_structure(groups["all"], output_dir, args.rf_fold.resolve())
    modify_reactivities_files(groups, output_dir)


if __name__ == "__main__":
    main()
