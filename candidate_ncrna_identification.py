#!/usr/bin/env python3
"""
Usage: python candidate_ncrna_identification.py --genome-fasta dm6 genome/dm6.fa
       --genome-fasta hg38 genome/hg38.fa
FASTA requires a .fai index.

Identify and filter ncRNA candidates from the BAM files produced by flareseq_process_data.py.

Python: 3.9.18
pandas: 2.0.3
PyDESeq2: 0.4.12
ncrnaScan: in-house tool
featureCounts: v2.1.1
bedtools: v2.31.1
"""

import argparse
import csv
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

END_TYPES = ("1P", "3P", "CAP", "OH")
CELL_LINES = {"dm6": ("Kc167", "S2"), "hg38": ("HepG2", "U87")}
BED_COLUMNS = ["#chrom", "start", "end", "name", "score", "strand"]

STATS_COLUMNS = ["supportSamples", "supportSampleNum", "groupScoreSum", "representativeScoreSum",
                 "representativeScoreBalance", "representativeLoc", "representativeReason", "groupIsoNum"]


def scan_sample(sample, genome_fasta):
    """Run ncrnaScan and label identical candidate sequences across genomic loci."""
    directory, prefix = sample["directory"], sample["prefix"]
    directory.mkdir(parents=True, exist_ok=True)
    base = directory / f"{prefix}.mispriming.ncrnaScan"
    with base.open("w") as handle:
        subprocess.run([
            "ncrnaScan", "--norm", "--cdSeq", "--best", "--fold", "-n", "2", "-p", "1.5", "-e", "1.5",
            "-E", "1.5", "-i", "18", "-l", "20", "-S", "1500", "-c", "18", "-m", "2",
            "--primer", "TGGAA", "-u", "4", "-P", "-r", "150", "--fa", genome_fasta,
            "--fai", f"{genome_fasta}.fai", "--bam", sample["bam"],
        ], stdout=handle, check=True)

    # Label identical sequences together, retaining all genomic copies as in groupNcrnaScan.
    data = pd.read_table(base, keep_default_na=False)
    data.columns = BED_COLUMNS + list(data.columns[6:])
    sequence_groups = []
    for number, (_, group) in enumerate(data.groupby("clusterSeq", sort=True), 1):
        group = group.copy()
        group["name"] = [f"ncrnaScan_group{number}" + (f"-iso{i}" if len(group) > 1 else "")
                         for i in range(1, len(group) + 1)]
        sequence_groups.append(group)

    data = pd.concat(sequence_groups, ignore_index=True) if sequence_groups else data
    grouped_path = Path(f"{base}.groupNcRNA")

    # Preserve the old annotated input's name order for tied representative scores.
    data.sort_values("name", kind="stable").to_csv(grouped_path, sep="\t", index=False)
    return grouped_path


def group_candidates(samples, prefix, output_dir):
    """Group calls by reciprocal overlap and select the original representative loci."""
    records = []

    for sample in samples:
        data = pd.read_table(sample["candidates"], keep_default_na=False)
        columns = list(data.columns)
        selected = data[(data["startReadNum"] >= 7) & (data["endReadNum"] >= 7) & (data["clusterLength"] >= 0)]
        for row in selected.to_dict("records"):
            records.append({"row": row, "sample": sample["sample_name"], "score": float(row["startReadNum"])})

    records.sort(key=lambda r: (r["row"]["strand"], r["row"]["#chrom"], r["row"]["start"],
                                -(r["row"]["end"] - r["row"]["start"])))
    groups, current, anchor = [], [], None

    # The longest member is the anchor; this deliberately preserves the original greedy grouping order.
    for record in records:
        row = record["row"]
        overlap = min(row["end"], anchor["end"]) - max(row["start"], anchor["start"]) if anchor else 0
        matched = (anchor is not None and row["#chrom"] == anchor["#chrom"] and row["strand"] == anchor["strand"]
                   and overlap > 0 and overlap / (row["end"] - row["start"]) >= 0.1
                   and overlap / (anchor["end"] - anchor["start"]) >= 0.1)
        if matched:
            current.append(record)
            if row["end"] - row["start"] > anchor["end"] - anchor["start"]:
                anchor = row
        else:
            if current:
                groups.append(current)
            current, anchor = [record], row
    if current:
        groups.append(current)

    representatives, isoforms = [], []
    for group in groups:
        support = sorted({r["sample"] for r in group})
        group_score = sum(r["score"] for r in group)
        if len(support) < 2 or group_score < 20:
            continue

        loci = defaultdict(list)
        for record in group:
            loci[(record["row"]["start"], record["row"]["end"])].append(record)

        supported_loci, passed = [], []
        for locus, members in loci.items():
            if len({r["sample"] for r in members}) < 2:
                continue
            scores = [r["score"] for r in members]
            score_sum, balance = sum(scores), max(scores) / min(scores)
            supported_loci.append(locus)
            if score_sum >= 20 and balance <= 10:
                passed.append((score_sum, balance, locus, members))

        if passed:
            score_sum, balance, locus, members = sorted(passed, key=lambda item: (-item[0], item[1]))[0]
            representative = max(members, key=lambda r: r["score"])
            reason = "passed_multi_rep_exact_loc"
        else:
            representative = max(group, key=lambda r: r["score"])
            score_sum, balance = representative["score"], 1
            locus = (representative["row"]["start"], representative["row"]["end"])
            reason = "no_passed_multi_rep_exact_loc" if supported_loci else "no_multi_rep_exact_loc"

        stats = dict(zip(STATS_COLUMNS, [",".join(support), len(support), group_score, score_sum, balance,
                                         f"{locus[0]}-{locus[1]}", reason, len(group)]))
        group_name = f"ncrnaGroup{len(representatives) + 1}"
        representatives.append({**representative["row"], "name": f"{group_name}_{representative['row']['name']}",
                                "sample": representative["sample"], **stats})
        for i, record in enumerate(sorted(group, key=lambda r: -r["score"]), 1):
            isoforms.append({**record["row"], "name": f"{group_name}-{i}_{record['row']['name']}",
                             "sample": record["sample"], **stats})

    directory = output_dir / f"{prefix}_2samples_startEnd_7reads_ncRNA"
    directory.mkdir(parents=True, exist_ok=True)
    columns += ["sample"] + STATS_COLUMNS
    isoform_path = directory / f"{prefix}.groupNcrna.isoforms.txt"
    pd.DataFrame(isoforms, columns=columns).to_csv(isoform_path, sep="\t", index=False)
    return pd.DataFrame(representatives, columns=columns), directory


def quantify_candidates(data, samples, genome, prefix, directory, by_cell_line, threads):
    """Count fragments and retain replicate-supported groups."""

    gtf = directory.parent / f"{prefix}.filtered_grouped_ncRNA.gtf"
    with gtf.open("w") as handle:
        for chrom, start, end, name, score, strand in data[BED_COLUMNS].itertuples(index=False, name=None):
            handle.write(f'{chrom}\tbed\texon\t{start + 1}\t{end}\t0\t{strand}\t.\tgene_id "{name}";\n')

    count_path = directory / f"{prefix}.filtered_grouped_ncRNA.featureCounts.txt"
    command = ["featureCounts", "-a", gtf, "-t", "exon", "-g", "gene_id", "-p", "--countReadPairs",
               "-s", "1", "-o", count_path, "-T", str(threads), "-M"]

    # Cell-line matrices must contain integer counts for DESeq2; end-type matrices use fractional counts.
    if not by_cell_line:
        command.append("--fraction")
    bams = [str(sample["sorted_bam"]) for sample in samples]
    subprocess.run(command + bams, check=True)

    counts = pd.read_table(count_path, comment="#")[["Geneid"] + bams]
    counts = counts.rename(columns={"Geneid": "name", **dict(zip(bams, [s["count_name"] for s in samples]))})
    data = data.merge(counts, on="name", how="left")

    subsets = [samples] if by_cell_line else [[s for s in samples if s["cell_line"] == c] for c in CELL_LINES[genome]]
    keep = pd.Series(False, index=data.index)
    for subset in subsets:
        values = data[[sample["count_name"] for sample in subset]]
        keep |= (values.gt(0).sum(axis=1) >= 2) & (values.sum(axis=1) >= 20)

    data = data.loc[keep].copy()
    data.to_csv(directory / f"{prefix}.groupNcrna.txt", sep="\t", index=False)
    return data


def run_deseq2(data, samples, cell_line, output_dir, threads):
    """Compare 3P, CAP and OH with 1P using the cell-line integer count matrix."""
    samples = sorted(samples, key=lambda s: (END_TYPES.index(s["end_type"]), int(s["rep"].removeprefix("rep"))))
    names = [sample["count_name"] for sample in samples]
    counts = data.set_index("name")[names].astype("int64").T
    counts.index.name = "sample"
    metadata = pd.DataFrame({"condition": [sample["end_type"] for sample in samples]}, index=names)

    inference = DefaultInference(n_cpus=threads)
    dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors="condition", refit_cooks=True,
                       inference=inference)
    dds.deseq2()

    results = {}
    directory = output_dir / "deg"
    directory.mkdir(parents=True, exist_ok=True)
    for end_type in END_TYPES[1:]:
        stats = DeseqStats(dds, contrast=["condition", end_type, "1P"], inference=inference, alpha=1)
        stats.summary()

        result = stats.results_df.copy()
        result.index.name = "name"
        result.to_csv(directory / f"{cell_line}.1P_{end_type}.pydeseq2_stat_res.txt", sep="\t")
        results[end_type] = result[["log2FoldChange", "padj"]].to_dict("index")
    return results


def best_overlaps(query, target, directory, label, min_coverage):
    """Choose the same-strand hit with greatest minimum reciprocal coverage, then overlap length."""
    query_bed, target_bed = query[BED_COLUMNS].copy(), target[BED_COLUMNS].copy()
    query_bed["name"] = [f"q{i}" for i in range(len(query))]
    target_bed["name"] = [f"t{i}" for i in range(len(target))]

    qpath, tpath, opath = [directory / f"{label}.{suffix}" for suffix in ("query.bed", "target.bed", "overlap.tsv")]
    query_bed.to_csv(qpath, sep="\t", index=False, header=False)
    target_bed.to_csv(tpath, sep="\t", index=False, header=False)
    with opath.open("w") as handle:
        subprocess.run(["bedtools", "intersect", "-wo", "-s", "-a", qpath, "-b", tpath], stdout=handle, check=True)

    best = {}
    with opath.open() as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            qi, ti, overlap = int(fields[3][1:]), int(fields[9][1:]), int(fields[12])
            qcov = overlap / (int(fields[2]) - int(fields[1]))
            tcov = overlap / (int(fields[8]) - int(fields[7]))
            hit = (min(qcov, tcov), overlap, target.iloc[ti]["name"], qcov, tcov)
            if qi not in best or hit[:2] > best[qi][:2]:
                best[qi] = hit
    return {qi: {"name": h[2], "qcov": h[3], "tcov": h[4]} for qi, h in best.items()
            if min_coverage is None or min(h[3], h[4]) >= min_coverage}


def filter_candidates(grouped, cell_grouped, de_results, prefixes, output_dir):
    """Apply the original enrichment and partial-overlap rules to each end type."""
    directory = output_dir / f"{prefixes['1P']}_2samples_startEnd_7reads_ncRNA"
    shutil.copyfile(directory / f"{prefixes['1P']}.groupNcrna.txt",
                    directory / f"{prefixes['1P']}.filtered.groupNcrna.txt")
    for end_type in END_TYPES[1:]:
        prefix = prefixes[end_type]
        directory = output_dir / f"{prefix}_2samples_startEnd_7reads_ncRNA"
        data = grouped[end_type].reset_index(drop=True).copy()
        control_hits = best_overlaps(data, grouped["1P"], directory, "1P", None)
        for field in ("name", "qcov", "tcov"):
            data[f"similar_1P_{field}"] = [control_hits.get(i, {}).get(field, pd.NA) for i in range(len(data))]

        for cell_line, reference in cell_grouped.items():
            hits = best_overlaps(data, reference, directory, cell_line, 0.8)
            annotations = []
            for i in range(len(data)):
                hit = hits.get(i)
                if hit is None:
                    annotations.append([pd.NA] * 5 + [False])
                else:
                    result = de_results[cell_line][end_type][hit["name"]]
                    fold, padj = result["log2FoldChange"], result["padj"]
                    enriched = bool(fold >= 1 and pd.notna(padj) and padj < 0.05)
                    annotations.append([hit["name"], hit["qcov"], hit["tcov"], fold, padj, enriched])
            columns = [f"{cell_line}_{s}" for s in ("ncrnaFam", "qcov", "tcov", "log2FoldChange", "padj", "enriched")]
            data[columns] = pd.DataFrame(annotations, columns=columns, index=data.index)

        keep, reasons = [], []
        for i in range(len(data)):
            hit = control_hits.get(i)
            enriched = any(data.loc[i, f"{cell_line}_enriched"] for cell_line in cell_grouped)
            partial = hit is not None and min(hit["qcov"], hit["tcov"]) < 0.8
            keep.append(partial or enriched)

            # Partial 1P overlap is retained even without enrichment, matching the original locus-support rule.
            reason = "partial_overlap_locus_supported" if partial else "filtered"
            if not partial and enriched:
                reason = "no_control_overlap_but_enriched" if hit is None else "high_reciprocal_coverage_and_enriched"
            reasons.append(reason)

        data["keep"], data["keep_reason"] = keep, reasons
        output = directory / f"{prefix}.filtered.groupNcrna.txt"
        data.loc[data["keep"].astype(bool)].to_csv(output, sep="\t", index=False)


def main():
    """Scan BAMs, group and quantify candidates, then filter each species independently."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample-info", type=Path, default=Path(__file__).with_name("sample_info.tsv"))
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).with_name("group"))
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--genome-fasta", nargs=2, action="append", required=True, metavar=("GENOME", "PATH"))
    args = parser.parse_args()

    sample_info, output_dir = args.sample_info.resolve(), args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    genome_fastas = {genome: Path(path).resolve() for genome, path in args.genome_fasta}

    with sample_info.open(newline="") as handle:
        samples = [row for row in csv.DictReader(handle, delimiter="\t") if row["assay"] == "FLAREseq"]

    for sample in samples:
        directory = sample_info.parent / sample["output_dir"]
        prefix = f"{sample['genome']}_{sample['sample_id']}_{sample['sample_name']}"
        sample.update(directory=directory, prefix=prefix, bam=directory / f"{prefix}.EndToEnd.STAR.Aligned.out.bam",
                      sorted_bam=directory / f"{prefix}.EndToEnd.STAR.Aligned.out.sorted.bam",
                      count_name=f"{sample['cell_line']}_{sample['end_type']}_FLAREseq_{sample['rep']}")
        print(f"Scanning {sample['sample_id']}", flush=True)
        sample["candidates"] = scan_sample(sample, genome_fastas[sample["genome"]])

    for genome in sorted({s["genome"] for s in samples}):
        species_samples = [s for s in samples if s["genome"] == genome and s["cell_line"] in CELL_LINES[genome]]
        prefixes = {end: end if genome == "dm6" else f"{genome}_{end}" for end in END_TYPES}
        grouped, cell_grouped, de_results = {}, {}, {}
        for end_type in END_TYPES:
            subset = [s for s in species_samples if s["end_type"] == end_type]
            prefix = prefixes[end_type]
            print(f"Grouping {genome}/{end_type}", flush=True)
            data, directory = group_candidates(subset, prefix, output_dir)
            ordered = sorted(subset, key=lambda s: (s["cell_line"], int(s["rep"].removeprefix("rep"))))
            grouped[end_type] = quantify_candidates(data, ordered, genome, prefix, directory, False, args.threads)

        for cell_line in CELL_LINES[genome]:
            subset = [s for s in species_samples if s["cell_line"] == cell_line]
            print(f"Grouping and testing {cell_line}", flush=True)
            data, directory = group_candidates(subset, cell_line, output_dir)
            ordered = sorted(subset, key=lambda s: (END_TYPES.index(s["end_type"]), int(s["rep"].removeprefix("rep"))))
            cell_grouped[cell_line] = quantify_candidates(data, ordered, genome, cell_line, directory,
                                                           True, args.threads)
            de_results[cell_line] = run_deseq2(cell_grouped[cell_line], ordered, cell_line, output_dir, args.threads)

        filter_candidates(grouped, cell_grouped, de_results, prefixes, output_dir)


if __name__ == "__main__":
    main()
