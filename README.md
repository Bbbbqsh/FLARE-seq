# FLARE-seq and FLARE-SHAPE-MaP analysis

This directory contains the core scripts for read preprocessing, ncRNA candidate identification, and RNA structure analysis. Raw sequencing data can be downloaded from GEO.

The workflow starts with shared preprocessing, then follows two branches:

```text
preprocess.py
├── FLARE-seq → flareseq_process_data.py → candidate_ncrna_identification.py
└── FLARE-SHAPE-MaP → flareshapemap_process_data.py
```

## 1. Read preprocessing — `preprocess.py`

Processes paired-end reads from both assays using cutadapt, collapseRead, and barcode trimming. FLARE-seq and FLARE-SHAPE-MaP are handled by separate functions, removing 4 nt and 6 nt barcodes from both ends of each mate, respectively.

Output: preprocessed paired-end FASTQ files ending in `.cutadapt.collapse.cutBarcodes.fq.gz`.

## 2. FLARE-seq alignment and signal tracks — `flareseq_process_data.py`

Builds STAR genome indexes, aligns preprocessed reads, and sorts and indexes the resulting BAM files. The in-house bamToWig tool and UCSC wigToBigWig generate six RPM-normalized tracks: coverage, 5′ start, and 3′ end signals for each strand.

Output: BAM files and six bigWig files per sample.

## 3. Candidate ncRNA identification — `candidate_ncrna_identification.py`

Uses ncrnaScan to identify candidates from FLARE-seq BAM files, then groups candidates across samples and selects representative loci. Candidates require terminal-read and replicate support. featureCounts provides quantification, and DESeq2 enrichment tests together with overlap rules filter candidates relative to the 1P control.

Output: `*.filtered.groupNcrna.txt` for 1P, 3P, CAP, and OH, processed separately for dm6 and hg38.

## 4. FLARE-SHAPE-MaP reactivity and structure analysis — `flareshapemap_process_data.py`

Builds a Bowtie2 index for the target transcriptome and uses RNAFramework to map preprocessed reads, count mutations and coverage, and normalize reactivities from matched DMSO/NAIN3 samples. Reactivities are combined for each cell-line source and for all replicates, then exported as three-column SHAPE files.

Structure prediction uses the separate normalized replicates from all sources in multi-experiment mode for every RNA. Results are written to `all_multi_experiment`.

This step requires a target transcript table with `name` and `seq` columns, supplied via `--reference-table`; the filtered candidate table is not used directly as this input.

## Inputs and usage

- `sample_info.tsv` contains metadata for both assays. Edit `r1`, `r2`, and `output_dir` to match your data locations. Relative paths in the table are resolved against the table's directory.
- Processing scripts read preprocessed FASTQs or BAMs from each sample's `output_dir`.
- Supply reference files and tool locations through the script arguments. Genome FASTA files require a matching `.fai` index.
- Install the dependencies listed in each script's opening docstring and make command-line tools available on `PATH`. In-house tools such as collapseRead, bamToWig, and ncrnaScan must be supplied separately.
- Use `python SCRIPT.py --help` to view the available arguments.
