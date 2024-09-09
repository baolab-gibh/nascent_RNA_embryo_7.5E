#!/usr/bin/env python3

from collections import OrderedDict
from zipfile import ZipFile
from pathlib import Path

import pandas as pd
import numpy as np
import polars as pls
import matplotlib.pyplot as plt

project_dir = Path("~/Documents/projects/wp_vasaseq").expanduser()
all_batches = ["240409_Lib_embryo", "240612_Lib_28region", "240620_Lib_38region", "240703_Lib_32region", "240710_Lib_37region", "240717_Lib_28region"]

def fetch_total_sequence(qc_archive):
    counts = None
    qc_report = ZipFile(qc_archive, mode="r")
    names_list = [x for x in qc_report.namelist() if x.endswith("/fastqc_data.txt")]
    if len(names_list) == 0: return 0
    qcdata, *_ = names_list
    for per_line in qc_report.read(qcdata).split(b"\n"):
        per_line = per_line.decode()
        if per_line.startswith("Total Sequence"):
            _, counts = per_line.split("\t")
            break
    try:
        qc_report.close()
    except Exception as e:
        print(f"Error while closing the QC archive: {e}")

    return counts


def fetch_alignment_stats(qc_archive):
    report_dict = OrderedDict({"reads_aligned": 0, "exonic_reads": 0, "intronic_reads": 0, "intergenic_reads": 0, "sj_reads": 0})
    with open(qc_archive, "r") as qc_report:
        for line in qc_report:
            if "=" not in line: continue
            _, counts = line.strip().split("=")
            if " reads aligned" in line:
                report_dict["reads_aligned"] = int(counts.strip().replace(",", ""))
            elif " exonic" in line:
                report_dict["exonic_reads"] = int(counts.strip().split(" ")[0].replace(",", ""))
            elif " intronic" in line:
                report_dict["intronic_reads"] = int(counts.strip().split(" ")[0].replace(",", ""))
            elif " intergenic" in line:
                report_dict["intergenic_reads"] = int(counts.strip().split(" ")[0].replace(",", ""))
            elif " overlapping exon" in line:
                report_dict["sj_reads"] = int(counts.strip().split(" ")[0].replace(",", ""))

    return list(report_dict.values())

report = []
for per_batch in all_batches:
    per_batch_qc_dir = project_dir / "outputs/analysis/preprocessing/quality_control" / per_batch
    for per_sample_qc_dir in per_batch_qc_dir.iterdir():
        per_sample = per_sample_qc_dir.stem

        # Original number of reads
        r1_path, r2_path = [per_sample_qc_dir / (per_sample + f".{x}_fastqc.zip") for x in ["R1", "R2"]]
        r1_counts_raw = fetch_total_sequence(r1_path)
        r2_counts_raw = fetch_total_sequence(r2_path)

        # Number of fragments unmerged
        r1_path_unmerged, r2_path_unmerged = [per_sample_qc_dir / (per_sample + f".clean_unmerged.{x}_fastqc.zip") for x in ["R1", "R2"]]
        r1_counts_umerged = fetch_total_sequence(r1_path_unmerged)
        r2_counts_umerged = fetch_total_sequence(r2_path_unmerged)

        # Number of fragments merged by overlapping.
        rx_path = per_sample_qc_dir / f"{per_sample}.clean_merged_fastqc.zip"
        rx_counts_merged = fetch_total_sequence(rx_path)

        # Numbner of alignments
        per_batch_alqc_dir = project_dir / "outputs/analysis/preprocessing/alignment_statistics" / per_batch / per_sample / "rnaseq_qc_results.txt"
        aligned_segs, exonic_segs, intronic_segs, intergenic_segs, sj_segs, = fetch_alignment_stats(per_batch_alqc_dir)

        report.append([per_batch, per_sample, r1_counts_raw, r2_counts_raw, r1_counts_umerged, r2_counts_umerged, rx_counts_merged, aligned_segs, exonic_segs, intronic_segs, intergenic_segs, sj_segs])

new_names = dict(zip([f"column_{x}" for x in range(12)], ["Batch", "Sample", "R1_raw", "R2_raw", "R1_unmerged", "R2_unmerged", "Rx_merged", "Aligned", "Exonic", "Intronic", "Intergenic", "Splicing_junction"]))
ttl_count_report = (pls
    .DataFrame(np.asarray(report))
    .rename(new_names)
    .with_columns(pls.all().exclude("Batch", "Sample").cast(pls.Int32, strict=True))
    .with_columns((pls.col("Sample") + pls.col("Batch")).alias("Sample_id"))
    .with_columns(pls.sum_horizontal("Exonic", "Intronic", "Intergenic").alias("Reads_primary_aligned"))
    .with_columns(pls.sum_horizontal("Rx_merged", "R1_unmerged").alias("Clean"))
    .with_columns((pls.col("Clean") - pls.col("Aligned")).alias("Clean_unaligned"))
    .with_columns((pls.col("R1_raw") - pls.col("Clean")).alias("Contamination"))
    .with_columns((pls.col("Aligned") / pls.col("R1_raw")).alias("Clean_aligned_perc"))
    .with_columns((pls.col("Clean_unaligned") / pls.col("R1_raw")).alias("Clean_unaligned_perc"))
    .with_columns((pls.col("Contamination") / pls.col("R1_raw")).alias("Contamination_perc"))
    .with_columns((pls.col("Exonic") / pls.col("Reads_primary_aligned")).alias("Exonic_perc"))
    .with_columns((pls.col("Intronic") / pls.col("Reads_primary_aligned")).alias("Intronic_perc"))
    .with_columns((pls.col("Intergenic") / pls.col("Reads_primary_aligned")).alias("Intergenic_perc"))
    .with_columns((pls.col("Splicing_junction") / pls.col("Reads_primary_aligned")).alias("Splicing_junction_perc"))
)


# Contamination
ttl_count_report = ttl_count_report.sort("Batch", "R1_raw", "Clean", "Aligned", maintain_order=True)
fig_size = (15, 8)
gridspec_kws = dict(height_ratios=[1, 4])
axe_keys = [['Counts'], ['Percentage']]
fig, axe_dict = plt.subplot_mosaic(axe_keys, gridspec_kw=gridspec_kws, figsize=fig_size, sharex=True)
sample_id = ttl_count_report["Sample_id"]
for axe_name, axe_perse in axe_dict.items():
    axe_perse.xaxis.set_visible(False)
    axe_perse.spines["top"].set_visible(False)
    axe_perse.spines["right"].set_visible(False)
    axe_perse.spines["bottom"].set_visible(False)
    if axe_name == "Counts":
        # sample_values = ttl_count_report["Clean_unaligned", "Clean", "R1_raw"].to_numpy().T
        sample_values = ttl_count_report["Aligned", "Clean_unaligned", "Contamination"].to_numpy().T
        axe_perse.stackplot(sample_id, sample_values)
        axe_perse.set_ylabel("Reads per categoty (#)")
        axe_perse.set_ylim([0, 1.5e7])

    if axe_name == "Percentage":
        # sample_values = ttl_count_report["Exonic_perc", "Intronic_perc", "Intergenic_perc"].to_numpy().T
        sample_values = ttl_count_report["Clean_aligned_perc", "Clean_unaligned_perc", "Contamination_perc"].to_numpy().T
        axe_perse.stackplot(sample_id, sample_values, labels=("UMI", "Clean", "Contamination"))
        axe_perse.set_ylim((-0.025, 1.025))
        axe_perse.set_ylabel("Reads per category (f)")
        axe_perse.legend(loc="upper right", bbox_to_anchor=(.0, .0, 1.1, 0.9))

fig.subplots_adjust(wspace=0, hspace=0.05)
fig.align_labels()
fig.savefig(project_dir / "outputs/analysis/overview/genomic_contamination.pdf")
fig.clear()


# Intronic, exonic, integenic
ttl_count_report = ttl_count_report.sort("Batch", "Reads_primary_aligned", "Exonic", "Intronic", maintain_order=True)
fig_size = (15, 8)
gridspec_kws = dict(height_ratios=[3, 8, 2])
axe_keys = [['Counts'], ['Percentage'], ['SJ_per']]
fig, axe_dict = plt.subplot_mosaic(axe_keys, gridspec_kw=gridspec_kws, figsize=fig_size, sharex=True)
sample_id = ttl_count_report["Sample_id"]
for axe_name, axe_perse in axe_dict.items():
    axe_perse.xaxis.set_visible(False)
    axe_perse.spines["top"].set_visible(False)
    axe_perse.spines["right"].set_visible(False)
    axe_perse.spines["bottom"].set_visible(False)
    if axe_name == "Counts":
        sample_values = ttl_count_report["Exonic", "Intronic", "Intergenic"].to_numpy().T
        axe_perse.stackplot(sample_id, sample_values)
        axe_perse.set_ylabel("Aligned reads (#)")

    if axe_name == "Percentage":
        sample_values = ttl_count_report["Exonic_perc", "Intronic_perc", "Intergenic_perc"].to_numpy().T
        axe_perse.stackplot(sample_id, sample_values, labels=("Exonic", "Intronic", "Intergenic"))
        axe_perse.set_ylim((-0.025, 1.025))
        axe_perse.set_ylabel("Reads in genomic regions (f)")
        axe_perse.legend(loc="upper right", bbox_to_anchor=(.0, .0, 1.1, 0.9))

    if axe_name == "SJ_per":
        sample_values = ttl_count_report["Splicing_junction_perc"]
        axe_perse.bar(sample_id, sample_values, color="0.1")
        axe_perse.set_ylim((0.175, 0))
        axe_perse.set_ylabel("SJ reads (f)")

fig.subplots_adjust(wspace=0, hspace=0.05)
fig.align_labels()
fig.savefig(project_dir / "outputs/analysis/overview/reads_genomic_distribution.pdf")
fig.clear()

ttl_count_report.write_csv(project_dir / "outputs/analysis/overview/sequencing_statistics.csv")



# Peng etal 2019 Nature
report = []
per_batch_qc_dir = project_dir / "outputs/analysis/preprocessing/2019_peng_etal_nature/quality_control"
for per_sample_qc_dir in per_batch_qc_dir.iterdir():
    per_sample = per_sample_qc_dir.stem

    # Original number of reads
    r1_path, r2_path = [per_sample_qc_dir / (per_sample + f".{x}_fastqc.zip") for x in ["R1", "R2"]]
    r1_counts_raw = fetch_total_sequence(r1_path)
    r2_counts_raw = fetch_total_sequence(r2_path)

    # Number of fragments unmerged
    r1_path_unmerged, r2_path_unmerged = [per_sample_qc_dir / (per_sample + f".{x}_fastqc.zip") for x in ["R1", "R2"]]
    r1_counts_umerged = fetch_total_sequence(r1_path_unmerged)
    r2_counts_umerged = fetch_total_sequence(r2_path_unmerged)

    # Numbner of alignments
    per_batch_alqc_dir = project_dir / "outputs/analysis/preprocessing/2019_peng_etal_nature/mapping_reports" / per_sample / "rnaseq_qc_results.txt"
    aligned_segs, exonic_segs, intronic_segs, intergenic_segs, sj_segs, = fetch_alignment_stats(per_batch_alqc_dir)

    report.append([per_sample, r1_counts_raw, r2_counts_raw, r1_counts_umerged, r2_counts_umerged, aligned_segs, exonic_segs, intronic_segs, intergenic_segs, sj_segs])

new_names = dict(zip([f"column_{x}" for x in range(12)], ["Sample", "R1_raw", "R2_raw", "R1_unmerged", "R2_unmerged", "Aligned", "Exonic", "Intronic", "Intergenic", "Splicing_junction"]))
ttl_count_report = (pls.DataFrame(np.asarray(report)).rename(new_names)
    .with_columns(pls.all().exclude("Sample").cast(pls.Int32, strict=True))
    .with_columns(pls.sum_horizontal("Exonic", "Intronic", "Intergenic").alias("Reads_primary_aligned"))
    .with_columns(pls.sum_horizontal("R1_unmerged").alias("Clean"))
    .with_columns((pls.col("Clean") - pls.col("Aligned")).alias("Clean_unaligned"))
    .with_columns((pls.col("R1_raw") - pls.col("Clean")).alias("Contamination"))
    .with_columns((pls.col("Aligned") / pls.col("R1_raw")).alias("Clean_aligned_perc"))
    .with_columns((pls.col("Clean_unaligned") / pls.col("R1_raw")).alias("Clean_unaligned_perc"))
    .with_columns((pls.col("Contamination") / pls.col("R1_raw")).alias("Contamination_perc"))
    .with_columns((pls.col("Exonic") / pls.col("Reads_primary_aligned")).alias("Exonic_perc"))
    .with_columns((pls.col("Intronic") / pls.col("Reads_primary_aligned")).alias("Intronic_perc"))
    .with_columns((pls.col("Intergenic") / pls.col("Reads_primary_aligned")).alias("Intergenic_perc"))
    .with_columns((pls.col("Splicing_junction") / pls.col("Reads_primary_aligned")).alias("Splicing_junction_perc"))
)


ttl_count_report = ttl_count_report.sort("Reads_primary_aligned", "Exonic", "Intronic", maintain_order=True)
ttl_count_report.write_csv(project_dir / "outputs/analysis/overview/2019_peng_etal_nature.sequencing_statistics.csv")
fig_size = (10, 6)
gridspec_kws = dict(height_ratios=[3, 8, 2])
axe_keys = [['Counts'], ['Percentage'], ['SJ_per']]
fig, axe_dict = plt.subplot_mosaic(axe_keys, gridspec_kw=gridspec_kws, figsize=fig_size, sharex=True)
sample_id = ttl_count_report["Sample"]
for axe_name, axe_perse in axe_dict.items():
    axe_perse.xaxis.set_visible(False)
    axe_perse.spines["top"].set_visible(False)
    axe_perse.spines["right"].set_visible(False)
    axe_perse.spines["bottom"].set_visible(False)
    if axe_name == "Counts":
        sample_values = ttl_count_report["Exonic", "Intronic", "Intergenic"].to_numpy().T
        axe_perse.stackplot(sample_id, sample_values)
        axe_perse.set_ylabel("Aligned reads (#)")

    if axe_name == "Percentage":
        sample_values = ttl_count_report["Exonic_perc", "Intronic_perc", "Intergenic_perc"].to_numpy().T
        axe_perse.stackplot(sample_id, sample_values, labels=("Exonic", "Intronic", "Intergenic"))
        axe_perse.set_ylim((-0.025, 1.025))
        axe_perse.set_ylabel("Reads in genomic regions (f)")
        axe_perse.legend(loc="upper right", bbox_to_anchor=(.0, .0, 1.125, 0.9))

    if axe_name == "SJ_per":
        sample_values = ttl_count_report["Splicing_junction_perc"]
        axe_perse.bar(sample_id, sample_values, color="0.1")
        axe_perse.set_ylim((0.175, 0))
        axe_perse.set_ylabel("SJ reads (f)")

fig.subplots_adjust(wspace=0, hspace=0.05)
fig.align_labels()
fig.savefig(project_dir / "outputs/analysis/overview/2019_peng_etal_nature.reads_genomic_distribution.pdf")
fig.clear()
