#!/usr/bin/env python
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
CGVIEW_SCRIPT = SCRIPT_DIR / "cgviewBuilderPy.py"
SCIFR_SCRIPT = SCRIPT_DIR / "generateSCIFR.py"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate representative and per-cluster BLITSFR reports from preclustered queries."
    )
    parser.add_argument("--compiled-results", required=True, help="Compiled BLAST results TSV")
    parser.add_argument("--coverage", required=True, help="BLAST coverage TSV")
    parser.add_argument("--clusters", required=True, help="Cluster assignments TSV")
    parser.add_argument("--representatives", required=True, help="Cluster representatives TSV")
    parser.add_argument("--cluster-manifest", required=True, help="Cluster manifest TSV")
    parser.add_argument("--non-aligning", required=True, help="Non-aligning query TSV")
    parser.add_argument("--reference", required=True, help="Reference fasta")
    parser.add_argument("--features", help="Reference features GFF3")
    parser.add_argument("--metadata", help="Optional metadata TSV")
    parser.add_argument("--template", required=True, help="BLITSFR HTML template")
    parser.add_argument("--params-json", required=True, help="Pipeline params JSON")
    parser.add_argument("--pipeline-version", required=True, help="Pipeline version string")
    parser.add_argument("--report-title", required=True, help="Base report title")
    parser.add_argument("--save-intermediate", action="store_true", default=False)
    return parser.parse_args()


def normalize_optional_path(path_value):
    if not path_value:
        return None
    path = Path(path_value)
    if path.name == "[]":
        return None
    return path


def subset_metadata(metadata_path, query_ids, output_path):
    if metadata_path is None:
        return None

    metadata_df = pd.read_csv(metadata_path, sep="\t")
    if "id" not in metadata_df.columns:
        raise ValueError("Metadata file must contain an 'id' column.")

    subset_df = metadata_df.loc[metadata_df["id"].isin(query_ids)].copy()
    if subset_df.empty:
        return None

    subset_df.to_csv(output_path, sep="\t", index=False)
    return output_path


def run_command(command):
    subprocess.run(command, check=True)


def build_report(
    report_dir,
    report_name,
    compiled_subset,
    coverage_subset,
    metadata_subset,
    args,
):
    cgview_json = report_dir / "cgview.json"
    html_output = report_dir / "blitsfr.html"
    json_output = report_dir / "blitsfr.json"

    cgview_command = [
        sys.executable,
        str(CGVIEW_SCRIPT),
        "--name",
        report_name,
        "--sequence",
        str(args.reference),
        "--blast",
        str(compiled_subset),
        "--output",
        str(cgview_json),
    ]
    if args.features:
        cgview_command.extend(["--features", str(args.features)])
    run_command(cgview_command)

    scifr_command = [
        sys.executable,
        str(SCIFR_SCRIPT),
        "-c",
        str(cgview_json),
        "-d",
        str(coverage_subset),
        "-t",
        str(args.template),
        "-o",
        str(html_output),
        "-j",
        str(json_output),
        "-l",
        str(args.params_json),
        "-p",
        args.pipeline_version,
    ]
    if metadata_subset is not None:
        scifr_command.extend(["-m", str(metadata_subset)])
    if args.save_intermediate:
        scifr_command.append("--save_intermediate")
    run_command(scifr_command)


def write_subset_tables(compiled_df, coverage_df, query_ids, output_dir):
    compiled_subset = output_dir / "compiled_results.tsv"
    coverage_subset = output_dir / "blast_coverage.tsv"

    compiled_filtered = compiled_df.loc[compiled_df["query_file"].isin(query_ids)].copy()
    coverage_filtered = coverage_df.loc[coverage_df["query_file"].isin(query_ids)].copy()

    compiled_filtered.to_csv(compiled_subset, sep="\t", index=False)
    coverage_filtered.to_csv(coverage_subset, sep="\t", index=False)

    return compiled_subset, coverage_subset


def main():
    args = parse_args()
    args.reference = Path(args.reference)
    args.template = Path(args.template)
    args.params_json = Path(args.params_json)
    args.features = normalize_optional_path(args.features)
    metadata_path = normalize_optional_path(args.metadata)

    compiled_df = pd.read_csv(args.compiled_results, sep="\t")
    coverage_df = pd.read_csv(args.coverage, sep="\t")
    clusters_df = pd.read_csv(args.clusters, sep="\t")
    representatives_df = pd.read_csv(args.representatives, sep="\t")
    cluster_manifest_df = pd.read_csv(args.cluster_manifest, sep="\t")
    non_aligning_df = pd.read_csv(args.non_aligning, sep="\t")

    if "query_file" not in compiled_df.columns:
        raise ValueError("Compiled BLAST results must contain a 'query_file' column.")
    if "query_file" not in coverage_df.columns:
        raise ValueError("Coverage TSV must contain a 'query_file' column.")

    report_manifest_rows = []

    representative_ids = representatives_df["representative_query_id"].tolist()
    representative_compiled, representative_coverage = write_subset_tables(
        compiled_df,
        coverage_df,
        representative_ids,
        Path("."),
    )
    representative_metadata = subset_metadata(
        metadata_path,
        representative_ids,
        Path("metadata.tsv"),
    )
    build_report(
        Path("."),
        f"{args.report_title} [representatives]",
        representative_compiled,
        representative_coverage,
        representative_metadata,
        args,
    )
    report_manifest_rows.append(
        {
            "report_id": "representatives",
            "report_type": "representatives",
            "cluster_id": "",
            "query_count": len(representative_ids),
            "report_html": "blitsfr.html",
        }
    )

    clusters_root = Path("clusters")
    clusters_root.mkdir(exist_ok=True)

    for cluster_row in cluster_manifest_df.itertuples(index=False):
        cluster_id = str(cluster_row.cluster_id)
        if cluster_id == "non_aligning":
            continue

        member_ids = clusters_df.loc[clusters_df["cluster_id"] == cluster_id, "query_id"].tolist()
        if not member_ids:
            continue

        cluster_dir = clusters_root / cluster_id
        cluster_dir.mkdir(parents=True, exist_ok=True)

        cluster_compiled, cluster_coverage = write_subset_tables(
            compiled_df,
            coverage_df,
            member_ids,
            cluster_dir,
        )
        cluster_metadata = subset_metadata(
            metadata_path,
            member_ids,
            cluster_dir / "metadata.tsv",
        )
        build_report(
            cluster_dir,
            f"{args.report_title} [{cluster_id}]",
            cluster_compiled,
            cluster_coverage,
            cluster_metadata,
            args,
        )
        report_manifest_rows.append(
            {
                "report_id": cluster_id,
                "report_type": "cluster_members",
                "cluster_id": cluster_id,
                "query_count": len(member_ids),
                "report_html": str(cluster_dir / "blitsfr.html"),
            }
        )

    if not non_aligning_df.empty:
        non_aligning_dir = Path("non_aligning")
        non_aligning_dir.mkdir(exist_ok=True)
        shutil.copyfile(args.non_aligning, non_aligning_dir / "non_aligning_queries.tsv")
        report_manifest_rows.append(
            {
                "report_id": "non_aligning",
                "report_type": "non_aligning",
                "cluster_id": "non_aligning",
                "query_count": int(len(non_aligning_df)),
                "report_html": "",
            }
        )

    pd.DataFrame(report_manifest_rows).to_csv("report_manifest.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
