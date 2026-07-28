#!/usr/bin/env python
"""
Draft for pre-clusetring steps

Last checked by: BP on going dev
"""

import argparse
import io
import shlex
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


OUTPUT_FILES = {
    "query_features": "query_features.tsv",
    "clusters": "clusters.tsv",
    "representatives": "representatives.tsv",
    "cluster_manifest": "cluster_manifest.tsv",
    "non_aligning": "non_aligning_queries.tsv",
    "tree": "hdbscan_tree.tsv",
}
DEFAULT_MIN_CLUSTER_SIZE = 5
DEFAULT_SKANI_BATCH_SIZE = 1000


def load_hdbscan():
    import hdbscan

    return hdbscan


def parse_args():
    parser = argparse.ArgumentParser(
        description="Precluster query assemblies against a reference using skani features and HDBSCAN."
    )
    parser.add_argument("--reference", required=True, help="Reference fasta file")
    parser.add_argument("--query-list", required=True, help="Text file with one query fasta path per line")
    parser.add_argument("--threads", type=int, default=1, help="Threads to pass to skani")
    parser.add_argument(
        "--skani-args",
        default="",
        help="Additional raw arguments to pass to skani dist",
    )
    return parser.parse_args()


def read_query_paths(query_list_path):
    # load one query path per line
    query_paths = []
    with open(query_list_path, "r", encoding="utf-8") as handle:
        for line in handle:
            query_path = line.strip()
            if query_path:
                query_paths.append(query_path)
    if not query_paths:
        raise ValueError("No query paths were provided for preclustering.")
    return query_paths


def get_query_id(query_path):
    return Path(query_path).stem


def chunk_iterable(items, chunk_size):
    # split large query lists into smaller skani calls
    for index in range(0, len(items), chunk_size):
        yield items[index:index + chunk_size]


def run_skani(reference_path, query_paths, threads, skani_args):
    # run skani in batches so large query sets do not create huge commands
    extra_args = shlex.split(skani_args) if skani_args else []
    chunk_frames = []

    for chunk in chunk_iterable(query_paths, DEFAULT_SKANI_BATCH_SIZE):
        command = [
            "skani",
            "dist",
            "-t",
            str(max(1, threads)),
            "-r",
            reference_path,
            "-q",
            *chunk,
            *extra_args,
        ]
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
        stdout = result.stdout.strip()
        if not stdout:
            continue
        chunk_frames.append(pd.read_csv(io.StringIO(stdout), sep="\t"))

    if not chunk_frames:
        return pd.DataFrame(
            columns=[
                "Ref_file",
                "Query_file",
                "ANI",
                "Align_fraction_ref",
                "Align_fraction_query",
                "Ref_name",
                "Query_name",
            ]
        )

    skani_df = pd.concat(chunk_frames, ignore_index=True)
    return skani_df


def minmax_scale(series):
    # keep all features on the same scale for clustering
    min_value = series.min()
    max_value = series.max()
    if pd.isna(min_value) or pd.isna(max_value) or np.isclose(min_value, max_value):
        return pd.Series(np.zeros(len(series), dtype=float), index=series.index)
    return (series - min_value) / (max_value - min_value)


def build_feature_table(skani_df, query_paths):
    # start from the full query set so non-aligning queries can be tracked explicitly
    all_queries_df = pd.DataFrame(
        {
            "query_path": query_paths,
            "query_id": [get_query_id(path) for path in query_paths],
        }
    )

    if skani_df.empty:
        aligned_df = pd.DataFrame(
            columns=[
                "query_path",
                "query_id",
                "ANI",
                "AF_ref",
                "AF_query",
                "ANI_scaled",
                "AF_ref_scaled",
                "AF_query_scaled",
            ]
        )
        non_aligning_df = all_queries_df.copy()
        non_aligning_df["reason"] = "no_skani_output"
        return aligned_df, non_aligning_df

    aligned_df = skani_df.rename(
        columns={
            "Query_file": "query_path",
            "ANI": "ANI",
            "Align_fraction_ref": "AF_ref",
            "Align_fraction_query": "AF_query",
        }
    ).copy()
    aligned_df["query_id"] = aligned_df["query_path"].map(get_query_id)

    aligned_df = (
        # keep the strongest skani hit per query
        aligned_df.sort_values(
            by=["query_path", "ANI", "AF_ref", "AF_query"],
            ascending=[True, False, False, False],
        )
        .drop_duplicates(subset=["query_path"], keep="first")
        .loc[:, ["query_path", "query_id", "ANI", "AF_ref", "AF_query"]]
    )

    aligned_df["ANI_scaled"] = minmax_scale(aligned_df["ANI"].astype(float))
    aligned_df["AF_ref_scaled"] = minmax_scale(aligned_df["AF_ref"].astype(float))
    aligned_df["AF_query_scaled"] = minmax_scale(aligned_df["AF_query"].astype(float))

    non_aligning_df = all_queries_df.loc[
        ~all_queries_df["query_path"].isin(aligned_df["query_path"])
    ].copy()
    non_aligning_df["reason"] = "no_skani_output"

    return aligned_df, non_aligning_df


def find_selected_descendants(cluster_id, children_map, selected_cluster_ids, cache):
    # cache selected descendants so repeated tree walks stay cheap
    if cluster_id in cache:
        return cache[cluster_id]

    descendants = set()
    for child_id in children_map.get(cluster_id, []):
        if child_id in selected_cluster_ids:
            descendants.add(child_id)
        descendants.update(
            find_selected_descendants(child_id, children_map, selected_cluster_ids, cache)
        )
    cache[cluster_id] = descendants
    return descendants


def assign_noise_by_tree(clusterer):
    hdbscan = load_hdbscan()
    # use the condensed tree first so noise reassignment follows hdbscan structure
    tree_df = clusterer.condensed_tree_.to_pandas().copy()
    tree_df.to_csv(OUTPUT_FILES["tree"], sep="\t", index=False)

    prediction_data = getattr(clusterer, "prediction_data_", None)
    cluster_map = getattr(prediction_data, "cluster_map", {}) or {}
    selected_cluster_ids = set(cluster_map.keys())
    cluster_stability = {
        cluster_id: clusterer.cluster_persistence_[label]
        for cluster_id, label in cluster_map.items()
        if label < len(clusterer.cluster_persistence_)
    }

    internal_edges = tree_df.loc[tree_df["child_size"] > 1, ["parent", "child"]]
    children_map = {}
    for row in internal_edges.itertuples(index=False):
        children_map.setdefault(int(row.parent), []).append(int(row.child))

    descendant_cache = {}
    labels = clusterer.labels_.copy()
    assignment_method = np.full(len(labels), "hdbscan_flat", dtype=object)
    membership_vectors = None

    for point_index, label in enumerate(labels):
        if label != -1:
            continue

        # walk up the tree from each noise point until a selected cluster can be resolved
        point_rows = tree_df.loc[tree_df["child"] == point_index].sort_values(
            by="lambda_val",
            ascending=False,
        )
        assigned_label = None

        for row in point_rows.itertuples(index=False):
            parent_id = int(row.parent)
            if parent_id in cluster_map:
                assigned_label = int(cluster_map[parent_id])
                assignment_method[point_index] = "condensed_tree_parent"
                break

            descendant_clusters = find_selected_descendants(
                parent_id,
                children_map,
                selected_cluster_ids,
                descendant_cache,
            )
            if len(descendant_clusters) == 1:
                selected_cluster_id = next(iter(descendant_clusters))
                assigned_label = int(cluster_map[selected_cluster_id])
                assignment_method[point_index] = "condensed_tree_descendant"
                break
            if len(descendant_clusters) > 1:
                ranked_candidates = sorted(
                    descendant_clusters,
                    key=lambda cluster_id: cluster_stability.get(cluster_id, 0.0),
                    reverse=True,
                )
                if ranked_candidates:
                    assigned_label = int(cluster_map[ranked_candidates[0]])
                    assignment_method[point_index] = "condensed_tree_stability"
                    break

        if assigned_label is None:
            # fall back to soft membership only if the tree walk is ambiguous
            if membership_vectors is None:
                membership_vectors = hdbscan.all_points_membership_vectors(clusterer)
            if membership_vectors.shape[1] > 0:
                assigned_label = int(np.argmax(membership_vectors[point_index]))
                assignment_method[point_index] = "soft_membership"
            else:
                assigned_label = 0
                assignment_method[point_index] = "single_cluster_fallback"

        labels[point_index] = assigned_label

    return labels, assignment_method, tree_df


def fit_clusters(feature_df, min_cluster_size):
    hdbscan = load_hdbscan()
    tree_df = pd.DataFrame(columns=["parent", "child", "lambda_val", "child_size"])
    if feature_df.empty:
        return np.array([], dtype=int), np.array([], dtype=float), np.array([], dtype=object), tree_df

    if len(feature_df) < max(2, min_cluster_size):
        # keep tiny aligned sets together instead of forcing hdbscan on too few points
        labels = np.zeros(len(feature_df), dtype=int)
        probabilities = np.ones(len(feature_df), dtype=float)
        assignment_method = np.full(len(feature_df), "small_dataset_fallback", dtype=object)
        return labels, probabilities, assignment_method, tree_df

    features = feature_df[["ANI_scaled", "AF_ref_scaled", "AF_query_scaled"]].to_numpy(dtype=float)
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=None,
        metric="euclidean",
        cluster_selection_method="eom",
        allow_single_cluster=True,
        prediction_data=True,
    )
    clusterer.fit(features)

    labels, assignment_method, tree_df = assign_noise_by_tree(clusterer)
    probabilities = clusterer.probabilities_.copy()

    if np.all(labels == -1):
        labels = np.zeros(len(feature_df), dtype=int)
        assignment_method = np.full(len(feature_df), "all_noise_fallback", dtype=object)

    return labels, probabilities, assignment_method, tree_df


def assign_cluster_ids(feature_df, labels):
    # convert numeric labels into stable user-facing cluster ids
    label_order = sorted(pd.unique(labels))
    label_to_cluster_id = {
        label: f"cluster_{index:04d}"
        for index, label in enumerate(label_order, start=1)
    }
    feature_df["cluster_label"] = labels
    feature_df["cluster_id"] = feature_df["cluster_label"].map(label_to_cluster_id)
    return label_to_cluster_id


def select_representatives(feature_df):
    # pick one representative per cluster using the smallest mean feature-space distance
    representatives = []
    for cluster_id, group in feature_df.groupby("cluster_id", sort=True):
        if len(group) == 1:
            representative = group.iloc[0]
            mean_distance = 0.0
        else:
            coords = group[["ANI_scaled", "AF_ref_scaled", "AF_query_scaled"]].to_numpy(dtype=float)
            pairwise_distances = np.linalg.norm(
                coords[:, None, :] - coords[None, :, :],
                axis=2,
            )
            mean_distances = pairwise_distances.mean(axis=1)
            candidate_order = np.lexsort((group["query_id"].to_numpy(), mean_distances))
            representative = group.iloc[int(candidate_order[0])]
            mean_distance = float(mean_distances[int(candidate_order[0])])

        representatives.append(
            {
                "cluster_id": cluster_id,
                "cluster_size": int(len(group)),
                "representative_query_id": representative["query_id"],
                "representative_query_path": representative["query_path"],
                "representative_mean_distance": mean_distance,
                "mean_ANI": float(group["ANI"].mean()),
                "mean_AF_ref": float(group["AF_ref"].mean()),
                "mean_AF_query": float(group["AF_query"].mean()),
            }
        )
    representatives_df = pd.DataFrame(representatives).sort_values(
        by=["cluster_size", "cluster_id"],
        ascending=[False, True],
    )
    return representatives_df


def write_empty_tree():
    # emit an empty tree file so downstream contracts stay stable
    pd.DataFrame(columns=["parent", "child", "lambda_val", "child_size"]).to_csv(
        OUTPUT_FILES["tree"],
        sep="\t",
        index=False,
    )


def main():
    args = parse_args()
    query_paths = read_query_paths(args.query_list)

    skani_df = run_skani(
        reference_path=args.reference,
        query_paths=query_paths,
        threads=args.threads,
        skani_args=args.skani_args,
    )

    aligned_df, non_aligning_df = build_feature_table(skani_df, query_paths)

    if aligned_df.empty:
        # all queries fall into the special non-aligning bucket
        clusters_df = pd.DataFrame(
            columns=[
                "query_id",
                "query_path",
                "cluster_label",
                "cluster_id",
                "hdbscan_probability",
                "assignment_method",
                "ANI",
                "AF_ref",
                "AF_query",
                "ANI_scaled",
                "AF_ref_scaled",
                "AF_query_scaled",
            ]
        )
        representatives_df = pd.DataFrame(
            columns=[
                "cluster_id",
                "cluster_size",
                "representative_query_id",
                "representative_query_path",
                "representative_mean_distance",
                "mean_ANI",
                "mean_AF_ref",
                "mean_AF_query",
            ]
        )
        cluster_manifest_df = pd.DataFrame(
            [
                {
                    "cluster_id": "non_aligning",
                    "cluster_size": int(len(non_aligning_df)),
                    "representative_query_id": "",
                    "representative_query_path": "",
                }
            ]
        )
        write_empty_tree()
    else:
        # cluster aligned queries and assign every point to exactly one cluster
        labels, probabilities, assignment_method, tree_df = fit_clusters(
            aligned_df,
            min_cluster_size=DEFAULT_MIN_CLUSTER_SIZE,
        )
        if tree_df.empty:
            write_empty_tree()

        assign_cluster_ids(aligned_df, labels)
        aligned_df["hdbscan_probability"] = probabilities
        aligned_df["assignment_method"] = assignment_method

        representatives_df = select_representatives(aligned_df)
        representative_map = representatives_df.set_index("cluster_id")[
            "representative_query_id"
        ].to_dict()

        # attach representative ids to each cluster member for downstream reporting
        clusters_df = aligned_df.copy()
        clusters_df["representative_query_id"] = clusters_df["cluster_id"].map(representative_map)
        cluster_manifest_df = representatives_df.loc[
            :,
            [
                "cluster_id",
                "cluster_size",
                "representative_query_id",
                "representative_query_path",
            ],
        ].copy()
        if not non_aligning_df.empty:
            cluster_manifest_df = pd.concat(
                [
                    cluster_manifest_df,
                    pd.DataFrame(
                        [
                            {
                                "cluster_id": "non_aligning",
                                "cluster_size": int(len(non_aligning_df)),
                                "representative_query_id": "",
                                "representative_query_path": "",
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )

    aligned_df.to_csv(OUTPUT_FILES["query_features"], sep="\t", index=False)
    clusters_df.to_csv(OUTPUT_FILES["clusters"], sep="\t", index=False)
    representatives_df.to_csv(OUTPUT_FILES["representatives"], sep="\t", index=False)
    cluster_manifest_df.to_csv(OUTPUT_FILES["cluster_manifest"], sep="\t", index=False)
    non_aligning_df.to_csv(OUTPUT_FILES["non_aligning"], sep="\t", index=False)


if __name__ == "__main__":
    main()
