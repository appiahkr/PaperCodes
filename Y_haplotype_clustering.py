#!/usr/bin/env python3
"""
Y-linked haplotype clustering and visualization pipeline.

This script:
1. Reads a list of Y-linked SNP positions
2. Extracts corresponding genotypes from a Y chromosome VCF
3. Filters to male samples with complete data
4. Computes pairwise Hamming distances between haplotypes
5. Performs hierarchical clustering and selects the optimal number of clusters
6. Visualizes haplotype structure using PCA and MDS
7. Identifies shared SNPs among clusters and visualizes overlap with an UpSet plot

"""

import sys
import pandas as pd
import allel
import numpy as np

from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.metrics import silhouette_score

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster

from plotnine import *
import matplotlib.pyplot as plt
import seaborn as sns

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from upsetplot import from_contents, UpSet


# ---------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------
def main():
    """
    Parse command-line arguments and run the haplotype analysis.
    """
    usage = 'usage:' + sys.argv[1] + "<Y_linked_snps_file> <chrYSDR.recode.vcf>"
    if len(sys.argv) != 3:
        print(usage)
        sys.exit(1)

    Y_linked_snps_file = sys.argv[1]
    chrYSDR_vcf = sys.argv[2]

    compute_hamming_distance(Y_linked_snps_file, chrYSDR_vcf)


# ---------------------------------------------------------------------
# Analysis function
# ---------------------------------------------------------------------
def compute_hamming_distance(Y_linked_snps_file, chrYSDR_vcf):
    """
    Compute pairwise Hamming distances between Y-linked haplotypes,
    cluster samples, and generate PCA/MDS visualizations.

    Parameters
    ----------
    Y_linked_snps_file : str
        Tab-delimited file containing Y-linked SNP positions (must include 'POS').
    chrYSDR_vcf : str
        VCF file containing Y chromosome genotypes.
    """

    # -----------------------------------------------------------------
    # Load Y-linked SNP positions
    # -----------------------------------------------------------------
    y_snps = pd.read_csv(Y_linked_snps_file, sep='\t')
    y_positions = set(y_snps['POS'].tolist())

    # -----------------------------------------------------------------
    # Read VCF (only required fields for efficiency)
    # -----------------------------------------------------------------
    callset = allel.read_vcf(
        chrYSDR_vcf,
        fields=['samples', 'variants/CHROM', 'variants/POS', 'calldata/GT']
    )

    samples = callset['samples']

    # Identify male samples (assumes "female" appears in female sample IDs)
    male_samples = [s for s in samples if 'female' not in s.lower()]
    male_indices = [i for i, s in enumerate(samples) if s in male_samples]

    chroms = callset['variants/CHROM']
    positions = callset['variants/POS']

    # Filter to Y chromosome SNPs that are in the predefined Y-linked list
    is_y_and_linked = [
        (ch == 'Y' and pos in y_positions)
        for ch, pos in zip(chroms, positions)
    ]

    # Extract genotypes for selected variants and male samples
    genotypes = callset['calldata/GT'][is_y_and_linked, :, :]
    filtered_positions = positions[is_y_and_linked]
    genotypes_male = genotypes[:, male_indices, :]

    # Convert diploid GT to haplotypes (Y is haploid → take first allele)
    haplotypes = genotypes_male[:, :, 0].T  # shape: (samples, variants)

    # -----------------------------------------------------------------
    # Remove samples with any missing genotypes
    # -----------------------------------------------------------------
    mask_no_missing = np.all(haplotypes != -1, axis=1)
    haplotypes = haplotypes[mask_no_missing, :]
    male_samples = np.array(male_samples)[mask_no_missing]

    print(f"Kept {len(male_samples)} male samples with no missing data.")

    # -----------------------------------------------------------------
    # Optional: convert numeric alleles to bases (useful for FASTA output)
    # -----------------------------------------------------------------
    allele_map = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    haplotypes_bases = np.vectorize(allele_map.get)(haplotypes)

    records = [
        SeqRecord(Seq(''.join(row)), id=sample, description='')
        for sample, row in zip(male_samples, haplotypes_bases)
    ]

    # -----------------------------------------------------------------
    # Compute pairwise Hamming distance matrix
    # -----------------------------------------------------------------
    def hamming_proportion(u, v):
        """Proportion of sites at which two haplotypes differ."""
        return np.sum(u != v) / len(u)

    dist_condensed = pdist(haplotypes, metric=hamming_proportion)

    # -----------------------------------------------------------------
    # Hierarchical clustering
    # -----------------------------------------------------------------
    linkage_matrix = linkage(dist_condensed, method='average')

    # Choose optimal number of clusters using silhouette score
    best_score = -1
    best_k = 0

    for k in range(2, 10):
        clusters_k = fcluster(linkage_matrix, t=k, criterion='maxclust')
        score = silhouette_score(haplotypes, clusters_k, metric='hamming')

        if score > best_score:
            best_score = score
            best_k = k

    print(f"Best number of clusters based on silhouette score: {best_k}")

    clusters = fcluster(linkage_matrix, t=best_k, criterion='maxclust')

    # -----------------------------------------------------------------
    # PCA visualization
    # -----------------------------------------------------------------
    sample_names_clean = [name.replace("male_", "") for name in male_samples]

    pca = PCA(n_components=2)
    pca_coords = pca.fit_transform(haplotypes)

    pca_df = pd.DataFrame({
        'Sample': sample_names_clean,
        'PC1': pca_coords[:, 0],
        'PC2': pca_coords[:, 1],
        'Cluster': clusters.astype(str)
    })

    pca_df.to_csv("pca_clusters_includingMissing.csv", index=False)

    pca_plot = (
        ggplot(pca_df, aes(x='PC1', y='PC2', color='Cluster', fill='Cluster'))
        + geom_point(size=3, alpha=0.8)
        + stat_ellipse(geom='polygon', alpha=0.2, show_legend=False)
        + theme_classic()
        + labs(title='PCA of Y-linked Haplotypes (Males Only)')
        + theme(figure_size=(8, 6), legend_position='right')
    )

    pca_plot.save("Y_Haplotype_PCA.pdf", dpi=300)
    pca_plot.save("Y_Haplotype_PCA.svg", dpi=300)

    # -----------------------------------------------------------------
    # MDS visualization (based on Hamming distances)
    # -----------------------------------------------------------------
    dist_matrix = squareform(dist_condensed)
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    mds_coords = mds.fit_transform(dist_matrix)

    mds_df = pd.DataFrame({
        'MDS1': mds_coords[:, 0],
        'MDS2': mds_coords[:, 1],
        'Haplogroup': clusters.astype(str)
    })

    colors = {
        '1': '#ef3b2c',
        '2': 'green',
        '3': '#67a9cf',
        '4': '#756bb1'
    }

    mds_plot = (
        ggplot(mds_df, aes(x='MDS1', y='MDS2',
                           color='Haplogroup', fill='Haplogroup'))
        + geom_point(size=3, alpha=0.8)
        + stat_ellipse(geom='polygon', alpha=0.2)
        + scale_color_manual(values=colors)
        + scale_fill_manual(values=colors)
        + theme_classic()
        + theme(figure_size=(6, 4), legend_position='right')
    )

    mds_plot.save("Y_Haplotype_MDS_v2.pdf", dpi=300)
    mds_plot.save("Y_Haplotype_MDS_v2.svg", dpi=300)

    # -----------------------------------------------------------------
    # Identify shared SNPs among clusters and plot UpSet diagram
    # -----------------------------------------------------------------
    cluster_ids = np.unique(clusters)
    cluster_snps = {}

    for cid in cluster_ids:
        idx = np.where(clusters == cid)[0]
        snps_in_cluster = set(
            filtered_positions[np.any(haplotypes[idx, :], axis=0)]
        )
        cluster_snps[f'HG{cid}'] = snps_in_cluster

    upset_data = from_contents(cluster_snps)

    upset = UpSet(upset_data, subset_size='count', show_counts=True)
    upset.plot()

    plt.title("Shared SNPs Among Y-linked Haplogroup Clusters")
    plt.savefig("venndiagram.pdf")
    plt.show()

main()
