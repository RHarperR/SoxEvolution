import re
from collections import defaultdict
import os

def cluster(l, gap):
    l.sort()
    clusters = []
    cur = []
    for i in l:
        if not cur or i - cur[-1] > gap:
            if cur:
                clusters.append(cur)
            cur = [i]
        else:
            cur.append(i)
    if cur:
        clusters.append(cur)
    return clusters

def sox_cluster_find(infile, outfile_all, outfile_soxBAXYZ, outfile_soxBAXYZCD, outfile_without_clusters, log_file):
    blastres = defaultdict(str)
    all_genomes = set()
    genomes_with_clusters = set()
    for l in open(infile):
        a = re.split(r'\t', l.strip())
        query, subj = a[1], a[0]
        genome, _ = re.split(r'___', subj)
        all_genomes.add(genome)
        if subj in blastres:
            if query not in blastres[subj]:
                blastres[subj] += query.replace('sox', '')
        else:
            blastres[subj] = query
    classified = {subj: queries for subj, queries in blastres.items()}
    operon = defaultdict(lambda: defaultdict(list))
    for i in classified:
        genome, locus = re.split(r'___', i)
        index = re.split(r'_', locus)[-1]
        contig = re.sub(r'_%s$' % index, r'', locus)
        operon[genome][contig].append(int(index))
    total_clusters = 0
    soxBAXYZ_clusters = 0
    soxBAXYZCD_clusters = 0
    with open(outfile_all, 'w') as fout_all, open(outfile_soxBAXYZ, 'w') as fout_soxBAXYZ, open(outfile_soxBAXYZCD, 'w') as fout_soxBAXYZCD, open(outfile_without_clusters, 'w') as fout_without, open(log_file, 'a') as log:
        for g in operon:
            for contig in operon[g]:
                clusters = cluster(operon[g][contig], gap)
                for c in clusters:
                    genes_in_cluster = [classified[f'{g}___{contig}_{m}'] for m in c]
                    if 'soxB' in genes_in_cluster and any(
                            gene in genes_in_cluster for gene in ['soxA', 'soxX', 'soxAX']) and any(
                            gene in genes_in_cluster for gene in ['soxY', 'soxZ', 'soxYZ']):
                        acc = g.replace('_protein.faa', '')
                        total_clusters += 1
                        genomes_with_clusters.add(g)  # 记录含有完整 cluster 的基因组
                        for m in c:
                            fout_all.write(
                                f'{total_clusters}\t{acc}\t{contig}_{m}\t{g}___{contig}_{m}\t{classified[g + "___" + contig + "_" + str(m)]}\n')
                        if 'soxY' in genes_in_cluster and 'soxZ' in genes_in_cluster and (
                                'soxA' in genes_in_cluster or 'soxAX' in genes_in_cluster):
                            soxBAXYZ_clusters += 1
                            for m in c:
                                fout_soxBAXYZ.write(
                                    f'{total_clusters}\t{acc}\t{contig}_{m}\t{g}___{contig}_{m}\t{classified[g + "___" + contig + "_" + str(m)]}\n')
                            if 'soxC' in genes_in_cluster and 'soxD' in genes_in_cluster:
                                soxBAXYZCD_clusters += 1
                                for m in c:
                                    fout_soxBAXYZCD.write(
                                        f'{total_clusters}\t{acc}\t{contig}_{m}\t{g}___{contig}_{m}\t{classified[g + "___" + contig + "_" + str(m)]}\n')
        genomes_without_clusters = all_genomes - genomes_with_clusters
        for genome in genomes_without_clusters:
            fout_without.write(f"{genome}\n")
        log.write(f"Total clusters: {total_clusters}\n")
        log.write(f"Clusters with soxBAXYZ: {soxBAXYZ_clusters}\n")
        log.write(f"Clusters with soxBAXYZCD: {soxBAXYZCD_clusters}\n")
        log.write(f"Genomes without sox clusters: {len(genomes_without_clusters)}\n")

def parse_genome_for_multiple_clusters(outfile1, outfile_multiple, log_file):
    acc_cluster_map = defaultdict(set)
    cluster_gene_map = defaultdict(lambda: defaultdict(list))
    with open(outfile1, 'r') as fin:
        for line in fin:
            cluster, acc, _, _, genes = line.strip().split('\t')
            acc_cluster_map[acc].add(cluster)
            cluster_gene_map[cluster][acc].append(genes)  # Collect genes for each cluster
    # Write to output file for accessions with multiple clusters and include gene compositions
    with open(outfile_multiple, 'w') as fout, open(log_file, 'a') as log:
        for acc, clusters in acc_cluster_map.items():
            if len(clusters) > 1:  # Only keep accessions with multiple clusters
                gene_compositions = []
                for cluster in clusters:
                    genes = cluster_gene_map[cluster][acc]  # Get genes for this cluster
                    gene_compositions.append(f"{cluster}: {', '.join(genes)}")
                fout.write(f'{acc}\t{" | ".join(gene_compositions)}\n')
        log.write(f"Number of acc with multiple clusters: {len([acc for acc in acc_cluster_map if len(acc_cluster_map[acc]) > 1])}\n")

# 参数设置
gap = 4
infile = "00_sox_all.txt"
basename = os.path.splitext(infile)[0]
outfile1 = re.sub(r'(\d+)_', r'01_', basename) + ".operon.all.txt"  # 所有潜在执行功能的 sox cluster
outfile2 = re.sub(r'(\d+)_', r'01_', basename) + ".operon.soxBAXYZ.txt"  # 所有有完整 soxBA(X)YZ 的 sox cluster
outfile3 = re.sub(r'(\d+)_', r'01_', basename) + ".operon.soxBAXYZCD.txt"  # 所有有完整 soxBA(X)YZCD 的 sox cluster
outfile_multiple = re.sub(r'(\d+)_', r'01_', basename) + ".operon.all.multiple.txt"  # 只包含 multiple sox cluster 的情况
outfile_without_clusters = re.sub(r'(\d+)_', r'01_', basename) + ".operon.without_clusters.txt"  # 不含有完整 sox cluster 的 genome
log_file = re.sub(r'(\d+)_', r'01_', basename) + ".overview.txt"  # 日志文件路径

# 执行函数并记录输出
sox_cluster_find(infile, outfile1, outfile2, outfile3, outfile_without_clusters, log_file)
parse_genome_for_multiple_clusters(outfile1, outfile_multiple, log_file)
