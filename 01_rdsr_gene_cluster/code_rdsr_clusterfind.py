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


def rdsr_cluster_find(infile, outfile_all, outfile_rdsrAB_dsrEFH, outfile_without_clusters,
                      log_file):
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
                blastres[subj] += query + ","
        else:
            blastres[subj] = query + ","

    classified = {subj: queries.rstrip(',').split(',') for subj, queries in blastres.items()}
    operon = defaultdict(lambda: defaultdict(list))

    for i in classified:
        genome, locus = re.split(r'___', i)
        index = re.split(r'_', locus)[-1]
        contig = re.sub(r'_%s$' % index, r'', locus)
        operon[genome][contig].append(int(index))

    total_clusters = 0
    rdsrAB_dsrEFH_clusters = 0

    with open(outfile_all, 'w') as fout_all, \
            open(outfile_rdsrAB_dsrEFH, 'w') as fout_rdsr, \
            open(outfile_without_clusters, 'w') as fout_without, \
            open(log_file, 'a') as log:

        for g in operon:
            for contig in operon[g]:
                clusters = cluster(operon[g][contig], gap)  # 使用固定gap值4
                for c in clusters:
                    genes_in_cluster = []
                    for m in c:
                        key = f"{g}___{contig}_{m}"
                        genes_in_cluster.extend(classified[key])

                    # 优化条件检查：单次遍历同时检查两类基因
                    has_rdsr = False
                    has_dsr = False
                    for gene in genes_in_cluster:
                        if gene in {'dsrA', 'dsrB'}:
                            has_rdsr = True
                        elif gene in {'dsrE', 'dsrF', 'dsrH'}:
                            has_dsr = True
                        if has_rdsr and has_dsr:
                            break  # 提前退出循环

                    if has_rdsr and has_dsr:
                        total_clusters += 1
                        genomes_with_clusters.add(g)
                        acc = g.replace('_protein.faa', '')

                        # 合并写入操作，减少循环次数
                        for m in c:
                            key = f"{g}___{contig}_{m}"
                            line = f"{total_clusters}\t{acc}\t{contig}_{m}\t{key}\t{','.join(classified[key])}\n"
                            fout_all.write(line)
                            fout_rdsr.write(line)  # 同时写入两个文件

                        rdsrAB_dsrEFH_clusters += 1

        # 记录没有cluster的基因组
        genomes_without_clusters = all_genomes - genomes_with_clusters
        for genome in genomes_without_clusters:
            fout_without.write(f"{genome}\n")

        log.write(f"Total clusters: {total_clusters}\n")
        log.write(f"Clusters with rdsrAB + dsrEFH: {rdsrAB_dsrEFH_clusters}\n")
        # 移除不存在的rdsrAB_dsrEFHCD_clusters统计
        log.write(f"Genomes without clusters: {len(genomes_without_clusters)}\n")


def parse_genome_for_multiple_clusters(outfile1, outfile_multiple, log_file):
    acc_cluster_map = defaultdict(set)
    cluster_gene_map = defaultdict(lambda: defaultdict(list))

    with open(outfile1, 'r') as fin:
        for line in fin:
            parts = line.strip().split('\t')
            cluster_id, acc = parts[0], parts[1]
            acc_cluster_map[acc].add(cluster_id)
            cluster_gene_map[cluster_id][acc].append(parts[-1])

    with open(outfile_multiple, 'w') as fout, open(log_file, 'a') as log:
        for acc, clusters in acc_cluster_map.items():
            if len(clusters) > 1:
                compositions = []
                for cid in clusters:
                    # 直接获取基因列表，避免不必要的sum操作
                    genes = cluster_gene_map[cid][acc]
                    compositions.append(f"{cid}: {', '.join(set(genes))}")
                fout.write(f"{acc}\t{' | '.join(compositions)}\n")
        # 修正括号问题
        log.write(
            f"Genomes with multiple clusters: {len([acc for acc in acc_cluster_map if len(acc_cluster_map[acc]) > 1])}\n")


# 参数配置（修正缩进问题）
gap = 10
infile = "E:/Research/10 soxB stepwise/data_new/25_dsr_all.txt"
basename = os.path.splitext(infile)[0]

outfile_all = basename + ".operon.all.txt"
outfile_rdsrAB_dsrEFH = basename + ".operon.rdsrAB_dsrEFH.txt"
outfile_without = basename + ".operon.without.txt"
log_file = basename + ".overview.txt"
outfile_multiple = basename + ".multiple_clusters.txt"

# 运行主函数
rdsr_cluster_find(infile, outfile_all, outfile_rdsrAB_dsrEFH, outfile_without, log_file)
parse_genome_for_multiple_clusters(outfile_all, outfile_multiple, log_file)
