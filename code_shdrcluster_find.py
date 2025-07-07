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


def rdsr_cluster_find(infile, outfile_all, outfile_without_clusters, log_file):
    blastres = defaultdict(str)
    all_genomes = set()
    genomes_with_clusters = set()

    # 解析 BLAST 结果
    with open(infile, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            query, subj = parts[1], parts[0]
            genome, _ = subj.split('___', 1)
            all_genomes.add(genome)

            if subj in blastres:
                if query not in blastres[subj]:
                    blastres[subj] += query + ","
            else:
                blastres[subj] = query + ","

    classified = {subj: queries.rstrip(',').split(',') for subj, queries in blastres.items()}
    operon = defaultdict(lambda: defaultdict(list))

    for subject in classified:
        genome, locus = subject.split('___', 1)
        index = locus.split('_')[-1]
        contig = locus.rsplit('_', 1)[0]  # 去除最后的 _index
        operon[genome][contig].append(int(index))

    total_clusters = 0
    qualifying_clusters = 0

    with open(outfile_all, 'w') as fout_all, \
            open(outfile_without_clusters, 'w') as fout_without, \
            open(log_file, 'w') as log:  # 修改为 'w' 模式，避免追加

        for genome in operon:
            for contig in operon[genome]:
                clusters = cluster(operon[genome][contig], gap)
                for c in clusters:
                    genes_in_cluster = []
                    for idx in c:
                        key = f"{genome}___{contig}_{idx}"
                        genes_in_cluster.extend(classified[key])

                    # 判断基因簇是否符合要求
                    # 必须包含 sHdrA
                    has_sHdrA = 'sHdrA' in genes_in_cluster

                    # 必须包含 sHdrB1 或 sHdrB2 中的一个
                    has_sHdrB = any(g in {'sHdrB1', 'sHdrB2'} for g in genes_in_cluster)

                    # 必须包含 sHdrC1 或 sHdrC2 中的一个
                    has_sHdrC = any(g in {'sHdrC1', 'sHdrC2'} for g in genes_in_cluster)

                    if has_sHdrA and has_sHdrB and has_sHdrC:
                        total_clusters += 1
                        genomes_with_clusters.add(genome)
                        acc = genome.replace('_protein.faa', '')

                        for idx in c:
                            key = f"{genome}___{contig}_{idx}"
                            line = f"{total_clusters}\t{acc}\t{contig}_{idx}\t{key}\t{','.join(classified[key])}\n"
                            fout_all.write(line)

                        qualifying_clusters += 1

        # 记录没有 cluster 的基因组
        genomes_without_clusters = all_genomes - genomes_with_clusters
        for genome in genomes_without_clusters:
            fout_without.write(f"{genome}\n")

        log.write(f"Total clusters: {total_clusters}\n")
        log.write(f"Clusters with sHdrA + (sHdrB1 or sHdrB2) + (sHdrC1 or sHdrC2): {qualifying_clusters}\n")
        log.write(f"Genomes without clusters: {len(genomes_without_clusters)}\n")


def parse_genome_for_multiple_clusters(outfile_all, outfile_multiple, log_file):
    acc_cluster_map = defaultdict(set)
    cluster_gene_map = defaultdict(lambda: defaultdict(list))

    with open(outfile_all, 'r') as fin:
        for line in fin:
            parts = line.strip().split('\t')
            cluster_id, acc = parts[0], parts[1]
            acc_cluster_map[acc].add(cluster_id)
            cluster_gene_map[cluster_id][acc].append(parts[-1])

    with open(outfile_multiple, 'w') as fout, open(log_file, 'a') as log:
        multiple_cluster_count = 0
        for acc, clusters in acc_cluster_map.items():
            if len(clusters) > 1:
                multiple_cluster_count += 1
                compositions = [f"{cid}: {', '.join(set(cluster_gene_map[cid][acc]))}" for cid in clusters]
                fout.write(f"{acc}\t{' | '.join(compositions)}\n")

        log.write(f"Genomes with multiple clusters: {multiple_cluster_count}\n")


import sys

def main():
    # 参数配置
    gap = 3  # 允许的基因间隔
    
    # 方法1: 使用命令行参数
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file>")
        print("Example: python script.py /path/to/your/input.txt")
        sys.exit(1)
    
    infile = sys.argv[1]
    
    # 方法2: 如果想要交互式输入，可以替换上面的代码为：
    # infile = input("Please enter the input file path: ")
    
    # 检查文件是否存在
    if not os.path.exists(infile):
        print(f"Error: File '{infile}' not found!")
        sys.exit(1)
    
    basename = os.path.splitext(infile)[0]
    
    outfile_all = basename + ".operon.all.txt"
    outfile_without = basename + ".operon.without.txt"
    log_file = basename + ".overview.txt"
    outfile_multiple = basename + ".multiple_clusters.txt"
    
    print(f"Processing file: {infile}")
    print(f"Output files will be saved with prefix: {basename}")
    
    # 运行主函数
    rdsr_cluster_find(infile, outfile_all, outfile_without, log_file)
    parse_genome_for_multiple_clusters(outfile_all, outfile_multiple, log_file)
    
    print("Analysis completed!")

if __name__ == "__main__":
    main()