#!/usr/bin/env python3
import os
import subprocess
from ete3 import Tree
import time
import argparse

def read_processed_a(output_file):
    processed = set()
    if os.path.exists(output_file):
        with open(output_file, 'r') as f:
            # 跳过标题行
            header = f.readline()
            # 验证标题是否正确（可选）
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split('\t')
                    try:
                        a_val = int(parts[0])
                        processed.add(a_val)
                    except (IndexError, ValueError):
                        continue  # 忽略格式错误的行
    return processed

def main():
    parser = argparse.ArgumentParser(description="Extract MRCA node labels from chronogram trees.")
    parser.add_argument("chain_file", type=str, help="Path to a single chain file name (without extension)")
    parser.add_argument("--start", "-s", type=int, default=250, help="Start value for a (default: 250)")
    parser.add_argument("--end", "-e", type=int, default=600, help="End value for a (inclusive, default: 600)")
    parser.add_argument("--output", "-o", type=str, default="mrca_results.txt", help="Output file name (default: mrca_results.txt)")
    args = parser.parse_args()

    chain = args.chain_file
    start_a = args.start
    end_a = args.end
    output_file = args.output
    print(f"chain={chain}, start_a={start_a}, end_a={end_a}, output_file={output_file}")

    # 读取已处理的a值
    processed_as = read_processed_a(output_file)
    print(f"已处理的a值：{sorted(processed_as)}")

    # 确保输出文件存在并写入标题（如果不存在）
    if not os.path.exists(output_file):
        with open(output_file, "w") as f:
            f.write("a\tmrca_label_tsox\tmrca_label_tsox_rdsr\tmrca_label_tsox_soxCD\n")

    # 遍历所有a值，跳过已处理的
    for a in range(start_a, end_a + 1):
        if a in processed_as:
            print(f"跳过已处理的a值：{a}")
            continue

        b = a + 1
        cmd = f"readdiv -x {a} 1 {b} {chain}"
        print(f"执行命令: {cmd}")
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"执行readdiv命令时出错: {e}")
            continue
        
        chronogram_file = f"{chain}_sample.chronogram"
        wait_time = 0
        while not os.path.exists(chronogram_file) and wait_time < 30:
            time.sleep(1)
            wait_time += 1
        
        if not os.path.exists(chronogram_file):
            print(f"错误: {chronogram_file} 文件未生成")
            continue
        
        try:
            tree = Tree(chronogram_file, format=1)
            tsox_leaves = tree.search_nodes(name="acc336") + tree.search_nodes(name="acc158")
            tsox_rdsr_leaves = tree.search_nodes(name="acc259") + tree.search_nodes(name="acc374")
            tsox_soxCD_leaves = tree.search_nodes(name="acc210") + tree.search_nodes(name="acc412")
            
            if not tsox_leaves or not tsox_rdsr_leaves or not tsox_soxCD_leaves:
                print(f"警告: 在{chronogram_file}中未找到某些叶子节点")
                continue
            
            mrca_tsox = tree.get_common_ancestor(tsox_leaves[0], tsox_leaves[1])
            mrca_tsox_rdsr = tree.get_common_ancestor(tsox_rdsr_leaves[0], tsox_rdsr_leaves[1])
            mrca_tsox_soxCD = tree.get_common_ancestor(tsox_soxCD_leaves[0], tsox_soxCD_leaves[1])
            
            mrca_label_tsox = mrca_tsox.name
            mrca_label_tsox_rdsr = mrca_tsox_rdsr.name
            mrca_label_tsox_soxCD = mrca_tsox_soxCD.name
            
            with open(output_file, "a") as f:
                f.write(f"{a}\t{mrca_label_tsox}\t{mrca_label_tsox_rdsr}\t{mrca_label_tsox_soxCD}\n")
            
            print(f"a = {a}: MRCA标签已写入{output_file}")
        
        except Exception as e:
            print(f"处理{chronogram_file}时出错: {e}")
            continue

    print(f"所有操作完成，结果保存在 {output_file}")

if __name__ == "__main__":
    main()