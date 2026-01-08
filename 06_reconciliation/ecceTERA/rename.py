#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import argparse
from Bio import Phylo
from io import StringIO

def process_leaf_names(tree1_file: str, tree2_file: str, output_file: str):
    # ---------- 1. 读取 tree1 并建立映射 ----------
    tree1 = Phylo.read(tree1_file, "newick")
    tree1_map = {}
    for leaf in tree1.get_terminals():
        name = leaf.name
        if name and len(name) >= 16:
            tree1_map[name[:16]] = name
    print(f"从第一棵树构建完成，共 {len(tree1_map)} 条前 16 字符映射")

    # ---------- 2. 逐棵处理 tree2 ----------
    def cleaned_lookup(raw_name: str):
        if not raw_name:
            return raw_name
            
        processed = raw_name
        
        # 保留去除前缀的步骤
        if processed.startswith("GB_") or processed.startswith("RS_"):
            processed = processed[3:]
            
        # 提取前16个字符（用于匹配）和剩余部分
        match_part = processed[:15]
        remaining_part = processed[15:]
        
        # 将前16个字符中的下划线替换为两个点
        lookup_key = match_part.replace('_', '..')
        
        # 如果找到匹配项，使用tree1中的完整名称加上剩余部分
        if lookup_key in tree1_map:
            return tree1_map[lookup_key] + remaining_part
        else:
            # 如果没有找到匹配，保留原始名称
            return raw_name

    replaced_total = 0
    with open(output_file, "w") as out_fh:
        for tree_idx, tree2 in enumerate(Phylo.parse(tree2_file, "newick"), 1):
            replaced_count = 0
            for leaf in tree2.get_terminals():
                new_name = cleaned_lookup(leaf.name)
                if new_name != leaf.name:
                    leaf.name = new_name
                    replaced_count += 1
            # 把替换后的树写回文件（Newick 单行）
            out_fh.write(tree2.format("newick").strip() + "\n")
            replaced_total += replaced_count
            print(f"处理第 {tree_idx} 棵树，替换了 {replaced_count} 个叶子名称")
    
    print(f"全部处理完成，共替换 {replaced_total} 个叶子名称，结果已写入 {output_file}")

def main():
    """主函数，处理命令行参数"""
    parser = argparse.ArgumentParser(
        description="根据参考树中的完整名称重命名目标树中的叶子节点名称",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python script.py -r reference_tree.newick -t target_trees.newick -o output.newick
        """
    )
    
    parser.add_argument(
        '-r', '--reference',
        required=True,
        help='参考树文件路径（包含完整的标准名称）'
    )
    
    parser.add_argument(
        '-t', '--target',
        required=True,
        help='需要重命名的目标树文件路径'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='输出文件路径'
    )
    
    args = parser.parse_args()
    
    print(f"参考树文件: {args.reference}")
    print(f"目标树文件: {args.target}")
    print(f"输出文件: {args.output}")
    print("-" * 50)
    
    # 运行处理函数
    process_leaf_names(
        tree1_file=args.reference,
        tree2_file=args.target,
        output_file=args.output
    )

if __name__ == "__main__":
    main()