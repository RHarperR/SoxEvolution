
INIT_FILE=${1:-05_initial.txt}   # 初始文件
JOBS=${2:-8}                     # 第二个参数可指定并发核数，默认 8

# 生成所有 (D,L) 组合，用制表符分隔
for ((D=400; D<=1000; D+=100)); do
  for ((L=200; L>=25; L-=25)); do
    printf "%s\t%d\t%d\n" "$INIT_FILE" "$D" "$L"
  done
done | \
parallel --colsep '\t' -j "$JOBS" \
         'echo "== Running D={2} L={3} ==" && OptRoot.linux -i {1} -D {2} -T 600 -L {3} -o costfind.{2}.600.{3}'