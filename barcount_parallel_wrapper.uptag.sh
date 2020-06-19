num_procs=$1
num_jobs="\j"  # The prompt escape for number of jobs currently running
for f in */*up*.assembled.fastq; do
  while (( ${num_jobs@P} >= num_procs )); do
    wait -n
  done
  barcount --fastq "$f" --rmdup --flanking_left CAAGCTAAGATATC --flanking_right TTTAAATGCGAAGTAA --max_distance_flanks 1 --max_distance_barcode 3 --barcode_table up.Safe.3060.csv --verbose --out "${f}_Barcode_filter" > ${f}_Barcode_filter1.log 2>&1 & 
done