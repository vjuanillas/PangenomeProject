cd ~/asms/

# try two parameters of -s
mash triangle all_rice_O_mer.fa.gz -i > all_rice_O_mer.mash_triangle.default.txt

mash triangle all_rice_O_mer.fa.gz -s 10000 -i > all_rice_O_mer.mash_triangle.10k.txt

sed 1,1d all_rice_O_mer.mash_triangle.default.txt | tr '\t' '\n' | grep chr -v | LC_ALL=C sort -g -k 1nr | uniq | head -n 2

sed 1,1d all_rice_O_mer.mash_triangle.10k.txt | tr '\t' '\n' | grep chr -v | LC_ALL=C sort -g -k 1nr | uniq | head -n 2

# To Estimate the segment length, get the longets repeat length
sed 1,3d asm5.npb.fa.out | awk '{print $6, $7, $7-$6}' | awk '$3 > max {max = $3} END {print max}' | head

echo "--Check longest repeat length--"

# Partition the sequences into 12 communities/chromosomes
partition-before-pggb -i all_rice_O_mer.fa.gz -o output_all_O_mer_p90 -n 7 -t 16 -p 90 -s 50000 -V 'O_mer:100000'

# Run pggb
# To edit so we can change the parameters for p=[90-60] and k=[19-23]
/usr/bin/time -v bash -c '
seq 0 11 | while read i; do
pggb -i output_all_O_mer_p90/all_rice_O_mer.fa.gz.dac1d73.community.$i.fa \
	     -o output_all_O_mer_p90/all_rice_O_mer.fa.gz.dac1d73.community.$i.fa.out \
	     -s 50000 -l 250000 -p 90 -c 1 -K 19 -F 0.001 -g 30 \
	     -k 23 -f 0 -B 10M \
	     -n 7 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ \
	     -Y "#" -V O_mer:100000 --threads 16 --poa-threads 16
done
'


