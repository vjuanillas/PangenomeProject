# Sample flow and commands:
## STRUCTURAL VARIATIONS

1. Call structural variations using gfatools

~/tools/gfatools/./gfatools bubble asm5.nip.gfa > asm5.nip.bubble.tsv

2. Divide the output to useable files for next step:
2.1 Biallelic bubbles can be divided based on the number of nodes (3 or 4) and paths (2). Otherwise, these are multiallelic bubbles

awk '$5 == 2 {{ print $1,$2,$4,$5,$12 }}' asm5.nip.bubble.tsv > asm5.nip.biallelic.bubble.tsv

awk '$5>2 && $5 < 8 {{ print $1,$2,$4,$5,$12 }}' asm5.nip.bubble.tsv > asm5.nip.multiallelic.bubble.tsv 

awk '{{ print $1, $2, $2 + 1, $1"_"$2 }}' OFS="\t" asm5.nip.bubble.tsv > asm5.nip.bubble.bed

3. Use 05.1_get_bialsv.py to annotate the biallelic structural variations. 
python 05_get_bialsv.py -b biallelic.bubble.tsv -c combined_coverage.tsv > asm5.nip.biallelic_sv.tsv

4. Trace the paths through the bubbles per assembly using minigraph --call

for asm_file in IRGSP-1.0_genome.fasta nh232.renamed.fasta cw02.renamed.fasta 

nh236.renamed.fasta nh286.renamed.fasta nh273.renamed.fasta; do
  case $asm_file in
    IRGSP-1.0_genome.fasta)
      base_name="IRGSP"
      ;;
    nh232.renamed.fasta)
      base_name="nh232"
      ;;
    cw02.renamed.fasta)
      base_name="cw02"
      ;;
    nh236.renamed.fasta)
      base_name="nh236"
      ;;
    nh286.renamed.fasta)
      base_name="nh286"
      ;;
    nh273.renamed.fasta)
      base_name="nh273"
      ;;
  esac
  
  ~/tools/minigraph.0.21-r606/./minigraph -cxasm --call /opt/home/jong/oryza/asm/renamed_asm/00_mash_distance/asm5.nip.gfa /opt/home/jong/oryza/asm/renamed_asm/$asm_file > ${base_name}.alleles.bed
done

5. Parse the call output using 05.2_trace_path.py and trace each allele back to the bubbles.2_trace_path
python 05.2_trace_path.py -b /Users/jongpaduhilao/Desktop/LAB_Files/Initial_Pangenome_analysis/Trial_4/crysnanto_bubble/asm5.nip.biallelic_sv.tsv -a /Users/jongpaduhilao/Desktop/LAB_Files/Initial_Pangenome_analysis/bubble/04_bubble -o test_05_biallelic_asm_paths

6. Intersect the structural variations back to the reference annotation using bedtools intersect
bedtools intersect -a 05x4_biallelic_asm_paths.wID.bed -b ../IRGSP-1.0_representative/transcripts.gff -wb -loj > 05x4_biallelic_asm2annot.intersect.tsv

7. Parse the intersected output to see the estimate locations of the structural variations 
python Scripts/06.1_annot_sv.py -i Raw_data/data/05x4_biallelic_asm2annot.intersect.tsv -o 06_biallelic_sv_annot.tsv
