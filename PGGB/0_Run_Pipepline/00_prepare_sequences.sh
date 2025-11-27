working_dir=pwd

ASM_DIR="~/asms/"
FASTIX="~/.cargo/bin/fastix"

cd $ASM_DIR

sed '/^>/ s/^>\([^ ]*\).*/>out_chr\1/' O_mer.fa > O_mer.renamed.fa

mv cw02.renamed.fasta Or.renamed.fasta
mv nh232.renamed.fasta Osj.renamed.fasta
mv nh236.renamed.fasta Osi.renamed.fasta
mv nh286.renamed.fasta Ob.renamed.fasta
mv nh273.renamed.fasta Og.renamed.fasta

# fastix
for asm in *.renamed.fasta; do
basename=$(basename "$asm" .renamed.fasta);
fastix -p "${basename}#0#" $asm > "${basename}_prefixed.fa"; 
done 

export PATH="~/.cargo/bin:$PATH"


fastix -p "O_mer#0#" $asm > "O_mer_prefixed.fa"
done

fastix -p "Osj(IRGSP)#0#" $asm > "IRGSP.prefixed.fa" 
done

# merge all fasta files
cat O_mer_prefixed.fa IRGSP_prefixed.fa Or_prefixed.fa Osj_prefixed.fa Osi_prefixed.fa Ob_prefixed.fa Og_prefixed.fa> all_rice_O_mer.fa

bgzip -@ 4 all_rice_O_mer.fa
samtools faidx all_rice_O_mer.fa.gz # had to do manual local dir install 
