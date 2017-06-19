# vcf2fasta  
produces a fasta file from a VCF  
**compilation**  
g++ vcf2fasta.cpp -std=c++17 -O3 -o vcf2fasta  
  
**exemple**:  
./vcf2fasta subVCF_ama_11.vcf 10 Hmel201012 ama Hmel201012.fasta  
arg1: vcf file  
arg2: minimum number of reads to call one allele at one position (needs 2.cov for being homozygote for instance)  
arg3: name of the loci/contig/chromosome to print in sequences's ID of the final fasta file  
arg4: species name to print in sequence's ID  
arg5: name of the final fasta file  

