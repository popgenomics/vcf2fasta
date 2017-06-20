# vcf2fasta  
produces a fasta file from a VCF  
**compilation**  
g++ vcf2fasta.cpp -std=c++17 -O3 -o vcf2fasta  
  
**exemple**:  
./vcf2fasta subVCF_ama_11.vcf 10 ama  
**arg1:** vcf file  
**arg2:** minimum number of reads to call one allele at one position (needs 2.cov for being homozygote for instance)  
**arg3:** species name to print in sequence's ID  

