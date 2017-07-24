# vcf2fasta  
produces a fasta file from a VCF  
  
## compilation  
g++ vcf2fasta.cpp -std=c++17 -O3 -o vcf2fasta  
  
## example:  
### v1:  
./vcf2fasta subVCF_ama_11.vcf 8 ama  
**arg1:** vcf file  
**arg2:** minimum number of reads to call one allele at one position (needs 2.cov for being homozygote for instance)  
**arg3:** species name to print in sequence's ID  
  
### v2:  
./vcf2fasta subVCF_ama_11.vcf 8 30 ama  
**arg1:** vcf file  
**arg2:** minimum number of reads to call one allele at one position   
**arg3:** minimum phred quality score to call a genotype  
**arg4:** species name to print in sequence's ID  

