# varibase v1.0
This is the code to correct variation errors of SNPs/indles from genome assembly

You have a variation VCF file from GATK or samtools, but you want to change the bases in your reference sequence.

### Download and Compile:

    $ git clone  https://github.com/wtsi-hpag/varibase.git 
    $ cd varibase 
    $ make 

 
### For SNP base change
egrep "DP=" var.flt.vcf | awk '((length($4)==1)&&(length($5)==1)){print $1,$2,$4,$5,$6}' > snp.dat
./varibase -snp 1 snp.dat input.fasta output.fasta

### For indel base change

#egrep "DP=" var.flt.vcf | egrep INDEL | awk '{print $1,$2,$4,$5,$6}' > indel.dat
./varibase -indel 1 indel.dat input.fasta output.fasta 

### Indel length can be upto 1000bps

Note: please always do SNP corrections first and then indels. 

Any problems, contact Zemin Ning ( zn1@sanger.ac.uk )
