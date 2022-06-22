#!/bin/bash

mkdir stage_two && cd stage_two
mkdir raw_data && cd raw_data

#download normal and tumor datasets
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz

#download reference sequence
wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz

#unzip the reference sequence
gunzip hg19.chr5_12_17.fa.gz

#implement fastqc
mkdir qc_reports

#implement fastqc for all raw data
fastqc *.fastq.gz -o qc_reports/

#Implement multiQC (assembling quality control reports)
multiqc qc_reports -o qc_reports

#fastp

nano trim.sh

#!/bin/bash
mkdir trimmed_reads && cd trimmed_reads
mkdir qc_results

SAMPLES=(
  "SLGFSK-N_231335"
  "SLGFSK-N_231335"
  "SLGFSK-T_231336"
  "SLGFSK-T_231336"
)

for SAMPLE in "${SAMPLES[@]}"; do

  fastp \
    -i "$PWD/${SAMPLE}_r1_chr5_12_17.fastq.gz" \
    -I "$PWD/${SAMPLE}_r2_chr5_12_17.fastq.gz" \
    -o "trimmed_reads/${SAMPLE}_r1_chr5_12_17.fastq.gz" \
    -O "trimmed_reads/${SAMPLE}_r2_chr5_12_17.fastq.gz" \
    --html "trimmed_reads/${SAMPLE}_fastp.html"
done


#fastqc of trimmed_reads
cd \trimmed_reads
mkdir trimmed_qc_reports
fastqc *.fastq.gz -o trimmed_qc_reports/

#multiqc of trimmed_reads
multiqc trimmed_qc_reports -o trimmed_qc_reports/


#Read mapping
#copy reference file to reference directory

mkdir reference
cp hg19.chr5_12_17.fa reference

cd \reference

#Index reference file	
bwa index hg19.chr5_12_17.fa 

mkdir mapping

samtools view -b > mapping/SLGFSK-N_231335.bam
samtools view -b > mapping/SLGFSK-T_231336.bam

cd raw_data

#Perform alignment
bwa mem -R '@RG\tID:231335\tSM:Normal' reference/hg19.chr5_12_17.fa trimmed_reads/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz  trimmed_reads/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz > mapping/SLGFSK-N_231335.sam

bwa mem -R '@RG\tID:231336\tSM:Tumor' reference/hg19.chr5_12_17.fa trimmed_reads/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz trimmed_reads/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz > mapping/SLGFSK-T_231336.sam


cd mapping

#convert sam to bam
#samtools view -S -b filename.sam > filename.bam
samtools view -S -b SLGFSK-N_231335.sam > SLGFSK-N_231335.bam
samtools view -S -b SLGFSK-T_231336.sam > SLGFSK-T_231336.bam

#sort the bam files
samtools sort SLGFSK-N_231335.bam -o SLGFSK-N_231335.sorted.bam
samtools sort SLGFSK-T_231336.bam -o SLGFSK-T_231336.sorted.bam

#index sorted bam files
samtools index SLGFSK-N_231335.sorted.bam
samtools index SLGFSK-T_231336.sorted.bam

#mapped reads filtering
samtools view -q 1 -f 0x2 -F 0x8 -b SLGFSK-N_231335.sorted.bam > SLGFSK-N_231335.filtered1.bam
samtools view -q 1 -f 0x2 -F 0x8 -b SLGFSK-T_231336.sorted.bam > SLGFSK-T_231336.filtered1.bam

#To view the output of the results use 
samtools flagstat <bam file>

#duplicates removal
#collate-shuffle and group alignments by name
samtools collate SLGFSK-N_231335.filtered1.bam SLGFSK-N_231335.namecollate 
samtools collate SLGFSK-T_231336.filtered1.bam SLGFSK-T_231336.namecollate

#got an error running fixmate: Coordinate sorted, require grouped/sorted by queryname.
#so sorting the file first
samtools sort -n -o SLGFSK-N_231335.namecollate.bam SLGFSK-N_231335.namecollate.bam
samtools sort -n -o SLGFSK-T_231336.namecollate.bam SLGFSK-T_231336.namecollate.bam

samtools fixmate -m SLGFSK-N_231335.namecollate.bam SLGFSK-N_231335.fixmate.bam
samtools fixmate -m SLGFSK-T_231336.namecollate.bam SLGFSK-T_231336.fixmate.bam

samtools sort -@ 32 SLGFSK-N_231335.fixmate.bam -o SLGFSK-N_231335.positionsort.bam
samtools sort -@ 32 SLGFSK-T_231336.fixmate.bam -o SLGFSK-T_231336.positionsort.bam
 
samtools markdup -@32 -r SLGFSK-N_231335.positionsort.bam SLGFSK-N_231335.clean.bam
samtools markdup -@32 -r SLGFSK-T_231336.positionsort.bam SLGFSK-T_231336.clean.bam

#or
samtools rmdup SLGFSK-N_231335.sorted.bam  SLGFSK35.rdup 
samtools rmdup SLGFSK-T_231336.sorted.bam  SLGFSK36.rdup

cd ..

#left align bam
cat mapping/SLGFSK-N_231335.clean.bam  | bamleftalign -f reference/hg19.chr5_12_17.fa -m 5 -c > mapping/SLGFSK-N_231335.leftAlign.bam
cat mapping/SLGFSK-T_231336.clean.bam  | bamleftalign -f reference/hg19.chr5_12_17.fa -m 5 -c > mapping/SLGFSK-T_231336.leftAlign.bam

#Recalibrate read mapping qualities
samtools calmd -@ 32 -b mapping/SLGFSK-N_231335.leftAlign.bam reference/hg19.chr5_12_17.fa > mapping/SLGFSK-N_231335.recalibrate.bam
samtools calmd -@ 32 -b mapping/SLGFSK-T_231336.leftAlign.bam reference/hg19.chr5_12_17.fa > mapping/SLGFSK-T_231336.recalibrate.bam

#Refilter read mapping qualities
bamtools filter -in mapping/SLGFSK-N_231335.recalibrate.bam -mapQuality "<=254" > mapping/SLGFSK-N_231335.refilter.bam
bamtools filter -in mapping/SLGFSK-T_231336.recalibrate.bam -mapQuality "<=254" > mapping/SLGFSK-T_231336.refilter.bam

#Variant calling and classification 
#Convert data to pileup
mkdir variants
samtools mpileup -f reference/hg19.chr5_12_17.fa mapping/SLGFSK-N_231335.refilter.bam --min-MQ 1 --min-BQ 28 > variants/SLGFSK-N_231335.pileup
samtools mpileup -f reference/hg19.chr5_12_17.fa mapping/SLGFSK-T_231336.refilter.bam --min-MQ 1 --min-BQ 28 > variants/SLGFSK-T_231336.pileup

#copy VarScan.v2.3.9.jar from sample_tut to darwin
cp sample_tut/VarScan.v2.3.9.jar darwin
cp VarScan.v2.3.9.jar stage_two/raw_data

#Call variants
java -jar VarScan.v2.3.9.jar somatic variants/SLGFSK-N_231335.pileup variants/SLGFSK-T_231336.pileup variants/SLGFSK \ --normal-purity 1  --tumor-purity 0.5 --output-vcf 1

#merge vcf
bgzip variants/SLGFSK.snp.vcf > variants/SLGFSK.snp.vcf.gz
bgzip variants/SLGFSK.indel.vcf > variants/SLGFSK.indel.vcf.gz
tabix variants/SLGFSK.snp.vcf.gz
tabix variants/SLGFSK.indel.vcf.gz
cd variants
bcftools merge SLGFSK.snp.vcf.gz SLGFSK.indel.vcf.gz -o SLGFSK.vcf --force -samples




