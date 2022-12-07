#!/bin/bash

## Enter Data directory
cd /home/ccmb/Krishna_Projects/Test_GATK4  ##1
REF='/home/ccmb/Krishna_Projects/ref' ##2  ## The reference directory should contain all the reference files including the snp databases.
GATK4='/home/ccmb/Programs/gatk-4.3.0.0' ##3
BED='/home/ccmb/Krishna_Projects/truseq-exome-targeted-regions-manifest-v1-2.bed' ##4

## Reading files in folder
for i in $(ls *_R1_001.fastq.gz | cut -d "_" -f 1,2)   ##5
do
mkdir ${i}
echo "$i"
mv ${i}_R1_001.fastq.gz ${i}
mv ${i}_R2_001.fastq.gz ${i}
done

cd /home/ccmb/Krishna_Projects/Test_GATK4   ##6

ls | while read line
do 
echo $line

cd /home/ccmb/Krishna_Projects/Test_GATK4/$line/    ##7

echo "DATE"
date

## 1. FastqToSam Start
echo "############################################ FastqToSam Start ############################################"
java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar FastqToSam -F1 *_R1_001.fastq.gz -F2 *_R2_001.fastq.gz -O unmapped.bam -RG H0164.2 -SM ${line} -LB Solexa-272222 -PU H0164ALXX140820.2 -PL ILLUMINA
echo "############################################# FastqToSam End #############################################"

## 2. Adapter Marking Start
echo "############################################ Adapter Marking Start ############################################"  
java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar MarkIlluminaAdapters -I unmapped.bam -O unmapped_markilluminaadapters.bam -M unmapped_markilluminaadapters_metrics.txt
echo "############################################# Adapter Marking End #############################################"

## 3. Sam To Fastq Start
echo "############################################ Sam To Fastq Start ############################################" 
java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar SamToFastq -I unmapped_markilluminaadapters.bam --FASTQ unmapped_markilluminaadaptors_R1.fastq.gz --SECOND_END_FASTQ unmapped_markilluminaadaptors_R2.fastq.gz -CLIP_ATTR XT -CLIP_ACT 2 -NON_PF true
echo "############################################# Sam To Fastq End #############################################"

## 4. Bwa Alignment
echo "############################################ Bwa Alignment Start ############################################"  
bwa mem -t 80 -M $REF/ucsc.hg19.fasta unmapped_markilluminaadaptors_R1.fastq.gz unmapped_markilluminaadaptors_R2.fastq.gz -o Alignment.sam
echo "############################################# Bwa Alignment End #############################################"

## 5. Sorting & BAM conversion
echo "############################################ Sorting & BAM conversion Start ############################################"
java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar SortSam -I Alignment.sam -O Sorted_Alignment.bam -SO coordinate
echo "############################################# Sorting & BAM conversion End #############################################"

## 6. Alignment stats with Samtools
echo "############################################ Alignment stats with Samtools Start ############################################"  
samtools flagstat Sorted_Alignment.bam > Alignment_Stats.txt
echo "############################################# Alignment stats with Samtools End #############################################"

## 7. Merge BAM Start #Ref should have Seq dictionary (i.e Reference.dict file) and .fai file
echo "############################################ Merge BAM Start ############################################" 
java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar MergeBamAlignment -ALIGNED Sorted_Alignment.bam -UNMAPPED unmapped.bam -O Mapped.bam -R $REF/ucsc.hg19.fasta -MC true --CLIP_ADAPTERS false --CLIP_OVERLAPPING_READS true --INCLUDE_SECONDARY_ALIGNMENTS true --MAX_INSERTIONS_OR_DELETIONS -1 --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS
echo "############################################# Merge BAM End #############################################"

## 8. BamIndex
echo "############################################ BamIndex Start #############################################"
java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar BuildBamIndex -I Mapped.bam
echo "############################################# BamIndex End ##############################################"

## 9. Validate Sam
echo "############################################ Validate Sam Start #########################################"
java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar ValidateSamFile -I Mapped.bam --REFERENCE_SEQUENCE $REF/ucsc.hg19.fasta -M SUMMARY -O SamValidation_metrics.txt
echo "############################################ Validate Sam End #########################################"

## 10. MArkDuplicate
echo "############################################ MarkDuplicate Start ############################################"
java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar MarkDuplicates -I Mapped.bam -O Markdup.bam -M Markduplicate_metrics.txt --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500
echo "############################################# MarkDuplicate End #############################################"

## 11. BamIndex
echo "############################################ BamIndex Start ############################################" 
java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar BuildBamIndex -I Markdup.bam
echo "############################################# BamIndex End #############################################"

## 12. BASE RECALIBERATION ###########
## Base Quality Score Recalibration (BQSR) 1
echo "############################################ BQSR-1 Start ############################################" 
java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar BaseRecalibrator -R $REF/ucsc.hg19.fasta -I Markdup.bam -L $BED -ip 100 --known-sites $REF/dbsnp_138.hg19_sorted.vcf --known-sites $REF/Mills_and_1000G_gold_standard.indels.hg19.sites.sorted.vcf -O recal_data.table
echo "############################################# BQSR-1 End #############################################"

## Apply Recalibration 2
echo "############################################ Apply Recalibration-2 Start ############################################"  
java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar ApplyBQSR -bqsr recal_data.table -R $REF/ucsc.hg19.fasta -I Markdup.bam -O recal_reads.bam
echo "############################################# Apply Recalibration-2 End #############################################"

## 13. Index bam
echo "############################################ Bam Indexing Start ############################################" 
java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar BuildBamIndex -I recal_reads.bam
echo "############################################# Bam Indexing End #############################################"

## 14. Variant Calling
echo "############################################ Haplotypecaller Start ############################################" 
java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar HaplotypeCaller  -R $REF/ucsc.hg19.fasta -I recal_reads.bam -L $BED -ip 100 -stand-call-conf 30 -O ${line}_Raw_variants.vcf
echo "############################################# Haplotypecaller End #############################################"
# -ERC GVCF - Use this parameter for joint genotyping samples

## 15. Extract SNP and Idels
echo "############################################# SNP & Indel Extraction Start #############################################"
java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar SelectVariants -R $REF/ucsc.hg19.fasta -V ${line}_Raw_variants.vcf -select-type SNP -O Raw_snps.vcf

java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar SelectVariants -R $REF/ucsc.hg19.fasta -V ${line}_Raw_variants.vcf -select-type INDEL -O Raw_indels.vcf
echo "############################################## SNP & Indel Extraction End ##############################################"

## 16.1. SNP Filtering
echo "############################################# SNP Filtering Start #############################################"
java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar VariantFiltration -R $REF/ucsc.hg19.fasta -V Raw_snps.vcf --filter-expression "QD < 2.0" --filter-name "QualityByDepth2" --filter-expression "FS > 60.0" --filter-name "FisherStrand60" --filter-expression "MQ < 40.0" --filter-name "RMSMappingQuality40" --filter-expression "MQRankSum < -12.5" --filter-name "MappingQualityRankSumTest-12.5" --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSumTest-8" -O snps_filtered.vcf
echo "############################################## SNP Filtering End ##############################################"

## 16.2. Indel Filtering
echo "############################################# Indel Filtering Start #############################################"  
java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar VariantFiltration -R $REF/ucsc.hg19.fasta -V Raw_indels.vcf --filter-expression "QD < 2.0" --filter-name "QualityByDepth2" --filter-expression "FS > 200.0" --filter-name "FisherStrand200" --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O indels_filtered.vcf
echo "############################################## Indel Filtering End ##############################################"

## 17. Combine variants
echo "############################################# Combine Variants Start #############################################" 
java -Xmx400G -jar $GATK4/gatk-package-4.3.0.0-local.jar MergeVcfs -I snps_filtered.vcf -I indels_filtered.vcf -O ${line}_Final_Filtered_Variants.vcf
echo "############################################## Combine Variants End ##############################################"

echo "############################################# File Deletion Start #############################################" 
rm unmapped_markilluminaadaptors_R1.fastq.gz unmapped_markilluminaadaptors_R2.fastq.gz Alignment.sam Markdup.bai Markdup.bam recal_data.table Sorted_Alignment.bam unmapped.bam unmapped_markilluminaadapters.bam Raw_snps.vcf Raw_indels.vcf
echo "############################################## File Deletion End ##############################################"

echo "DATE"
date

done

## 18. Annovar
  # VCF to annovar input format conversion
#perl /data/applications/annovar/convert2annovar.pl --format vcf sample.vcf  --includeinfo --outfile KHGLBS38.avinput  --withzyg
  
  # Annotation
#perl ../../annovar_2020/table_annovar.pl 14111_S35.avinput ../../../humandb/ -buildver hg19 -out 14111_S35 -remove -protocol refGene,phastConsElements46way,genomicSuperDups,1000g2015aug_all,esp6500siv2_all,exac03_modified,gnomad211_exome_modified,gnomad211_genome_modified,cg69,Generic_In-house_Variants,Generic_All-Indian_Variants,snp138NonFlagged,avsnp138,avsnp150,clinvar_20200316,ljb26_sift,ljb26_pp2hdiv,ljb26_pp2hvar,ljb26_mt,caddindel,cadd13,mcap13,dbnsfp31a_interpro,intervar_20180118 -operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -genericdbfile hg19_Generic_In-house_Variants.txt -genericdbfile hg19_Generic_All-Indian_Variants.txt --nastring NA -csvout -otherinfo --nopolish

##### Comments ######
## 1. This script takes sample.fastq.gz (Both F & R reads), BWA indexed genome file, database files(Ex: dbsnp, 1000G etc) sorted and indexed according to the genome file.
## 2. This script when run, starts with fastq.gz files and gives out .vcf which is the final filtered variants file which must lateron be annotated using annovar.
## 3. The user can see labels like '##1 ##2 ##3 ##4 ##5 ##6 ##7' above which are crucial for the successful execution of this script.
## 4. ##1 - Give the path of the working directory.
      ##2 - Path of the reference directory which should contain bwa indexed genome files (along with .dict & .fai file) and snp database files
      #     (ex. dbsnp, 1000g etc) sorted according to the genome .dict file.
      ##3 - GATK4 path in your local computer or server.
      ##4 - The bed file path used for exome capture.
      ##5 - Make sure that your fastq files are named as 'Sample*_S*_R1_001.fastq.gz' (ex: 14591_S35_R1_001.fastq.gz).
      ##    In the fastq filename, 'R1/R2_001' is mandatory. (If your fastq files does not have R1/R2_001, 
      ##    please rename your fastq files accordingly instead of editing the script.
      ##6 - Path of the working directory where all the fastq files are copied.
      ##7 - Same path as in ##6 with '/$line/' extension must be given.
