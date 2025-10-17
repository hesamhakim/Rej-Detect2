library(SICtools)

hs37d5="/Users/hhakimjavadi/reference_data/hs37d5/hs37d5.fa"

bam_806_recipent <-"/Users/hhakimjavadi/Downloads/mixture_vcf/bam/CPM00023806-B-D_20220926-dragen.bam"
bam_805_donor <-"/Users/hhakimjavadi/Downloads/mixture_vcf/bam/CPM00023805-B-D_20220926-dragen.bam"
bam_804_cfDNA <-"/Users/hhakimjavadi/Downloads/mixture_vcf/bam/CPM00023804-PL-D_20220926-dragen.bam"

min_BaseQuality = 13
min_MapQuality = 10
pValue_CutOff = 1
base_DistCutOff=0.1

donor_vs_cfDNA<-snpDiff(bam1=bam_805_donor, bam2=bam_804_cfDNA, 
                       refFsa=hs37d5, regChr="MT", 
                       regStart=1, regEnd=16500, 
                       minBaseQuality = min_BaseQuality, minMapQuality = min_MapQuality, 
                       pValueCutOff = pValue_CutOff,
                       baseDistCutOff=base_DistCutOff,
                       nCores = 4)
donor_vs_cfDNA$comparison="donor_vs_cfDNA"
write.csv(donor_vs_cfDNA,"/Users/hhakimjavadi/Dropbox (UFL)/mixture_dna/tables/donor_vs_cfDNA.csv")


recipient_vs_cfDNA<-snpDiff(bam1=bam_806_recipent, bam2=bam_804_cfDNA, 
                       refFsa=hs37d5, regChr="MT", 
                       regStart=1, regEnd=16500, 
                       minBaseQuality = min_BaseQuality, minMapQuality = min_MapQuality, 
                       pValueCutOff = pValue_CutOff,
                       baseDistCutOff=base_DistCutOff,
                       nCores = 4)
recipient_vs_cfDNA$comparison="recipient_vs_cfDNA"
write.csv(recipient_vs_cfDNA,"/Users/hhakimjavadi/Dropbox (UFL)/mixture_dna/tables/recipient_vs_cfDNA.csv")


donor_vs_recipient<-snpDiff(bam1=bam_805_donor, bam2=bam_806_recipent, 
                            refFsa=hs37d5, regChr="MT", 
                            regStart=1, regEnd=16500, 
                            minBaseQuality = min_BaseQuality, minMapQuality = min_MapQuality, 
                            pValueCutOff = pValue_CutOff,
                            baseDistCutOff=base_DistCutOff,
                            nCores = 4)
donor_vs_recipient$comparison="donor_vs_recipient"
write.csv(donor_vs_recipient,"/Users/hhakimjavadi/Dropbox (UFL)/mixture_dna/tables/donor_vs_recipient.csv")

all_comparisons=rbind(recipient_vs_cfDNA, donor_vs_cfDNA, donor_vs_recipient)

all_comparisons= dplyr::mutate(all_comparisons, 
                               bam1_total_reads=(A1+C1+G1+T1),
                               bam2_total_reads=(A2+C2+G2+T2),
                               "bam1:A AF"=A1/bam1_total_reads,
                               "bam1:C AF"=C1/bam1_total_reads,
                               "bam1:G AF"=G1/bam1_total_reads,
                               "bam1:T AF"=T1/bam1_total_reads,
                               "bam2:A AF"=A2/bam2_total_reads,
                               "bam2:C AF"=C2/bam2_total_reads,
                               "bam2:G AF"=G2/bam2_total_reads,
                               "bam2:T AF"=T2/bam2_total_reads
                               )
all_comparisons[all_comparisons$pos==11719,] %>% view()
write.csv(all_comparisons,"/Users/hhakimjavadi/Dropbox (UFL)/mixture_dna/tables/all_comparisons.csv")
