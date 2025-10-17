
# this script generates simulated base compositions based on mitochondrial reference genome
#simulated samples have germline variants and noise variants

setwd("/Users/hhakimjavadi/Dropbox (UFL)/mixture_dna")
##to generate mt reference genome base composition
mito_ref_orig=read.table("~/Dropbox (UFL)/T_ALL_mtDNA/data/annotations/mt_annotations/MAESTER_mt_variants/mito_reference.txt", stringsAsFactors = F) # https://www.genome.jp/dbget-bin/www_bget?-f+refseq+NC_012920
mito_ref=NULL
for(i in 1:nrow(mito_ref_orig)){
  mito_ref=paste0(mito_ref, mito_ref_orig[i,])
}
mito_ref=toupper(mito_ref)
mito_ref=unlist(strsplit(mito_ref, split=""))

mt_length=length(mito_ref)
############
#mt_length=20
#mito_ref=mito_ref[1:20]
##########

mt_ref_tbl<-matrix(nrow = mt_length, ncol = 2) %>% as.data.frame()
mt_ref_tbl[,"POS"]=1:mt_length
mt_ref_tbl[,"REF"]=mito_ref[1:mt_length]
mt_ref_tbl[,"DP"]=runif(mt_length, 12, 50) %>% as.integer()
mt_ref_tbl[, c("A", "C", "T", "G")]=0
mt_ref_tbl[mt_ref_tbl$REF=="A", "A"]=mt_ref_tbl[mt_ref_tbl$REF=="A", "DP"]
mt_ref_tbl[mt_ref_tbl$REF=="C", "C"]=mt_ref_tbl[mt_ref_tbl$REF=="C", "DP"]
mt_ref_tbl[mt_ref_tbl$REF=="G", "G"]=mt_ref_tbl[mt_ref_tbl$REF=="G", "DP"]
mt_ref_tbl[mt_ref_tbl$REF=="T", "T"]=mt_ref_tbl[mt_ref_tbl$REF=="T", "DP"]


####we modify the reference base compositon by variants and noises

####function to simulate genome with germline variant

generate_germline<-function(var_count, var_min_dp, var_max_dp,
                            noise_count, noise_max_dp, 
                            genome_name, index) {
  germline_variant_pos=runif(var_count, 1, mt_length) %>% as.integer()
  noise_pos=runif(noise_count, 1, mt_length) %>% as.integer()
  var_pos_ref_allele<-mito_ref[germline_variant_pos]
  noise_pos_ref_allele<-mito_ref[noise_pos]
  
  ##to create two alternative genome: 1) with germline variants 2) with noise variants
  alt_vars=c()
  noise_vars=c()
  for (A in 1:var_count) {
    alt_var=setdiff(c("A", "C", "T", "G"), mito_ref[germline_variant_pos[A]])[as.integer(runif(1,1,4))]
    alt_vars=c(alt_vars, alt_var)
  }
  
  for (A in 1:noise_count) {
    noise_var=setdiff(c("A", "C", "T", "G"), mito_ref[noise_pos[A]])[as.integer(runif(1,1,4))]
    noise_vars<-c(noise_vars, noise_var)
  }
  
  ##create alternative base compositions
  mito_genome<- mito_ref
  noise_genome<-mito_ref
  mito_genome[germline_variant_pos]=alt_vars
  noise_genome[noise_pos]=noise_vars
  
  mt_tbl<-matrix(nrow = mt_length, ncol = 2) %>% as.data.frame()
  mt_tbl[,"POS"]=1:mt_length
  mt_tbl[, c("A", "C", "T", "G")]=0
  
  mt_var_tbl<-mt_tbl
  mt_var_tbl[,"REF"]=mito_ref[1:mt_length]
  mt_var_tbl[,"Allele"]=mito_genome
  mt_var_tbl[,"DP"]=runif(mt_length, var_min_dp, var_max_dp) %>% as.integer()
  
  
  mt_var_tbl[mt_var_tbl$Allele=="A", "A"]=mt_var_tbl[mt_var_tbl$Allele=="A", "DP"]
  mt_var_tbl[mt_var_tbl$Allele=="C", "C"]=mt_var_tbl[mt_var_tbl$Allele=="C", "DP"]
  mt_var_tbl[mt_var_tbl$Allele=="G", "G"]=mt_var_tbl[mt_var_tbl$Allele=="G", "DP"]
  mt_var_tbl[mt_var_tbl$Allele=="T", "T"]=mt_var_tbl[mt_var_tbl$Allele=="T", "DP"]
  mt_var_tbl<-mt_var_tbl[,-c(1:2)]
  
  ##here we introduce noises
  mt_noise_tbl<-mt_tbl
  mt_noise_tbl[,"REF"]=mito_ref[1:mt_length]
  mt_noise_tbl[, "Noise_Allele"]=noise_genome
  mt_noise_tbl[mt_noise_tbl$Noise_Allele=="A", "A"]=as.integer(runif(1, 1, noise_max_dp))
  mt_noise_tbl[mt_noise_tbl$Noise_Allele=="C", "C"]=as.integer(runif(1, 1, noise_max_dp))
  mt_noise_tbl[mt_noise_tbl$Noise_Allele=="G", "G"]=as.integer(runif(1, 1, noise_max_dp))
  mt_noise_tbl[mt_noise_tbl$Noise_Allele=="T", "T"]=as.integer(runif(1, 1, noise_max_dp))
  mt_noise_tbl<-mt_noise_tbl[,-c(1:2)]
  
  ###create sample that include both germline and noise variants
  mt_samp<-mt_var_tbl
  mt_samp[,c(2:5)]=mt_var_tbl[,c(2:5)]+mt_noise_tbl[,c(2:5)]
  mt_samp[, "Noise_Allele"]=noise_genome
  
  mt_samp[noise_pos, "noise_var_id"]=NA
  mt_samp[noise_pos,] <-mt_samp[noise_pos,] %>% dplyr::mutate(noise_var_id=paste(REF,POS,Noise_Allele, sep = ""))
  mt_samp[germline_variant_pos, "var_id"]=NA
  mt_samp[germline_variant_pos,] <-mt_samp[germline_variant_pos,] %>% dplyr::mutate(var_id=paste(REF,POS,Allele, sep = ""))
  

  #reorder table
  mt_samp_ordered<-mt_samp
  mt_samp_ordered$CHR="MT"
  mt_samp_ordered<-mt_samp_ordered[,c("CHR","POS","REF","Allele","Noise_Allele","A", "C","T","G",
                                      "DP","var_id", "noise_var_id")]

  write.csv(mt_samp_ordered, paste("/Users/hhakimjavadi/Dropbox (UFL)/mixture_dna/simulation/simulated_base_compositions/sim_", 
                           genome_name,"_", index,".csv", sep = ""),row.names = FALSE, na="")
}


var_count=15
var_min_dp=30
var_max_dp=50
noise_count=50
noise_max_dp=10
genome_name="set1"


for (I in 1:50) {
  generate_germline(var_count=var_count, var_min_dp=var_min_dp, var_max_dp = var_max_dp,
                    noise_count = noise_count,noise_max_dp =noise_max_dp ,
                    genome_name=genome_name, index = I)
}

