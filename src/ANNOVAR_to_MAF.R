annotated_files<-list.files(path ="/Users/hhakimjavadi/Dropbox (UFL)/mixture_dna/data/case1/annotated", full.names = TRUE)

var.tbl<-data.frame()
for (annot in annotated_files[2:4]) {
  #annovar.tbl = annovarToMaf(annovar = annot, refBuild = 'hg19', tsbCol = 'Tumor_Sample_Barcode', table = 'ensGene')
  #annovar.tbl = annovarToMaf(annovar = annot, refBuild = 'hg19', table = 'ensGene', ens2hugo=TRUE)
  annovar.tbl = annovarToMaf(annovar = annot, refBuild = 'hg19') %>% filter(Otherinfo10=="PASS") %>%
    filter(Chromosome=="MT")
  annovar.tbl$AF<-as.numeric(annovar.tbl$AF)
  annovar.tbl<-filter(annovar.tbl ,AF<0.05 | AF==".") #this is a soft filter to avoid SNP infiltration and reduce the file size
  var.tbl<-rbind(var.tbl, annovar.tbl)
}
var.tbl<-data.frame()

annovar.tbl = annovarToMaf(annovar = annotated_files[4], refBuild = 'hg19_MT') %>% filter(Otherinfo10=="PASS") %>%
  filter(Chromosome=="MT")

annovar.tbl = annovarToMaf(annovar = annot, refBuild = 'hg19_MT', table = 'ensGene', ens2hugo=TRUE) %>% filter(Chromosome=="MT")
annovar.tbl$AF<-as.numeric(annovar.tbl$AF)
var.tbl<-rbind(var.tbl, annovar.tbl)
