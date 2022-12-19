files<-dir('./orthologs_for_paramecium')

files<-files[-c(1, length(files))]

ortho<-integer()
for(i in seq_along(files)){
  dat<-read.table(paste0('./orthologs_for_paramecium/', files[i]), sep=',', header=TRUE, quote = "")
  orthologs<-length( unique(dat$gene[dat$ortholog=='LDO']) )
  
  ortho<-c(ortho, orthologs)
}


geneNames<-read.table( dir('./orthologs_for_paramecium', full.names = TRUE)[7], sep=',', header=TRUE )


ortho.counts<-data.frame(files, ortho, percent = ortho/nrow(geneNames)*100)

round(mean(ortho.counts$percent), 2)

round(sd(ortho.counts$percent), 2)

overlap<-matrix(rep(NA, nrow(ortho.counts)^2), ncol=nrow(ortho.counts))

sets<-list()

for(i in seq_along(files)){
  ref<-read.table(paste0('./orthologs_for_paramecium/', files[i]), sep=',', header=TRUE, quote = "")
  ref<-unique(ref$gene[ref$ortholog=='LDO'])
  sets[[i]]<-ref
  
  for(j in seq_along(files)[-i]){
    comp<-read.table(paste0('./orthologs_for_paramecium/', files[j]), sep=',', header=TRUE, quote = "")
    comp<-unique(comp$gene[comp$ortholog=='LDO'])
    
    num<-sum(!is.na(match(ref, comp)) )
    
    overlap[i,j]<-num
    
  }
  
}

overlap<-data.frame(overlap)
names(overlap)<-c("fly" , "amoeba" , "human", "mouse", "yeast")
row.names(overlap)<-c("fly" , "amoeba" , "human", "mouse", "yeast")


if(!require(VennDiagram)){
  install.packages('VennDiagram')
}

# Load library
library(VennDiagram)

# Chart
venn.diagram(
  x = sets,
  category.names = names(overlap),
  filename = 'test_venn_diagramm.png',
  output=TRUE,
)





################
# sets2<-sets
# 
# for(k in seq_along(sets2)){
#   
#   sets2[[k]]<-unlist(lapply(sets2[[k]], function(x){
#   substr(strsplit(x, "=")[[1]][2], 1, 17) 
# }))  
# 
# }
# 
# sets2[[6]]<-geneNames$GSPATT 
#   
#   
# venn.diagram(
#   x = sets2,
#   category.names = c(names(overlap), 'Paramecium'),
#   filename = 'test2_venn_diagramm.png',
#   output=TRUE
# )

