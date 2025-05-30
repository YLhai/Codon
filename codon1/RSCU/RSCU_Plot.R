# remove(list = ls()) #清除 Global Environment
setwd("D:/pythonProject/Codon/RSCU")

library(tidyverse)
seqinr::read.fasta("Post.cds") %>% 
  unlist() %>% 
  seqinr::uco(index="rscu") %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  magrittr::set_colnames(c("codon","rscu"))

codon.table=data.frame(read.csv('codon.table.csv'))
if(0){
data.frame(
  acids=c("Isoleucine","Leucine","Valine","Phenylalanine",
          "Methionine","Cysteine","Alanine","Glycine",
          "Proline","Threonine","Serine",
          "Tyrosine","Tryptophan","Glutamine",
          "Asparagine","Histidine","Glutamic acid",
          "Aspartic acid","Lysine","Arginine","Stop codons"),
  slc=c("I","L","V","F","M","C","A","G","P","T",
        "S","Y","W","Q","N","H","E","D","K","R","Stop"),
  codon=c("ATT, ATC, ATA","CTT, CTC, CTA, CTG, TTA, TTG",
          "GTT, GTC, GTA, GTG","TTT, TTC",
          "ATG","TGT, TGC",
          "GCT, GCC, GCA, GCG","GGT, GGC, GGA, GGG",
          "CCT, CCC, CCA, CCG","ACT, ACC, ACA, ACG",
          "TCT, TCC, TCA, TCG, AGT, AGC",
          "TAT, TAC","TGG","CAA, CAG","AAT, AAC","CAT, CAC","GAA, GAG","GAT, GAC",
          "AAA, AAG","CGT, CGC, CGA, CGG, AGA, AGG",
          "TAA, TAG, TGA")
  
) %>% 
  separate(codon,paste0("col",1:6),sep=", ") %>% 
  pivot_longer(!c(acids,slc)) %>% 
  na.omit() %>% 
  select(-name) %>% 
  magrittr::set_colnames(c("acids","abbr","codon"))-> codon.table
}

rscu.dat<-read.csv('Pcit.RSCU.csv')

if(0){
seqinr::read.fasta("Post.cds") %>% 
  unlist() %>% 
  seqinr::uco(index="rscu") %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  magrittr::set_colnames(c("codon","rscu")) %>% 
  mutate(codon=str_to_upper(codon)) %>% 
  left_join(codon.table,by=c("codon"="codon")) %>% 
  mutate(acids=str_sub(acids,1,3)) %>% 
  arrange(acids,rscu)%>% 
  group_by(acids) %>% 
  mutate(id=row_number()) %>% 
  dplyr::select(codon,acids,rscu,id) %>% 
  mutate(codon=factor(codon,level=codon),
         id=-id)  -> rscu.dat
}

rscu.dat

colnames(rscu.dat)<-paste0("V",1:4)



ggplot(rscu.dat,aes(fill=as.character(V4),x=V2,y=V3))+
  geom_bar(position = "stack",stat="identity")+
  theme_bw()+scale_y_continuous(expand=c(0,0),
                                limits = c(0,6.2))+
  theme(legend.position = "none")+labs(y="RSCU",x="")+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("#999999","#56b4e8","#d55e00",
                               "#009f73","#cc79a7","#0072b1"))+
  theme(plot.margin = unit(c(0.1,0.1,0,0.1),'mm'))-> rscu.p1

ggplot(rscu.dat,aes(x=V2,y=V4))+
  geom_label(aes(label=V1,fill=as.character(V4)),
             size=3)+
  labs(x="",y="")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values = c("#999999","#56b4e8","#d55e00",
                               "#009f73","#cc79a7","#0072b1"))+
  theme(plot.margin = unit(c(0,0.1,0.1,0.1),'mm'))+
  coord_cartesian(clip = "off") -> rscu.p2

library(patchwork)

cairo_pdf("RSCU_barplot.pdf",
          width = 9.4,height = 4)
rscu.p1+rscu.p2 + plot_layout(ncol=1,heights = c(2,1))
dev.off()