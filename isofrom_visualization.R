

#R中尝试绘制一个基因，不同转录本
#获取一个基于所有转录本
library(ggplot2)
library(ggpubr)
# 
# tmp=mm_gtf %>% filter(genesymbol=='Egfr') %>%
#   mutate(exon_L=abs(start-end)) %>%
#   group_by(transcriptid) %>%
#   summarise(N=n(),
#             N_start=min(start),
#             transcipt_L=sum(exon_L)) %>%
#   mutate(y1=row_number(),
#          y2=y1+0.5,
#          y3=y1+0.25,
#          transcriptid_label=paste0(transcriptid,":",transcipt_L,"bp")) %>% ungroup()
# p1=mm_gtf %>% filter(genesymbol=='Egfr') %>%
#   left_join(.,
#             tmp,
#             by=c('transcriptid'='transcriptid')) %>%
#   ggplot()+
#   geom_hline(yintercept = c(tmp$y3),lty=1,lwd=1.2)+
#   geom_rect(aes(xmin=log1p(start),xmax=log1p(end),ymin=y1,ymax=y2),fill='lightblue')+
#   geom_text(aes(x=(start+end)/2,y=y2,label=exon))+
#   #ggtitle(label = 'Trp53')+
#   ylim(c(0,7))
#   #theme_void()
# # p2=tmp %>% ggplot(aes(y=y3,x=1))+geom_text(aes(label=transcriptid_label))+ylim(c(0,7))+theme_void()
# # ggarrange(p2,NULL,p1,widths = c(2,-0.2,3),nrow = 1)

# mm_gene_list=unique(mm_gtf$genesymbol)
# hm_gene_list=unique(hm_gtf$genesymbol)
# saveRDS(mm_gene_list,'datasets/mm_gene_list.rds')
# saveRDS(hm_gene_list,'datasets/hm_gene_list.rds')
# saveRDS(mm_gtf,'datasets/mm_gtf.rds')
# saveRDS(hm_gtf,'datasets/hm_gtf.rds')

exons_in_isoforms<-function(gene,species='Human'){
  mm_gene_list=readRDS('./datasets/mm_gene_list.rds')
  mm_gtf=readRDS('./datasets/mm_gtf.rds')
  
  hm_gene_list=readRDS('./datasets/hm_gene_list.rds')
  hm_gtf=readRDS('./datasets/hm_gtf.rds')
  
  require(ggplot2)
  require(ggpubr)
  require(dplyr)
  #gene_list: mm_gene_list, hm_gene_list
  #gtf: mm_gtf,hm_gtf
  if(species=='Human'){
    gene_list=hm_gene_list
    gtf=hm_gtf
  }else{
    gene_list=mm_gene_list
    gtf=mm_gtf
  }
  
  if(gene %in% gene_list){
    tmp=gtf %>% filter(genesymbol== gene) %>% 
      mutate(exon_L=abs(start-end)) %>% 
      group_by(transcriptid) %>% 
      summarise(N=n(),
                N_start=min(start),
                transcipt_L=sum(exon_L)) %>% 
      mutate(y1=row_number(),
             y2=y1+0.5,
             y3=y1+0.25,
             transcriptid_label=paste0(transcriptid,": ",transcipt_L,"bp")) %>% 
      ungroup()
    
    N_isoforms=nrow(tmp)+1
    p1=gtf %>% filter(genesymbol== gene) %>% 
      left_join(.,
                tmp,
                by=c('transcriptid'='transcriptid')) %>% 
      ggplot()+
      geom_hline(yintercept = c(tmp$y3),lty=1,lwd=1.2)+
      geom_rect(aes(xmin=start,xmax=end,ymin=y1,ymax=y2),fill='lightblue')+
      geom_text(aes(x=(start+end)/2,y=y2,label=exon))+
      ylim(c(0,N_isoforms))+
      theme_void()
    
    p2=tmp %>% ggplot(aes(y=y3,x=1))+
      geom_text(aes(label=transcriptid_label))+
      ylim(c(0,N_isoforms))+
      theme_void()
    
    ps=ggarrange(p2,NULL,p1,widths = c(2,-0.2,3),nrow = 1)
  }else{
    ps=NULL
  }
  
  return(ps)
}
exons_in_isoforms('TET2')
exons_in_isoforms('Kras',species = 'Mus')
#如果我要标记一个位置呢？

