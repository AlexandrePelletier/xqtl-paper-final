setwd('/adpelle1/xqtl-paper-final/main_text/5_AD_xQTL_genes_cis_trans/staging/gene_priorization_table/')

#PARAMETERS/ENV####
out<-'.'
source('gene_prio_utils.R')
#install.packages('ggtext')
library(ggtext)
library(patchwork)
library(RColorBrewer)

context_groups=fread('context_groups_summfigure.tsv')

contexts_order=unique(context_groups$context_group)

#INPUT####
gwas_ncases_controls='/adpelle1/xqtl-resources/data/GWAS/gwas_n_cases_control.tsv'
long_table=fp(out,'xQTL_all_methods_overlap_with_AD_loci_unified_cs95orColocs_Pval1e5_noTrans.csv.gz') #generated in complete_ADlocus_level_summary.R
gwas_table=fp(out,'AD_loci_unified_cs95orColocs.csv.gz') #generated in complete_ADlocus_level_summary.R



#FUNCTIONS####
SummPanel<-function(res_adxqtloc,res_adxloc_gwas,ndatasets_range=NULL,palette=brewer.pal(8, "Pastel1")){
  #filter the gwas and annot with the gene too
  res_adxloc_gwasf<-merge(res_adxloc_gwas,
                          unique(res_adxqtloc[,.(locus_index,locus_gene_3,locus_gene)]),by='locus_index',allow.cartesian = TRUE)
  
  res_adxqtloc[,n_study_group.locus:=as.integer(n_study_group.locus)]
  p<-ggplot(res_adxqtloc)+
    geom_point(aes(y=locus_gene_3,x=context_group,
                   size=n_study_group.locus,
                   col=confidence_cat_group))+
    facet_grid(chr~'',scales = 'free',space = 'free')+theme_minimal()+
    scale_x_discrete(guide = guide_axis(angle = 90))+
    theme(strip.text.x = element_text(angle = 90))+
    labs(size='# datasets',col='Confidence level')+
    scale_color_manual(values = c('brown1','deepskyblue4','darkseagreen3'))+ 
    theme(strip.text.x = element_text(angle = 90),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())+
    scale_y_discrete(labels=function(x) str_remove(x, "_[0-9]+"))
  
  if(!is.null(ndatasets_range)){
    p<-p+scale_radius(limits = ndatasets_range,range = c(1.5,5),breaks =ceiling(as.vector(quantile(1:ndatasets_range[2],c(0.05,0.25,0.5,0.75,0.95)))) )
  }
  
  g<-ggplot(res_adxloc_gwasf)+
    geom_point(aes(y=locus_gene_3,x=gwas_short2,
                   col=gwas_sig),size=1,shape=15)+
    facet_grid(chr~'',scales = 'free',space = 'free')+
    # scale_size(range = c(0.5,2))+
    theme_minimal()+
    scale_x_discrete(guide = guide_axis(angle = 90))+
    scale_color_manual(values = c('grey','bisque3','orange3','brown4'))+
    theme(axis.text.x = element_markdown(size = 7),
          axis.text.y = element_blank(),
          axis.title.y =element_blank(),
          strip.text = element_blank(),   # remove facet text
          strip.background = element_blank() 
    )
  #add the locus number
  colors <- data.table(locus_index=unique(res_adxqtloc$locus_index),
                       color=palette)
  
  l<-ggplot(unique(res_adxqtloc,
                   by=c('locus_gene_3')),
            aes(x=0,y=locus_gene_3))+
    geom_tile(aes(fill=as.factor(locus_index)))+scale_fill_manual(values = setNames(colors$color,as.character(colors$locus_index)))+
    geom_text(aes(label=locus_index),size=2)+
    facet_grid(chr~'',scales = 'free',space = 'free')+theme_minimal()+
    theme(axis.text.y = element_markdown(size=6),
          axis.title.y = element_blank(), 
          strip.text = element_blank(),   # remove facet text
          strip.background = element_blank(),
          axis.text.x=element_blank(),axis.title.x = element_blank(),
          legend.position = "none")+
    scale_y_discrete(labels=function(x) str_remove(x, "[A-Za-z0-9]+_[2-9]|_[0-9]+"))
  
  
  return(list(l,g,p))
}

#MAIN####
res_adx<-fread(long_table)

#add the confidence level

res_adx[,top_confidence:=str_extract(xQTL_effects,'CL[0-9]')]


#xQTL evidence summary dot plot at locus level
#for each QTL evidence type put in 3 cat : absent, present, high, 
#for fp: high if cs95, coloc : high if npc > 0.95, TWAS: MR
table(res_adx$Method)
table(res_adx$context_short)
setdiff(res_adx$context_short,contexts_order)
table(res_adx$context)


res_adx[,locus_gene:=paste(locus_index,ifelse(gene_name=='','?',gene_name),sep='_')]


#annot
res_adx[,genomewide_sig_gene:=any(min_pval<5e-8,na.rm = T),by='gene_name']

res_adx[,n_gwas_gene:=strsplit(gwas_sources,'\\|')|>unlist()|>unique()|>length(),by=.(gene_name)]
res_adx[,n_gwas_locus:=strsplit(gwas_sources,'\\|')|>unlist()|>unique()|>length(),by=.(locus_index)]

res_adx[is.na(APOE_region),APOE_region:=FALSE]

res_adx[,APOE_region.gene:=any(APOE_region),by='gene_name']


res_adx[,chr:=seqid(ADlocusID[1],only_num = T),by=.(ADlocusID)]

res_adx[,MAPT_region:=chr==17&pos>45307498 &pos> 46836265]
res_adx[is.na(MAPT_region),MAPT_region:=FALSE]

res_adx[,MAPT_region.locus:=any(MAPT_region),by='locus_index']

#get at locus level and cis information only
res_adxloc<-unique(res_adx[Method!='trans_finemapping'][!context_short%in%c('B','J','K','W')][order(locus_gene,confidence_lvl,cV2F_rank)],by=c('locus_gene','context'))

#top gene per locus
res_adxloc[gene_name!='',n.study.gene.locus:=length(unique(context)),by=.(gene_name,locus_index)]

res_adxloc[gene_name!='',top_gene:=gene_name==gene_name[order(confidence_lvl,-n.study.gene.locus,cV2F_rank)][1],by=.(locus_index)]
res_adxloc[gene_name!='',top_genes:=gene_name%in%gene_name[confidence_lvl%in%c('CL1','CL2','CL3')|confidence_lvl%in%sort(confidence_lvl)[1]],by=.(locus_index)]


# #TopGenes per locus####
# res_adxlocge<-res_adxloc[context_short!=''][order(locus_index,cV2F_rank)][(top_genes)]
# unique(res_adxlocge$locus_index)|>length() #177
# unique(res_adxlocge$gene_name)|>length()#283>304/576
# 
# #
# # group some contexts
# 
# res_adxlocge[,context_group:=ifelse(context_short%in%c('bMono eQTL','bMac eQTL','bMic eQTL'),'Immune eQTL',
#                                      ifelse(context_short%in%c( 'bulk p-sQTL','bulk u-sQTL','bulk a-sQTL'),'bulk sQTL',
#                                             ifelse(context_short%in%c('bulk pQTL','bulk gpQTL'),'bulk (g)pQTL',
#                                                    ifelse(context_short%in%c('bulk haQTL','bulk mQTL'),'bulk epiQTL',as.character(context_short)))))]
# 
# 
# unique(res_adxlocge$context_group)
# 
# unique(res_adxlocge[order(confidence_lvl)],by='gene_name')$confidence_lvl|>table()
# 
# 
# # summarize using confidence score  and nstudy
# res_adxlocge[,n_study_group:=length(unique(context)),by=.(gene_name,context_group)]
# res_adxlocge[,confidence_lvl_group:=sort(confidence_lvl)[1],by=.(gene_name,context_group)]
# table(res_adxlocge$confidence_lvl_group)
# res_adxlocge[,n_loci_gene:=unique(locus_index)|>length(),by=.(gene_name)]
# res_adxlocge[,chr:=seqid(variant_ID[!is.na(variant_ID)][1],only_num = T),by=.(locus_index)]
# res_adxlocge[is.na(chr),chr:=`#chr`]
# table(res_adxlocge$confidence_lvl_group)
# 
# 
# res_adxlocge[,n_study_group.locus:=length(unique(context)),by=.(gene_name,locus_index,context_group)]
# res_adxlocge[,confidence_lvl_group.locus:=sort(confidence_lvl)[1],by=.(gene_name,locus_index,context_group)]
# 
# 
# 
# fwrite(res_adxlocge,fp(out,'res_summ_top_gene_by_locus.csv.gz'))
# 
# res_adxlocge<-fread(fp(out,'res_summ_top_gene_by_locus.csv.gz'))
# 
# #genome wide signifcant htis
# 
# contexts_order<-c('Exc eQTL','Exc sQTL',
#                   'Inh eQTL','Inh sQTL',
#                   'Oli eQTL','Oli sQTL',
#                   'OPC eQTL','OPC sQTL',
#                   'Ast eQTL','Ast caQTL','Ast sQTL',
#                   'Mic eQTL','Mic caQTL',
#                   'Immune eQTL',
#                   'bulk eQTL',
#                   'bulk sQTL',
#                   'bulk (g)pQTL',
#                   'bulk epiQTL')
# setdiff(res_adxlocge$context_group,contexts_order)
# res_adxlocge[,context_group:=factor(context_group,levels = contexts_order)]
# 
# 
# res_adxlocgef<-res_adxlocge[(genomewide_sig_gene)]
# 
# unique(res_adxlocgef$locus_index)|>length() #98
# unique(res_adxlocgef$gene_name)|>length()#79> 116 >134
# 
# unique(res_adxlocgef$locus_gene)|>length()#98> 138 >156
# unique(res_adxlocgef[!(APOE_region.gene)]$locus_gene)|>length()#86>122>140
# 
# 
# # extract each locus
# res_adxlocge_cont_top<-unique(res_adxlocgef,by=c('gene_name','context_group','locus_index'))[gene_name%in%gene_name[confidence_lvl_group.locus%in%c('CL1','CL2','CL3','CL4','CL5')]]
# 
# 
# #sep causal vs correlated and remove all C6
# #C6 is dropped
# res_adxlocge_cont_topf<-res_adxlocge_cont_top[confidence_lvl_group.locus!='CL6']
# #sep in 3categories
# res_adxlocge_cont_topf[,confidence_cat_group:=ifelse(confidence_lvl_group.locus%in%c('CL1',"CL2"),'Putative causal (CL1, CL2)',ifelse(confidence_lvl_group.locus=='CL3','Putative causal (CL3)','Associated (CL4, CL5)'))]
# res_adxlocge_cont_topf[,confidence_cat_group:=factor(confidence_cat_group,levels = c('Putative causal (CL1, CL2)','Putative causal (CL3)','Associated (CL4, CL5)'))]
# 
# 
# #add GWAS signals
# res_ad<-fread(fp(out,'AD_loci_unified_cs95orColocs.csv.gz'))
# res_ad<-merge(res_ad,unique(res_adx[,.(locus_index,ADlocusID)]))
# 
# res_ad[,gwas_pvalue:=getPval(gwas_zscore)]
# 
# res_adxloc_gwas<-unique(res_ad[order(gwas_pvalue)][gwas_source!=''],by=c('locus_index','gwas_source'))
# res_adxloc_gwas[,chr:=str_extract(chr,'[0-9]+')|>as.numeric()]
# 
# res_adxloc_gwas[,gwas_pvalue10:=ifelse(gwas_pvalue>1e-10,gwas_pvalue,1e-10)]
# res_adxloc_gwas[,gwas_sig:=ifelse(gwas_pvalue<5e-8,'p<5e-8',ifelse(gwas_pvalue<1e-6,'p<1e-6',ifelse(gwas_pvalue<1e-5,'p<1e-5',
#                                                                                                       'p>1e-5 but coloc')))]
# res_adxloc_gwas[is.na(gwas_sig),gwas_sig:='p>1e-5 but coloc']
# res_adxloc_gwas[,gwas_sig:=factor(gwas_sig,levels=c('p>1e-5 but coloc','p<1e-5','p<1e-6','p<5e-8'))]
# 
# 
# res_adxloc_gwas[,gwas_short:=str_extract(gwas_source,'Bellenguez|Jansen|Kunkle|Wightman')]
# table(res_adxloc_gwas$gwas_short)
# res_adxloc_gwas[str_detect(gwas_source,'EADB|EADI|UKB|23andMe'),gwas_short:=paste(gwas_short,str_extract(gwas_source,'EADB|EADI|UKB|23andMe'),sep='_')]
# res_adxloc_gwas[gwas_short=='Wightman_23andMe',gwas_short:='Wightman_no23andMe']
# res_adxloc_gwas[gwas_short=='Wightman_UKB',gwas_short:='Wightman_no23AndUKB']
# 
# #order by n AD cases (non AD by proxi first) 
# res_adxloc_gwas[,ad_by_proxi:=!gwas_source%in%c('AD_Bellenguez_EADB_2022','AD_Bellenguez_EADI_2022',
#                                           'AD_Wightman_ExcludingUKBand23andME_2021',
#                                           'AD_Kunkle_Stage1_2019'),by='locus_index']
# 
# gwmt<-fread(gwas_ncases_controls)
# res_adxloc_gwas<-merge(res_adxloc_gwas,gwmt[,.(gwas_source=study_id,n_case,n_control)])
# 
# 
# 
# #bold if not from  AD by proxy 
# res_adxloc_gwas[,gwas_short2:=ifelse(!ad_by_proxi,paste0("<b>",gwas_short,'</b>'),
#                                             gwas_short)]
# 
# res_adxloc_gwas[,gwas_short2:=factor(gwas_short2,levels = unique(gwas_short2[order(ad_by_proxi,-n_case)]))]
# levels(res_adxloc_gwas$gwas_short2)
# 
# 
# 
# #annot genes
# res_adxlocge_cont_topf[,locus_gene_2:=paste(gene_name,match(locus_gene,
#                                                             unique(locus_gene[order(chr,tss,min_pval)])),sep='_'),
#                        by=.(gene_name)]
# 
# res_adxlocge_cont_topf[,locus_gene_2:=paste0("<i>",locus_gene_2[1],'</i>'),
#                        by='locus_gene']
# #add asterisk if others genes in a CL4+ locus
# 
# res_adxlocge_cont_topf[,n_gene_cl:=length(unique(gene_name)),by=c('locus_index','top_confidence')]
# res_adxlocge_cont_topf[,locus_gene_2:=paste0(locus_gene_2,ifelse(n_gene_cl>1&!top_confidence%in%c("CL1",'CL2','CL3'),
#                                                                  '*','')),
#                        by='locus_gene']
# #add if from  AD by proxy only
# res_adxlocge_cont_topf[,locus_gene_2:=ifelse(!all(only_by_proxi,na.rm = T),paste0("<b>",locus_gene_2[1],'</b>'),
#                                              locus_gene_2[1]),
#                        by='locus_gene']
# 
# #order per gene tss and locus position
# res_adxlocge_cont_topf[,locus_gene_3:=factor(locus_gene_2,levels = unique(locus_gene_2[order(-chr,-tss,-as.numeric(str_extract(locus_gene_2,'_[0-9+]')|>str_remove('_')))]))]
# 
# #filter the gwas and annot with the gene too
# res_adxloc_gwasf<-merge(res_adxloc_gwas,unique(res_adxlocge_cont_topf[!(APOE_region.gene)&(top_gene|top_confidence%in%c("CL1",'CL2','CL3'))][,.(locus_index,locus_gene_3,locus_gene)]),by='locus_index')
# 
# 
# 
# #Plot all genes having causal evidence per loci for the CL1-3, and only top1 gene (based on number of datasets) for the CL4+
# 
# #chr1-9
# p1<-ggplot(res_adxlocge_cont_topf[chr%in%1:9][!(APOE_region.gene)&(top_gene|top_confidence%in%c("CL1",'CL2','CL3'))])+
#   geom_point(aes(y=locus_gene_3,x=context_group,
#                  size=n_study_group.locus,
#                  col=confidence_cat_group))+
#   facet_grid(chr~'',scales = 'free',space = 'free')+
#   scale_size(range = c(1.5,5),breaks = c(2,7,12))+theme_minimal()+
#   scale_x_discrete(guide = guide_axis(angle = 90))+
#   theme(strip.text.x = element_text(angle = 90))+
#   labs(size='# datasets',col='Confidence level')+
#   scale_color_manual(values = c('brown1','deepskyblue4','darkseagreen3'))+ 
#   theme(strip.text.x = element_text(angle = 90),
#         axis.text.y = element_blank(),
#         axis.title.y = element_blank())+
#   scale_y_discrete(labels=function(x) str_remove(x, "[A-Za-z0-9]+_[2-9]|_[0-9]+"))
# p1
# 
# g1<-ggplot(res_adxloc_gwasf[chr%in%1:9])+
#   geom_point(aes(y=locus_gene_3,x=gwas_short2,
#                  col=gwas_sig),size=1,shape=15)+
#   facet_grid(chr~'',scales = 'free',space = 'free')+
#   # scale_size(range = c(0.5,2))+
#   theme_minimal()+
#   scale_x_discrete(guide = guide_axis(angle = 90))+
#   scale_color_manual(values = c('grey','bisque3','orange3','brown4'))+
#   theme(axis.text.x = element_markdown(size = 7),
#         axis.text.y = element_markdown(size=7),axis.title.y =element_blank(),
#         strip.text = element_blank(),   # remove facet text
#         strip.background = element_blank() 
#   )+
#   scale_y_discrete(labels=function(x) str_remove(x, "[A-Za-z0-9]+_[2-9]|_[0-9]+"))
# 
# # P1<-g1+p1+plot_layout(guides = 'collect',widths = c(1.5,6))&
# #   theme(plot.margin = margin(0, 0, 0, 0)) 
# # P1
# 
# 
# #chr10-22
# p2<-ggplot(res_adxlocge_cont_topf[chr%in%10:22][!(APOE_region.gene)&(top_gene|top_confidence%in%c("CL1",'CL2','CL3'))])+
#   geom_point(aes(y=locus_gene_3,x=context_group,
#                  size=n_study_group.locus,
#                  col=confidence_cat_group))+
#   facet_grid(chr~'',scales = 'free',space = 'free')+
#   scale_size(range = c(1.5,5),breaks = c(2,7,12))+theme_minimal()+
#   scale_x_discrete(guide = guide_axis(angle = 90))+
#   theme(strip.text.x = element_text(angle = 90))+
#   labs(size='# datasets',col='Confidence level')+
#   scale_color_manual(values = c('brown1','deepskyblue4','darkseagreen3'))+ 
#   theme(strip.text.x = element_text(angle = 90),
#         axis.text.y = element_blank(),
#         axis.title.y = element_blank())+
#   scale_y_discrete(labels=function(x) str_remove(x, "[A-Za-z0-9]+_[2-9]|_[0-9]+"))
# p2
# 
# g2<-ggplot(res_adxloc_gwasf[chr%in%10:22])+
#   geom_point(aes(y=locus_gene_3,x=gwas_short2,
#                  col=gwas_sig),size=1,shape=15)+
#   facet_grid(chr~'',scales = 'free',space = 'free')+
#   # scale_size(range = c(0.5,2))+
#   theme_minimal()+
#   scale_x_discrete(guide = guide_axis(angle = 90))+
#   scale_color_manual(values = c('grey','bisque3','orange3','brown4'))+
#   theme(axis.text.x = element_markdown(size = 7),
#         axis.text.y = element_markdown(size=7),axis.title.y =element_blank(),
#         strip.text = element_blank(),   # remove facet text
#         strip.background = element_blank() 
#   )+
#   scale_y_discrete(labels=function(x) str_remove(x, "[A-Za-z0-9]+_[2-9]|_[0-9]+"))
# 
# 
# 
# g1+p1+g2+p2+plot_layout(guides = 'collect',widths = c(1.5,6,1.5,6))&
#   theme(plot.margin = margin(0, 0, 0, 0)) 
# ggsave(fp(out,'ADloci_xQTL_summary_genome_wide_signif_top_gene_per_locus_3categories_gwas_plot_top_CL_genes_all_causal_top1_associated.pdf'),height = 9,width = 10)
# ggsave(fp(out,'ADloci_xQTL_summary_genome_wide_signif_top_gene_per_locus_3categories_gwas_plot_top_CL_genes_all_causal_top1_associated.png'),height = 9,width = 10)
# 
# #add the locus number
# 
# g1<-g1+theme(axis.text.y = element_blank())
# g2<-g2+theme(axis.text.y = element_blank())
# 
# # display.brewer.all(type = "qual")
# 
# palette=brewer.pal(8, "Pastel1")
# 
# 
# cols1 <- data.table(locus_index=unique(res_adxlocge_cont_topf[chr%in%1:9][!(APOE_region.gene)&(top_gene|top_confidence%in%c("CL1",'CL2','CL3'))]$locus_index),
#                       color=palette)
# 
# l1<-ggplot(unique(res_adxlocge_cont_topf[chr%in%1:9][!(APOE_region.gene)&(top_gene|top_confidence%in%c("CL1",'CL2','CL3'))],
#                   by=c('locus_index','gene_name')),
#            aes(x=0,y=locus_gene_3))+
#   geom_tile(aes(fill=as.factor(locus_index)))+scale_fill_manual(values = setNames(cols1$color,as.character(cols1$locus_index)))+
#   geom_text(aes(label=locus_index),size=2)+
#   facet_grid(chr~'',scales = 'free',space = 'free')+theme_minimal()+
#   theme(axis.text.y = element_markdown(size=7),
#         axis.title.y = element_blank(), 
#         strip.text = element_blank(),   # remove facet text
#         strip.background = element_blank(),
#         axis.text.x=element_blank(),axis.title.x = element_blank(),
#         legend.position = "none")+
#   scale_y_discrete(labels=function(x) str_remove(x, "[A-Za-z0-9]+_[2-9]|_[0-9]+"))
# 
# 
# cols2 <- data.table(locus_index=unique(res_adxlocge_cont_topf[chr%in%10:22][!(APOE_region.gene)&(top_gene|top_confidence%in%c("CL1",'CL2','CL3'))]$locus_index),
#                       color=palette)
# 
# l2<-ggplot(unique(res_adxlocge_cont_topf[chr%in%10:22][!(APOE_region.gene)&(top_gene|top_confidence%in%c("CL1",'CL2','CL3'))],
#                   by=c('locus_index','gene_name')),
#            aes(x=0,y=locus_gene_3))+
#   geom_tile(aes(fill=as.factor(locus_index)))+scale_fill_manual(values = setNames(cols2$color,as.character(cols2$locus_index)))+
#   geom_text(aes(label=locus_index),size=2)+
#   facet_grid(chr~'',scales = 'free',space = 'free')+theme_minimal()+
#   theme(axis.text.y = element_markdown(size=7),
#         axis.title.y = element_blank(), strip.text = element_blank(),   # remove facet text
#         strip.background = element_blank(),
#         axis.text.x=element_blank(),axis.title.x = element_blank(),
#         legend.position = "none")+
#   scale_y_discrete(labels=function(x) str_remove(x, "[A-Za-z0-9]+_[2-9]|_[0-9]+"))
# 
# 
# l1+g1+p1+l2+g2+p2+plot_layout(guides = 'collect',widths = c(0.5,1.5,6,0.5,1.5,6))&
#   theme(plot.margin = margin(0, 0, 0, 0)) 
# 
# 
# #with all genes
# #filter the gwas and annot with the gene too
# res_adxlocge_cont_topf[,locus_gene_2:=str_remove(locus_gene_2,'\\*')]
# 
# res_adxlocge_cont_topf[,locus_gene_3:=factor(locus_gene_2,levels = unique(locus_gene_2[order(-chr,-tss,-as.numeric(str_extract(locus_gene_2,'_[0-9+]')|>str_remove('_')))]))]
# 
# res_adxloc_gwasf<-merge(res_adxloc_gwas,unique(res_adxlocge_cont_topf[!(APOE_region.gene)][,.(locus_index,locus_gene_3,locus_gene)]),by='locus_index')
# 
# #chr1-9
# p1<-ggplot(res_adxlocge_cont_topf[chr%in%1:9][!(APOE_region.gene)])+
#   geom_point(aes(y=locus_gene_3,x=context_group,
#                  size=n_study_group.locus,
#                  col=confidence_cat_group))+
#   facet_grid(chr~'',scales = 'free',space = 'free')+
#   scale_size(range = c(1.5,5),breaks = c(2,7,12))+theme_minimal()+
#   scale_x_discrete(guide = guide_axis(angle = 90))+
#   theme(strip.text.x = element_text(angle = 90))+
#   labs(size='# datasets',col='Confidence level')+
#   scale_color_manual(values = c('brown1','deepskyblue4','darkseagreen3'))+ 
#   theme(strip.text.x = element_text(angle = 90),
#         axis.text.y = element_blank(),
#         axis.title.y = element_blank())+
#   scale_y_discrete(labels=function(x) str_remove(x, "[A-Za-z0-9]+_[2-9]|_[0-9]+"))
# 
# g1<-ggplot(res_adxloc_gwasf[chr%in%1:9])+
#   geom_point(aes(y=locus_gene_3,x=gwas_short2,
#                  col=gwas_sig),size=1,shape=15)+
#   facet_grid(chr~'',scales = 'free',space = 'free')+
#   # scale_size(range = c(0.5,2))+
#   theme_minimal()+
#   scale_x_discrete(guide = guide_axis(angle = 90))+
#   scale_color_manual(values = c('grey','bisque3','orange3','brown4'))+
#   theme(axis.text.x = element_markdown(size = 7),
#         axis.text.y = element_markdown(size=7),axis.title.y =element_blank(),
#         strip.text = element_blank(),   # remove facet text
#         strip.background = element_blank() 
#   )+
#   scale_y_discrete(labels=function(x) str_remove(x, "[A-Za-z0-9]+_[2-9]|_[0-9]+"))
# 
# 
# #chr10-22
# p2<-ggplot(res_adxlocge_cont_topf[chr%in%10:22][!(APOE_region.gene)])+
#   geom_point(aes(y=locus_gene_3,x=context_group,
#                  size=n_study_group.locus,
#                  col=confidence_cat_group))+
#   facet_grid(chr~'',scales = 'free',space = 'free')+
#   scale_size(range = c(1.5,5),breaks = c(2,7,12))+theme_minimal()+
#   scale_x_discrete(guide = guide_axis(angle = 90))+
#   theme(strip.text.x = element_text(angle = 90))+
#   labs(size='# datasets',col='Confidence level')+
#   scale_color_manual(values = c('brown1','deepskyblue4','darkseagreen3'))+ 
#   theme(strip.text.x = element_text(angle = 90),
#         axis.text.y = element_blank(),
#         axis.title.y = element_blank())+
#   scale_y_discrete(labels=function(x) str_remove(x, "[A-Za-z0-9]+_[2-9]|_[0-9]+"))
# 
# g2<-ggplot(res_adxloc_gwasf[chr%in%10:22])+
#   geom_point(aes(y=locus_gene_3,x=gwas_short2,
#                  col=gwas_sig),size=1,shape=15)+
#   facet_grid(chr~'',scales = 'free',space = 'free')+
#   # scale_size(range = c(0.5,2))+
#   theme_minimal()+
#   scale_x_discrete(guide = guide_axis(angle = 90))+
#   scale_color_manual(values = c('grey','bisque3','orange3','brown4'))+
#   theme(axis.text.x = element_markdown(size = 7),
#         axis.text.y = element_markdown(size=7),axis.title.y =element_blank(),
#         strip.text = element_blank(),   # remove facet text
#         strip.background = element_blank() 
#   )+
#   scale_y_discrete(labels=function(x) str_remove(x, "[A-Za-z0-9]+_[2-9]|_[0-9]+"))
# 
# 
# #add the tile of the locus number
# g1<-g1+theme(axis.text.y = element_blank())
# g2<-g2+theme(axis.text.y = element_blank())
# 
# palette=brewer.pal(8, "Pastel1")
# cols1 <- data.table(locus_index=unique(res_adxlocge_cont_topf[chr%in%1:9][!(APOE_region.gene)]$locus_index),
#                     color=palette)
# 
# l1<-ggplot(unique(res_adxlocge_cont_topf[chr%in%1:9][!(APOE_region.gene)],
#                   by=c('locus_index','gene_name')),
#            aes(x=0,y=locus_gene_3))+
#   geom_tile(aes(fill=as.factor(locus_index)))+scale_fill_manual(values = setNames(cols1$color,as.character(cols1$locus_index)))+
#   geom_text(aes(label=locus_index),size=2)+
#   facet_grid(chr~'',scales = 'free',space = 'free')+theme_minimal()+
#   theme(axis.text.y = element_markdown(size=6),
#         axis.title.y = element_blank(), 
#         strip.text = element_blank(),   # remove facet text
#         strip.background = element_blank(),
#         axis.text.x=element_blank(),axis.title.x = element_blank(),
#         legend.position = "none")+
#   scale_y_discrete(labels=function(x) str_remove(x, "[A-Za-z0-9]+_[2-9]|_[0-9]+"))
# 
# 
# cols2 <- data.table(locus_index=unique(res_adxlocge_cont_topf[chr%in%10:22][!(APOE_region.gene)]$locus_index),
#                     color=palette)
# 
# l2<-ggplot(unique(res_adxlocge_cont_topf[chr%in%10:22][!(APOE_region.gene)],
#                   by=c('locus_index','gene_name')),
#            aes(x=0,y=locus_gene_3))+
#   geom_tile(aes(fill=as.factor(locus_index)))+scale_fill_manual(values = setNames(cols2$color,as.character(cols2$locus_index)))+
#   geom_text(aes(label=locus_index),size=2)+
#   facet_grid(chr~'',scales = 'free',space = 'free')+theme_minimal()+
#   theme(axis.text.y = element_markdown(size=6),
#         axis.title.y = element_blank(), strip.text = element_blank(),   # remove facet text
#         strip.background = element_blank(),
#         axis.text.x=element_blank(),axis.title.x = element_blank(),
#         legend.position = "none")+
#   scale_y_discrete(labels=function(x) str_remove(x, "[A-Za-z0-9]+_[2-9]|_[0-9]+"))
# 
# 
# l1+g1+p1+l2+g2+p2+plot_layout(guides = 'collect',widths = c(0.5,1.5,6,0.5,1.5,6))&
#   theme(plot.margin = margin(0, 0, 0, 0)) 
# 
# 
# 

#with all  genes####

#get GWAS signals
res_ad<-fread(gwas_table)
res_ad<-merge(res_ad,unique(res_adx[,.(locus_index,ADlocusID)]))

res_ad[,gwas_pvalue:=getPval(gwas_zscore)]

res_adxloc_gwas<-unique(res_ad[order(gwas_pvalue)][gwas_source!=''],by=c('locus_index','gwas_source'))
res_adxloc_gwas[,chr:=str_extract(chr,'[0-9]+')|>as.numeric()]

res_adxloc_gwas[,gwas_pvalue10:=ifelse(gwas_pvalue>1e-10,gwas_pvalue,1e-10)]
res_adxloc_gwas[,gwas_sig:=ifelse(gwas_pvalue<5e-8,'p<5e-8',ifelse(gwas_pvalue<1e-6,'p<1e-6',ifelse(gwas_pvalue<1e-5,'p<1e-5',
                                                                                                    'p>1e-5 but coloc')))]
res_adxloc_gwas[is.na(gwas_sig),gwas_sig:='p>1e-5 but coloc']
res_adxloc_gwas[,gwas_sig:=factor(gwas_sig,levels=c('p>1e-5 but coloc','p<1e-5','p<1e-6','p<5e-8'))]


res_adxloc_gwas[,gwas_short:=str_extract(gwas_source,'Bellenguez|Jansen|Kunkle|Wightman')]
table(res_adxloc_gwas$gwas_short)
res_adxloc_gwas[str_detect(gwas_source,'EADB|EADI|UKB|23andMe'),gwas_short:=paste(gwas_short,str_extract(gwas_source,'EADB|EADI|UKB|23andMe'),sep='_')]
res_adxloc_gwas[gwas_short=='Wightman_23andMe',gwas_short:='Wightman_no23andMe']
res_adxloc_gwas[gwas_short=='Wightman_UKB',gwas_short:='Wightman_no23AndUKB']

#order by n AD cases (non AD by proxi first) 
res_adxloc_gwas[,ad_by_proxi:=!gwas_source%in%c('AD_Bellenguez_EADB_2022','AD_Bellenguez_EADI_2022',
                                                'AD_Wightman_ExcludingUKBand23andME_2021',
                                                'AD_Kunkle_Stage1_2019'),by='locus_index']

gwmt<-fread(gwas_ncases_controls)
res_adxloc_gwas<-merge(res_adxloc_gwas,gwmt[,.(gwas_source=study_id,n_case,n_control)])
#bold if not from  AD by proxy 
res_adxloc_gwas[,gwas_short2:=ifelse(!ad_by_proxi,paste0("<b>",gwas_short,'</b>'),
                                     gwas_short)]

fwrite(res_adxloc_gwas,fp(out,'res_summ_all_gwas_by_locus.csv.gz'))

res_adxloc_gwas[,gwas_short2:=factor(gwas_short2,levels = unique(gwas_short2[order(ad_by_proxi,-n_case)]))]
levels(res_adxloc_gwas$gwas_short2)



# group some contexts
res_adxloc<-merge(res_adxloc,context_groups,by='context_short')

res_adxloc<-res_adxloc[!is.na(context_group)]
# summarize using confidence score  and nstudy
res_adxloc[,n_study_group:=length(unique(context)),by=.(gene_name,context_group)]
res_adxloc[,confidence_lvl_group:=sort(confidence_lvl)[1],by=.(gene_name,context_group)]
res_adxloc[,n_loci_gene:=unique(locus_index)|>length(),by=.(gene_name)]

res_adxloc[,n_study_group.locus:=length(unique(context)),by=.(gene_name,locus_index,context_group)]
res_adxloc[,confidence_lvl_group.locus:=sort(confidence_lvl)[1],by=.(gene_name,locus_index,context_group)]


fwrite(res_adxloc,fp(out,'res_summ_all_genes_by_locus.csv.gz'))

res_adxloc<-fread(fp(out,'res_summ_all_genes_by_locus.csv.gz'))



setdiff(res_adxloc$context_group,contexts_order)
res_adxloc[,context_group:=factor(context_group,levels = contexts_order)]

#genome wide signifcant hits####

res_adxlocf<-res_adxloc[(genomewide_sig_gene)]

unique(res_adxlocf$locus_index)|>length() #103
unique(res_adxlocf$gene_name)|>length()#79> 116 >134>292

unique(res_adxlocf$locus_gene)|>length()#98> 138 >156>359
unique(res_adxlocf[!(APOE_region.gene)]$locus_gene)|>length()#86>122>140>334

# extract each locus
res_adxlocf_cont_top<-unique(res_adxlocf[order(confidence_lvl_group.locus)],by=c('gene_name','context_group','locus_index'))[gene_name%in%gene_name[confidence_lvl_group.locus%in%c('CL1','CL2','CL3','CL4','CL5')]]


#sep causal vs correlated and remove all C6
#C6 is dropped
res_adxlocf_cont_topf<-res_adxlocf_cont_top[confidence_lvl_group.locus!='CL6']
#sep in 3categories
res_adxlocf_cont_topf[,confidence_cat_group:=ifelse(confidence_lvl_group.locus%in%c('CL1',"CL2"),'Putative causal (CL1, CL2)',ifelse(confidence_lvl_group.locus=='CL3','Putative causal (CL3)','Associated (CL4, CL5)'))]
res_adxlocf_cont_topf[,confidence_cat_group:=factor(confidence_cat_group,levels = c('Putative causal (CL1, CL2)','Putative causal (CL3)','Associated (CL4, CL5)'))]

#annot genes
res_adxlocf_cont_topf[,locus_gene_2:=paste(gene_name,match(locus_gene,
                                                            unique(locus_gene[order(chr,tss,min_pval)])),sep='_'),
                       by=.(gene_name)]

res_adxlocf_cont_topf[,locus_gene_2:=paste0("<i>",locus_gene_2[1],'</i>'),
                       by='locus_gene']

#add if from  AD by proxy only
res_adxlocf_cont_topf[,locus_gene_2:=ifelse(!all(only_by_proxi,na.rm = T),paste0("<b>",locus_gene_2[1],'</b>'),
                                             locus_gene_2[1]),
                       by='locus_gene']

#order per gene tss and locus position
res_adxlocf_cont_topf[,locus_gene_3:=factor(locus_gene_2,levels = unique(locus_gene_2[order(chr,-tss,-min_pval)]))]

ndatasets_range<-range(res_adxlocf_cont_topf$n_study_group.locus)
ps<-lapply(list(1:7,8:15,16:17,18:22), function(chrs){
  SummPanel(res_adxlocf_cont_topf[chr%in%chrs],res_adxloc_gwas,ndatasets_range=ndatasets_range)
})|>unlist(recursive = FALSE)

wrap_plots(ps)+plot_layout(guides = 'collect',widths = rep(c(0.5,1.5,6),4))&
  theme(plot.margin = margin(0, 0, 0, 0)) 


#excluding APOE MAPT
splits<-list(1:7,8:15,16:22)

ps<-lapply(splits, function(chrs){
  SummPanel(res_adxlocf_cont_topf[chr%in%chrs][!(MAPT_region.locus|APOE_region.gene)],res_adxloc_gwas,ndatasets_range=ndatasets_range)
})|>unlist(recursive = FALSE)

wrap_plots(ps)+plot_layout(guides = 'collect',widths = rep(c(0.5,1.5,6),length(splits)))&
  theme(plot.margin = margin(0, 0, 0, 0)) 



#suggestive signal####


res_adxlocf<-res_adxloc[!(genomewide_sig_gene)][gene_name%in%gene_name[confidence_lvl_group%in%c('CL1','CL2','CL3','CL4')&n_gwas_locus>2&n_study_group.locus>=2]]

res_adxlocf_cont<-unique(res_adxlocf,by=c('gene_name','context_group','locus_index'))


#sep causal vs correlated 
#add if frpm, APOE region, from AD by proxy only, or from more than one locus

res_adxlocf_cont[,confidence_cat_group:=ifelse(confidence_lvl_group.locus%in%c('CL1',"CL2"),'Putative causal (CL1, CL2)',ifelse(confidence_lvl_group.locus=='CL3','Putative causal (CL3)','Associated (CL4, CL5)'))]
res_adxlocf_cont[,confidence_cat_group:=factor(confidence_cat_group,levels = c('Putative causal (CL1, CL2)','Putative causal (CL3)','Associated (CL4, CL5)'))]

#annot genes
res_adxlocf_cont[,locus_gene_2:=paste(gene_name,match(locus_gene,
                                                           unique(locus_gene[order(chr,tss,min_pval)])),sep='_'),
                      by=.(gene_name)]

res_adxlocf_cont[,locus_gene_2:=paste0("<i>",locus_gene_2[1],'</i>'),
                      by='locus_gene']

#add if from  AD by proxy only
res_adxlocf_cont[,locus_gene_2:=ifelse(!all(only_by_proxi,na.rm = T),paste0("<b>",locus_gene_2[1],'</b>'),
                                            locus_gene_2[1]),
                      by='locus_gene']

#order per gene tss and locus position
res_adxlocf_cont[,locus_gene_3:=factor(locus_gene_2,levels = unique(locus_gene_2[order(chr,-tss,-min_pval)]))]


ps<-SummPanel(res_adxlocf_cont,res_adxloc_gwas)

wrap_plots(ps)+plot_layout(guides = 'collect',widths = c(0.5,1.5,6))&
  theme(plot.margin = margin(0, 0, 0, 0)) 


ggsave(fp(out,'ADloci_xQTL_summary_suggestive_loci_top_gene_per_locus_3categories.png'),height = 5,width = 6)


setdiff(res_adxlocge_cont_top2f$gene_name,res_adxlocge_cont_topf$gene_name)
#v "CCNT2"   "ICA1L"   "BLNK"    "PLEKHA1" "RITA1"   "CTSH" 



