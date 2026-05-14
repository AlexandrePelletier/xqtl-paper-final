source('../../utils/r_utils.R')
out<-'main_text/5_AD_xQTL_genes_cis_trans/staging/APOE/'
dir.create(out)
#APOE4 independant APOE cs are associated with neuropatho in ROSMAP ?

res_indep_cs<-readRDS('main_text/5_AD_xQTL_genes_cis_trans/staging/APOE/xqtl_only_APOE_all_cohorts_independent_sets_after_imputation_no_sQTL_only.rds')
length(res_indep_cs$variant_id)
indep_cs<-rbindlist(lapply(names(res_indep_cs$variant_id),function(n)data.table(cs_name=n,variant_id=res_indep_cs$variant_id[[n]],
                                                                                gwas_source=res_indep_cs$ad_gwas[[which(names(res_indep_cs$variant_id)==n)]],
                                                                                molecular_outcome=res_indep_cs$outcome[[n]])))
indep_cs[,cs_id:=paste0('set_',str_extract(cs_name,'[0-9]+$'))]
indep_cs[,chr:=seqid(variant_id)]
indep_cs[,pos:=pos(variant_id)]
indep_cs[,variant_id2:=str_replace_all(variant_id,':','_')|>str_replace('_',':')]

#extract genotype for those CS
fwrite(indep_cs[order(chr,pos)][,.(variant_id2)],fp(out,'apoe4_indep_snps.tsv'),sep = '\t',col.names=FALSE)

plink_file<-'/projectnb/tcwlab-adsp/member/adpelle1/projects/fungen-xqtl/ref-data/ROSMAP/ROSMAP_NIA_geno/ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc'

cmd=paste('plink --bfile',tools::file_path_sans_ext(plink_file),
          '--extract',fp(out,'apoe4_indep_snps.tsv'),
          '--recode vcf',
          '--out',fp(out,'ROSMAP_geno_apoe4_indep_snps'))

system(cmd)

#format the genotype
donors_apoecs<-fread(fp(out,'ROSMAP_geno_apoe4_indep_snps.vcf'))
dim(donors_apoecs)#216 1162

donors_apoecs<-melt(donors_apoecs[,-c(6:9)],
                    id.vars = c('#CHROM','POS','ID','REF','ALT'),
                    variable.name = 'IID',value.name = 'geno')
donors_apoecs[,dose:=str_count(geno,'1')]

donors_apoecs[,variant_id2:=ID]
donors_apoecs<-merge(donors_apoecs,indep_cs,by='variant_id2')


#add the clinical covariates
clin<-fread('/projectnb/tcwlab/ShareSpace/ROSMAP_Neuropatho/ROSMAP_clinical.csv')
bios<-fread('/projectnb/tcwlab/ShareSpace/ROSMAP_Neuropatho/ROSMAP_biospecimen_metadata.csv') #for WGS samples ID - IndividualID link
clinw<-merge(clin,bios[assay=='wholeGenomeSeq'][,.(individualID,specimenID)],by='individualID')
complete<-fread('/projectnb/tcwlab/ShareSpace/ROSMAP_Neuropatho/ROSMAP_xqtl_complete_samples_covariates_sex_death_pmi_study.csv')
setnames(complete,'sample_id','specimenID')
clinw<-merge(clinw[,-c('age_death','pmi','msex')],complete,by='specimenID')

donors_apoecs[,specimenID:=str_remove(IID,'0_')]

donors_apoecs<-merge(donors_apoecs,clinw,by='specimenID')

unique(donors_apoecs$individualID)|>length()#1116

#add neuropatho covs
#we test only a subset of the covariates
neuropath<-fread('/projectnb/tcwlab/ShareSpace/ROSMAP_Neuropatho/dataset_707_basic_02-10-2022.n3711.clean.plink.n2409.28FEB2022.txt')
dim(neuropath)
neuropathf<-neuropath[,.SD,.SDcols = c('projid','amyloid_sqrt','tangles_sqrt','cogng_path_slope','plaq_d','plaq_n_sqrt','nft_sqrt','tdp_cs_6reg','tdp_dn_6reg','tdp_st4','caa_4gp','caa_neo4')]

#add those phenotype to the genotypes 
donors_apoecs<-merge(donors_apoecs,neuropathf,all.x = TRUE,by='projid')

#neuropatho analysis correcting for APOE4 /2 gentoyep
donors_apoecs[,apoe4_dose:=str_count(apoe_genotype,'4')]
donors_apoecs[,apoe2_dose:=str_count(apoe_genotype,'2')]

#prep/reformat the phenotypes
donors_apoecs[,cogdx_binary:=ifelse(cogdx%in%c(4,5),'AD',ifelse(cogdx%in%1:3,'No AD',NA))] #binarize the cognitive status

donors_apoecs[,cogdx_merge:=ifelse(cogdx%in%c(4,5),'AD',ifelse(cogdx%in%2:3,'MCI',ifelse(cogdx==1,'NL',NA)))] # three level cognitive status
donors_apoecs[,cogdx_merge:=factor(cogdx_merge,levels = c('NL','MCI','AD'))]
donors_apoecs[cogdx!=6,cogdx_mergenum:=ifelse(cogdx_merge=='MCI',2,ifelse(cogdx_merge=='AD',3,1))]

donors_apoecs[,age_first_ad_dx_num:=ifelse(age_first_ad_dx=='90+',91,as.numeric(age_first_ad_dx))]
donors_apoecs[,other_dementia:=cogdx==6]
donors_apoecs[,other_CI:=cogdx%in%c(3,5)]
donors_apoecs[cogdx%in%c(1,4,5),ADvsNL:=cogdx%in%c(4,5)]
donors_apoecs[cogdx%in%c(1,2,3),MCIvsNL:=cogdx%in%c(2,3)]

donors_apoecs[,ceradsc_pos:=-ceradsc]
donors_apoecs[,braaksc06:=ifelse(braaksc==0,1,ifelse(braaksc==6,5,braaksc))]

donors_apoecs[,age_first_ad_dx_num:=ifelse(age_first_ad_dx=='90+',91,as.numeric(age_first_ad_dx))]
donors_apoecs[(cogdx_binary=='AD'),age_first_clinad_dx:=age_first_ad_dx_num]
donors_apoecs[(cogdx_binary=='AD'),cts_mmse30_lv_clinad:=cts_mmse30_lv]

donors_apoecs[,age_at_visit_num:=ifelse(age_at_visit_max=='90+',
                                      91,as.numeric(age_at_visit_max))]
donors_apoecs[dcfdx_lv!=6,dcfdx_lv_no6:=dcfdx_lv]

donors_apoecs[,apoe2_r:=cor(dose,apoe2_dose,use = 'pairwise.complete.obs'),by='ID']
donors_apoecs[,apoe4_r:=cor(dose,apoe4_dose,use = 'pairwise.complete.obs'),by='ID']

fwrite(donors_apoecs,fp(out,'ROSMAP_geno_apoecs_adpathology.csv.gz'))


#outcome to test with Lm
covs_lin<-c('ceradsc_pos','braaksc06',
            'dcfdx_lv_no6','cts_mmse30_lv','age_death',
            'cogdx_mergenum','apoe4_dose','apoe2_dose',
            'cts_mmse30_lv_clinad','age_first_ad_dx_num',
            'age_first_clinad_dx','amyloid_sqrt','tangles_sqrt',
            'cogng_path_slope',
            'plaq_d','plaq_n_sqrt','nft_sqrt','tdp_cs_6reg','tdp_dn_6reg','tdp_st4','caa_4gp','caa_neo4'
)

#outcome to test with logistic regression
covs_logis<-c('MCIvsNL','ADvsNL')

#outcome that are measured before death
atlastvisit=c('cts_mmse30_lv','cts_mmse30_lv_clinad')#outcome that need to be adjusted for age at last visit

beforedeath<-c('age_death','cts_mmse30_lv','cts_mmse30_lv_clinad','age_first_ad_dx_num','age_first_clinad_dx','cogng_path_slope')

#Adjusting for both APOE4 and APOE2####
apoecs_clin<-rbindlist(lapply(covs_lin,function(outc){
  message(outc)
  if(outc%in%atlastvisit){
    donors_apoecs[,summary(lm(unlist(.SD)~dose+msex+educ+age_at_visit_num+apoe4_dose+apoe2_dose))$coefficients|>data.table(keep.rownames = 'cov'),
                  by='variant_id2',.SDcols=outc][cov=='dose'][,outcome:=outc]

  }else if(outc%in%beforedeath){
    donors_apoecs[,summary(lm(unlist(.SD)~dose+msex+educ+apoe4_dose+apoe2_dose))$coefficients|>data.table(keep.rownames = 'cov'),
                  by='variant_id2',.SDcols=outc][cov=='dose'][,outcome:=outc]
  }else{

    donors_apoecs[,summary(lm(unlist(.SD)~dose+msex+educ+age_death+apoe4_dose+apoe2_dose))$coefficients|>data.table(keep.rownames = 'cov'),
                  by='variant_id2',.SDcols=outc][cov=='dose'][,outcome:=outc]
  }
  
}))
apoecs_clin[,zscore:=`t value`]
apoecs_clin[,pval:=`Pr(>|t|)`]

apoecs_clogit<-rbindlist(lapply(covs_logis,function(outc){
  message(outc)
  donors_apoecs[,summary(glm(unlist(.SD)~dose+msex+educ+age_death+apoe4_dose+apoe2_dose,family = binomial))$coefficients|>data.table(keep.rownames = 'cov'),
              by='variant_id2',.SDcols=outc][cov=='dose'][,outcome:=outc]
  
}))
apoecs_clogit[,zscore:=`z value`]
apoecs_clogit[,pval:=`Pr(>|z|)`]

apoecs_clin<-rbind(apoecs_clin,apoecs_clogit,fill=TRUE)

apoecs_clin<-merge(apoecs_clin,indep_cs,by='variant_id2')
apoecs_clin[,padj:=p.adjust(pval,n = length(unique(indep_cs$cs_name))),by=c('outcome','variant_id2')]
apoecs_clin[padj<0.05]


apoecs_clin[!outcome%in%c('apoe4_dose','apoe2_dose')][padj<0.05]

apoecs_clin[padj<0.05][!outcome%in%c('apoe4_dose','apoe2_dose')]$cs_id  |>unique()|>length()#3/19 

fwrite(apoecs_clin,fp(out,'res_lm_apoecs_adpathology_apoe4_2_correction.csv.gz'))


unique(apoecs_clin[padj<0.05][!outcome%in%c('apoe4_dose','apoe2_dose')][order(cs_id,outcome,padj)][,.(cs_id,outcome,padj,gwas_source,molecular_outcome)],by=c('cs_id','outcome'))

#     cs_id             outcome        padj                                                                                                              gwas_source
#    <char>              <char>       <num>                                                                                                                   <char>
# 1: set_27             caa_4gp 0.004708944 AD_Kunkle_Stage1_2019; AD_Wightman_Full_2021; AD_Wightman_Excluding23andMe_2021; AD_Wightman_ExcludingUKBand23andME_2021
# 2: set_27            caa_neo4 0.007830511 AD_Kunkle_Stage1_2019; AD_Wightman_Full_2021; AD_Wightman_Excluding23andMe_2021; AD_Wightman_ExcludingUKBand23andME_2021
# 3: set_27             tdp_st4 0.042915291 AD_Kunkle_Stage1_2019; AD_Wightman_Full_2021; AD_Wightman_Excluding23andMe_2021; AD_Wightman_ExcludingUKBand23andME_2021
# 4: set_31 age_first_ad_dx_num 0.020960086                                                                                                            AD_Bellenguez
# 5: set_57        tangles_sqrt 0.003434581                           AD_Bellenguez; AD_Kunkle_Stage1_2019; AD_Wightman_Full_2021; AD_Wightman_Excluding23andMe_2021
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     molecular_outcome
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <char>
# 1: AC_DeJager_eQTL_ENSG00000224916; AC_DeJager_eQTL_ENSG00000234906; BM_10_MSBB_eQTL_ENSG00000104853; Metabrain_Basalganglia_chr19_41840000_47960000_ENSG00000234906; Metabrain_Cerebellum_chr19_41840000_47960000_ENSG00000234906; Metabrain_Cortex_chr19_41840000_47960000_ENSG00000104853; Metabrain_Cortex_chr19_41840000_47960000_ENSG00000234906; Mic_DeJager_eQTL_ENSG00000130208; Mic_DeJager_eQTL_ENSG00000130203; Mic_mega_eQTL_ENSG00000130203; STARNET_eQTL_Mac_ENSG00000130208; STARNET_eQTL_Mac_ENSG00000130203; STARNET_eQTL_Mac_ENSG00000267467; STARNET_eQTL_Mac_ENSG00000224916; STARNET_eQTL_Mac_ENSG00000234906; DeJager_Mic_ENSG00000130203; Kellis_Mic_ENSG00000130203; DeJager_Mic_ENSG00000130208; Kellis_Mic_ENSG00000130208
# 2: AC_DeJager_eQTL_ENSG00000224916; AC_DeJager_eQTL_ENSG00000234906; BM_10_MSBB_eQTL_ENSG00000104853; Metabrain_Basalganglia_chr19_41840000_47960000_ENSG00000234906; Metabrain_Cerebellum_chr19_41840000_47960000_ENSG00000234906; Metabrain_Cortex_chr19_41840000_47960000_ENSG00000104853; Metabrain_Cortex_chr19_41840000_47960000_ENSG00000234906; Mic_DeJager_eQTL_ENSG00000130208; Mic_DeJager_eQTL_ENSG00000130203; Mic_mega_eQTL_ENSG00000130203; STARNET_eQTL_Mac_ENSG00000130208; STARNET_eQTL_Mac_ENSG00000130203; STARNET_eQTL_Mac_ENSG00000267467; STARNET_eQTL_Mac_ENSG00000224916; STARNET_eQTL_Mac_ENSG00000234906; DeJager_Mic_ENSG00000130203; Kellis_Mic_ENSG00000130203; DeJager_Mic_ENSG00000130208; Kellis_Mic_ENSG00000130208
# 3: AC_DeJager_eQTL_ENSG00000224916; AC_DeJager_eQTL_ENSG00000234906; BM_10_MSBB_eQTL_ENSG00000104853; Metabrain_Basalganglia_chr19_41840000_47960000_ENSG00000234906; Metabrain_Cerebellum_chr19_41840000_47960000_ENSG00000234906; Metabrain_Cortex_chr19_41840000_47960000_ENSG00000104853; Metabrain_Cortex_chr19_41840000_47960000_ENSG00000234906; Mic_DeJager_eQTL_ENSG00000130208; Mic_DeJager_eQTL_ENSG00000130203; Mic_mega_eQTL_ENSG00000130203; STARNET_eQTL_Mac_ENSG00000130208; STARNET_eQTL_Mac_ENSG00000130203; STARNET_eQTL_Mac_ENSG00000267467; STARNET_eQTL_Mac_ENSG00000224916; STARNET_eQTL_Mac_ENSG00000234906; DeJager_Mic_ENSG00000130203; Kellis_Mic_ENSG00000130203; DeJager_Mic_ENSG00000130208; Kellis_Mic_ENSG00000130208
# 4:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               AC_DeJager_eQTL_ENSG00000189114; DLPFC_DeJager_eQTL_ENSG00000189114; Metabrain_Cortex_chr19_41840000_47960000_ENSG00000189114; PCC_DeJager_eQTL_ENSG00000189114; ROSMAP_AC_ENSG00000189114; ROSMAP_DLPFC_ENSG00000189114; ROSMAP_PCC_ENSG00000189114
# 5:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Ast_10_Kellis_eQTL_ENSG00000160007


#set27 is the caa set of Fujita et al?
apoecs_clin[pos(variant_id2)==44946027]#yes



#adjusting for APOE4 only####
apoecs_clin2<-rbindlist(lapply(covs_lin,function(outc){
  message(outc)
  if(outc%in%atlastvisit){
    donors_apoecs[,summary(lm(unlist(.SD)~dose+msex+educ+age_at_visit_num+apoe4_dose))$coefficients|>data.table(keep.rownames = 'cov'),
                  by='variant_id2',.SDcols=outc][cov=='dose'][,outcome:=outc]
    
  }else if(outc%in%beforedeath){
    donors_apoecs[,summary(lm(unlist(.SD)~dose+msex+educ+apoe4_dose))$coefficients|>data.table(keep.rownames = 'cov'),
                  by='variant_id2',.SDcols=outc][cov=='dose'][,outcome:=outc]
  }else{
    
    donors_apoecs[,summary(lm(unlist(.SD)~dose+msex+educ+age_death+apoe4_dose))$coefficients|>data.table(keep.rownames = 'cov'),
                  by='variant_id2',.SDcols=outc][cov=='dose'][,outcome:=outc]
  }
  
}))
apoecs_clin2[,zscore:=`t value`]
apoecs_clin2[,pval:=`Pr(>|t|)`]

apoecs_clogit2<-rbindlist(lapply(covs_logis,function(outc){
  message(outc)
  donors_apoecs[,summary(glm(unlist(.SD)~dose+msex+educ+age_death+apoe4_dose,family = binomial))$coefficients|>data.table(keep.rownames = 'cov'),
                by='variant_id2',.SDcols=outc][cov=='dose'][,outcome:=outc]
  
}))
apoecs_clogit2[,zscore:=`z value`]
apoecs_clogit2[,pval:=`Pr(>|z|)`]

apoecs_clin2<-rbind(apoecs_clin2,apoecs_clogit2,fill=TRUE)

apoecs_clin2<-merge(apoecs_clin2,indep_cs,by='variant_id2')
apoecs_clin2[,padj:=p.adjust(pval,n = length(unique(indep_cs$cs_name))),by=c('outcome','variant_id2')]
apoecs_clin2[padj<0.05]


apoecs_clin2[!outcome%in%c('apoe4_dose','apoe2_dose')][padj<0.05]

apoecs_clin2[padj<0.05][!outcome%in%c('apoe4_dose','apoe2_dose')]$cs_id  |>unique()|>length()#3/19 

fwrite(apoecs_clin2,fp(out,'res_lm_apoecs_adpathology_apoe4_correction.csv.gz'))


unique(apoecs_clin2[padj<0.05][!outcome%in%c('apoe4_dose','apoe2_dose')][order(cs_id,cov,padj)][,.(cs_id,outcome,variant_id2,padj,gwas_source,molecular_outcome)],by=c('cs_id','outcome'))
#     cs_id             outcome              variant_id2        padj
#    <char>              <char>                   <char>       <num>
# 1: set_26 age_first_ad_dx_num       chr19:44918487_G_T 0.046608400
# 2: set_27             caa_4gp       chr19:44954310_T_C 0.002833505
# 3: set_27            caa_neo4       chr19:44954310_T_C 0.004877779
# 4: set_27             tdp_st4 chr19:44999110_TAAAA_TAA 0.042966548
# 5: set_31 age_first_ad_dx_num   chr19:45114770_ACCC_AC 0.020474889
# 6: set_57        tangles_sqrt       chr19:44840322_G_A 0.020962988
#                                                                                                                 gwas_source
#                                                                                                                      <char>
# 1:                                          AD_Kunkle_Stage1_2019; AD_Wightman_Full_2021; AD_Wightman_Excluding23andMe_2021
# 2: AD_Kunkle_Stage1_2019; AD_Wightman_Full_2021; AD_Wightman_Excluding23andMe_2021; AD_Wightman_ExcludingUKBand23andME_2021
# 3: AD_Kunkle_Stage1_2019; AD_Wightman_Full_2021; AD_Wightman_Excluding23andMe_2021; AD_Wightman_ExcludingUKBand23andME_2021
# 4: AD_Kunkle_Stage1_2019; AD_Wightman_Full_2021; AD_Wightman_Excluding23andMe_2021; AD_Wightman_ExcludingUKBand23andME_2021
# 5:                                                                                                            AD_Bellenguez
# 6:                           AD_Bellenguez; AD_Kunkle_Stage1_2019; AD_Wightman_Full_2021; AD_Wightman_Excluding23andMe_2021
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     molecular_outcome
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <char>
# 1:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               AC_DeJager_eQTL_ENSG00000130208; Metabrain_Cortex_chr19_41840000_47960000_ENSG00000130208; STARNET_eQTL_Mac_ENSG00000130208; ROSMAP_AC_ENSG00000130208; ROSMAP_DLPFC_ENSG00000130208
# 2: AC_DeJager_eQTL_ENSG00000224916; AC_DeJager_eQTL_ENSG00000234906; BM_10_MSBB_eQTL_ENSG00000104853; Metabrain_Basalganglia_chr19_41840000_47960000_ENSG00000234906; Metabrain_Cerebellum_chr19_41840000_47960000_ENSG00000234906; Metabrain_Cortex_chr19_41840000_47960000_ENSG00000104853; Metabrain_Cortex_chr19_41840000_47960000_ENSG00000234906; Mic_DeJager_eQTL_ENSG00000130208; Mic_DeJager_eQTL_ENSG00000130203; Mic_mega_eQTL_ENSG00000130203; STARNET_eQTL_Mac_ENSG00000130208; STARNET_eQTL_Mac_ENSG00000130203; STARNET_eQTL_Mac_ENSG00000267467; STARNET_eQTL_Mac_ENSG00000224916; STARNET_eQTL_Mac_ENSG00000234906; DeJager_Mic_ENSG00000130203; Kellis_Mic_ENSG00000130203; DeJager_Mic_ENSG00000130208; Kellis_Mic_ENSG00000130208
# 3: AC_DeJager_eQTL_ENSG00000224916; AC_DeJager_eQTL_ENSG00000234906; BM_10_MSBB_eQTL_ENSG00000104853; Metabrain_Basalganglia_chr19_41840000_47960000_ENSG00000234906; Metabrain_Cerebellum_chr19_41840000_47960000_ENSG00000234906; Metabrain_Cortex_chr19_41840000_47960000_ENSG00000104853; Metabrain_Cortex_chr19_41840000_47960000_ENSG00000234906; Mic_DeJager_eQTL_ENSG00000130208; Mic_DeJager_eQTL_ENSG00000130203; Mic_mega_eQTL_ENSG00000130203; STARNET_eQTL_Mac_ENSG00000130208; STARNET_eQTL_Mac_ENSG00000130203; STARNET_eQTL_Mac_ENSG00000267467; STARNET_eQTL_Mac_ENSG00000224916; STARNET_eQTL_Mac_ENSG00000234906; DeJager_Mic_ENSG00000130203; Kellis_Mic_ENSG00000130203; DeJager_Mic_ENSG00000130208; Kellis_Mic_ENSG00000130208
# 4: AC_DeJager_eQTL_ENSG00000224916; AC_DeJager_eQTL_ENSG00000234906; BM_10_MSBB_eQTL_ENSG00000104853; Metabrain_Basalganglia_chr19_41840000_47960000_ENSG00000234906; Metabrain_Cerebellum_chr19_41840000_47960000_ENSG00000234906; Metabrain_Cortex_chr19_41840000_47960000_ENSG00000104853; Metabrain_Cortex_chr19_41840000_47960000_ENSG00000234906; Mic_DeJager_eQTL_ENSG00000130208; Mic_DeJager_eQTL_ENSG00000130203; Mic_mega_eQTL_ENSG00000130203; STARNET_eQTL_Mac_ENSG00000130208; STARNET_eQTL_Mac_ENSG00000130203; STARNET_eQTL_Mac_ENSG00000267467; STARNET_eQTL_Mac_ENSG00000224916; STARNET_eQTL_Mac_ENSG00000234906; DeJager_Mic_ENSG00000130203; Kellis_Mic_ENSG00000130203; DeJager_Mic_ENSG00000130208; Kellis_Mic_ENSG00000130208
# 5:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               AC_DeJager_eQTL_ENSG00000189114; DLPFC_DeJager_eQTL_ENSG00000189114; Metabrain_Cortex_chr19_41840000_47960000_ENSG00000189114; PCC_DeJager_eQTL_ENSG00000189114; ROSMAP_AC_ENSG00000189114; ROSMAP_DLPFC_ENSG00000189114; ROSMAP_PCC_ENSG00000189114
# 6:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Ast_10_Kellis_eQTL_ENSG00000160007
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        Ast_10_Kellis_eQTL_ENSG00000160007
#set26 is APOE2?

apoecs_clin2[cs_id=='set_26'&pos(variant_id2)==44908822]#no

donors_apoecs[variant_id2=='chr19:44918487_G_T']$apoe2_r[1]#-0.21
donors_apoecs[variant_id2=='chr19:44918487_G_T']$apoe4_r[1]#-0.23

#CCL: We can validate 3 of the 19 independant variants sets independently explaining one neuropathology after adjusted for APOE4 and 2 dose (along with Sex, education and age covariates),
#and 4/19 sets once we do not adjust for APOE4. These results is for FDR<0.05, using BH correction with n=19 tests

# breakdown of the neuropathology associated sets:
#   set27 is the APOE Mic eQTL that have previously been found associated with cerebral amyloid angiopathy in Fujita et al, and we are able to recapitulate this association. The 3 others sets are set_26, 31 and 57, 
# respectively eQTL of APOC1 (in several brain regions and STARNET Mac), BLOC1S3( in several brain regions ), both explaining age at AD onset, and ARHGAP35 (in Ast.10 specifically), explaining Tau tangle level.
