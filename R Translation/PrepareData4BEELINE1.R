########################################################################
#         CICT for Single Cell RNA Seq
#      This code is an R translation of Causal Inference Using Composition of Transaction (CICT) 
#      For original paper and concepts, please refer to the paper :
#              Revisiting Causality Inference in Memory-less Transition Networks  at                      https://arxiv.org/abs/1608.02658
#       
#        Abbas Shojaee, 2016-2022
#        All rights reserved. Please do not use, reuse or distribute in any form
#        The code is provided only for use at Huang Lab for Genomic Analysis at NYU https://huanglab.rbind.io/
########################################################################

options(error=NULL)
  #Libraries -----
setwd("C:/E/HDDs/Work/NYU Research/0 Plant biology/R code")
#library(officer)

#library(mlbench)
library(caret)
#library(Spectrum)
library(Rtsne)
library(umap)
library(devtools)
library(PerformanceAnalytics)
library(sn)

library(tidyr)
library(infotheo)
library(mpmi)
library(stringr)
library(RCy3) #https://github.com/cytoscape/cytoscape-automation/wiki
library(tidyr)
library(infotheo)
library(mpmi)
library(stringr)
library(igraph)
library(DREAM4)
library(data.table)
library(plyr)

library(tidyverse)
library(dplyr)
select<-dplyr::select
mutate<-dplyr::mutate
arrange<-dplyr::arrange
filter<-dplyr::filter
group_by<-dplyr::group_by
require(zeallot)

#IDEAS:
#1- directed edge between each gene and all genes in previous time frame
#2- Use a second round RF on CICT predictions, keep most important influencers as causes
#3- replace mutual information with a better including a conditional on all others form
#4- MI (g1,gx.LAGk ) | MI(g1,gn.LAG1)  #MI(g1,gx ge and particular gx using all consequent time frames, conditioned on g1 and all gj in the previous time frame



  
  #NEW VERSION ========================================== 
  prepDataForBeeline.v1<-function(source.path, dest.path,datasets, topvaryingGN=1000,
                                  randomGN=0, includeTFs = T,ftype = 'tsv',genesInCols=F)
  {
  
    readfile = function(file,ftype,header=T,genesInCols = F){
      if(ftype == 'tsv') .m= read_tsv(thefile,col_names = header) #
      if(ftype == 'csv') .m=read.csv(thefile,header = header)
      #if(genesInCols == 'auto') if()
      if(genesInCols) .m = t(.m) 
      return(.m)
      
    }
    
    for(i in 1:nrow(datasets)){
      drow = datasets[i,]
      dbpath = paste0(dest.path,"/",drow$d.name)
      if(!dir.exists(dbpath)) dir.create(dbpath)
      #Realnames
      thefile = paste0(source.path,"/",drow$d.realnames)
      if(file.exists(thefile)) d.realnames = readfile(thefile,ftype,header = T) else {print(sprintf("Realnames file %s not found: ",drow$d.realnames ))} 
      
      thefile = paste0(source.path,"/",drow$d.exp)
      if(file.exists(thefile)) d.exp = readfile(thefile,ftype,header = T, genesInCols=genesInCols) else {print(sprintf("Expression file %s not found: ",drow$d.name )) }
      try({d.exp = d.exp %>% column_to_rownames('gene')})
      if(exists('d.realnames') ) rownames(d.exp) = mapvalues(rownames(d.exp),d.realnames$`#ID`,d.realnames$Name, warn_missing = F)
         #TOOD replace all gene names with the real names)
      
      thefile = paste0(source.path,"/",drow$d.gs)
      if(file.exists(thefile)) d.gs = readfile(thefile,ftype,header = T) else {print(sprintf("Ground truth file %s not found: ",drow$d.gs ))}
      if(ncol(d.gs)==3) 
        colnames(d.gs) <-c('Gene1','Gene2','Type') else {
          colnames(d.gs) <-c('Gene1','Gene2')
          d.gs$edgetruth = '1' #Add a third column
        }
      if(exists('d.realnames'))  d.gs = d.gs %>% mutate(Gene1= mapvalues(Gene1,d.realnames$`#ID`,d.realnames$Name, warn_missing = F),
                                                        Gene2= mapvalues(Gene2,d.realnames$`#ID`,d.realnames$Name, warn_missing = F))
         
      
      thefile = paste0(source.path,"/",drow$d.tf)
      if(file.exists(thefile)) d.tf = readfile(thefile,ftype) else {
        print(sprintf("Transcription file %s not found: ",drow$d.tf ));
        d.tf = d.gs %>% filter(edgetruth ==1) %>% select(src) %>% unique()
      }; colnames(d.tf)[1]<- 'tf'
      if(exists('d.realnames')) d.tf = d.tf %>% mutate(tf= mapvalues(tf,d.realnames$`#ID`,d.realnames$Name, warn_missing = F))
      
      
      #Convert numeric column labels = cells  to C1 C2 etc
      clnames = colnames(d.exp)
      
      #if(! clnames[-1] %>% as.numeric() %>% is.na()  %>% all())
      colnames(d.exp) <- ifelse(sapply(clnames,is.character,simplify = T),clnames,paste0('C',clnames))
      
      if(row.names(d.exp) %>% as.numeric() %>% is.na()  %>% all() &
         ! any( str_detect(colnames(d.exp),'(G|g)ene')) )
        all.dt = as.data.frame(d.exp) %>% rownames_to_column('gene') %>% select(gene,everything()) #rename(gene = X) # %>% select(any_of('vID',hseq.cellset1$V1))
      else {
        if(colnames(d.exp)[1]=='X' & class(d.exp[,1])=="character") colnames(d.exp)[1]<- 'gene'
        all.dt = as.data.frame(d.exp) %>% select(gene,everything())
      }
      rownames(all.dt)<- all.dt$gene 
      #hsec1 =  cbind(X = d.exp.sbst1$X, exp = rowMeans(d.exp.sbst1[,2:ncol(d.exp.sbst1)], dims = 1))
      
      thefile = paste0(source.path,"/",drow$d.order)
      if(file.exists(thefile)) geneorder = readfile(thefile,ftype) else {
        library(Seurat)
        library(cowplot)
        #pbmc.data <- Read10X(data.dir = "data/pbmc3k_filtered_gene_bc_matrices/hg19/")
        
        all.dt.s <- CreateSeuratObject(counts = all.dt, min.cells = 3, min.features  = 1, project = "dream5", assay = "RNA")
        all.dt.n <- NormalizeData(object = all.dt.s, normalization.method = "LogNormalize", scale.factor = 10000)
        all.dt.f <- FindVariableFeatures(object = all.dt.n, mean.function = ExpMean, dispersion.function = LogVMR, 
                                         x.low.cutoff = 0.0001, x.high.cutoff = 20,  nfeatures = nrow(all.dt)) #y.cutoff = 0.5,
        
        geneorder=  HVFInfo(object = all.dt.f) %>%  mutate(VGAMpValue = 0.001 , Variance = variance)  %>% 
          rownames_to_column('gene') %>%  select(gene,Variance,VGAMpValue)
        
      }
      
      if(colnames(geneorder)[1]=='X' & class(geneorder[,1])=="character") colnames(geneorder)[1]<- 'gene'
      finalgeneorder = geneorder %>% arrange(desc(Variance)) %>%
        filter(row_number() <=topvaryingGN | gene %in% d.tf$tf) %>% unique()
      
      try({write.csv(d.tf,paste0(dest.path,"/",drow$d.name,"/tf.csv"),row.names = F,quote = F)})
      write.csv(finalgeneorder,paste0(dest.path,"/",drow$d.name,"/GeneOrdering.csv"),row.names = F,quote = F)
        
      randomGenes = c()
      if(randomGN>0) randomGenes = geneorder %>% 
                          anti_join(finalgeneorder,by = 'gene') %>% 
                          sample_n(size = randomGN) %>% pull(gene)
      
      write.csv(d.gs,paste0(dest.path,"/",drow$d.name,"/refNetwork.csv"),row.names = F,quote = F)
      
      
      thefile = paste0(source.path,"/",drow$d.pseudo)
      if(file.exists(thefile)) file.copy(thefile,paste0(dest.path,"/",drow$d.name,"/PseudoTime.csv")) else{
        d.pseudo = data.frame(cellName = colnames(all.dt)[-1], PseudoTime = 0.1) #Create a pseudotime for cells all the same time
        write.csv(d.pseudo,paste0(dest.path,"/",drow$d.name,"/PseudoTime.csv"),row.names = F,quote = F)
      }
      
      finalExpData = all.dt %>% filter(gene %in% c(finalgeneorder$gene,randomGenes))
      write.csv(finalExpData,paste0(dest.path,"/",drow$d.name,"/ExpressionData.csv"),row.names = F,quote = F)
    }
    
  }
  

  #DREAM Datasets
  {
    dream.datasets=list(list('dream5_1', 'net1_expression_data.tsv','DREAM5_NetworkInference_GoldStandard_Network1 - in silico.tsv','net1_transcription_factors.tsv','net1_gene_ids.tsv','',''),
                        list('dream5_3','net3_expression_data.tsv','DREAM5_NetworkInference_GoldStandard_Network3 - E. coli.tsv','net3_transcription_factors.tsv','net3_gene_ids.tsv','',''),
                        list('dream5_4','net4_expression_data.tsv','DREAM5_NetworkInference_GoldStandard_Network4 - S. cerevisiae.tsv','net4_transcription_factors.tsv','net4_gene_ids.tsv','',''))
    dream.datasets = rbindlist(dream.datasets) ; 
    colnames(dream.datasets) <-c('d.name' ,'d.exp','d.gs','d.tf','d.realnames','d.order','d.pseudo')
  
    
    source.path = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/CICT Beeline config/data/Dream 5 Train and Ground Truth data"
    dest.path = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/CICT Beeline config/data/L1"
  
    prepDataForBeeline.v1(source.path, dest.path,dream.datasets,topvaryingGN=1000,randomGN=0,
                          includeTFs = T,ftype = 'tsv',genesInCols=F)
}    
    

  #BEELINE Datasets
  {
    beeline.datasets=list(list('hESC', 'BEELINE-data/inputs/scRNA-Seq/hESC/ExpressionData.csv',
                               'BEELINE-Networks/Networks/human/hESC-ChIP-seq-network.csv',
                               'BEELINE-Networks/human-tfs.csv', #Transcription factors
                               'NULL', #RealNames
                               'BEELINE-data/inputs/scRNA-Seq/hESC/GeneOrdering.csv',
                               'BEELINE-data/inputs/scRNA-Seq/hESC/PseudoTime.csv'),
                          
                          list('hHEP', 'BEELINE-data/inputs/scRNA-Seq/hHEP/ExpressionData.csv',
                               'BEELINE-Networks/Networks/human/HepG2-ChIP-seq-network.csv',
                               'BEELINE-Networks/human-tfs.csv', 'NULL', #RealNames
                               'BEELINE-data/inputs/scRNA-Seq/hHEP/GeneOrdering.csv',
                               'BEELINE-data/inputs/scRNA-Seq/hHEP/PseudoTime.csv'),
                          
                          list('mDC', 'BEELINE-data/inputs/scRNA-Seq/mDC/ExpressionData.csv',
                               'BEELINE-Networks/Networks/mouse/mDC-ChIP-seq-network.csv',
                               'BEELINE-Networks/mouse-tfs.csv', 'NULL', #RealNames
                               'BEELINE-data/inputs/scRNA-Seq/mDC/GeneOrdering.csv',
                               'BEELINE-data/inputs/scRNA-Seq/mDC/PseudoTime.csv'),
                          
                          list('mESC', 'BEELINE-data/inputs/scRNA-Seq/mESC/ExpressionData.csv',
                               'BEELINE-Networks/Networks/mouse/mESC-ChIP-seq-network.csv',
                               'BEELINE-Networks/mouse-tfs.csv', 'NULL', #RealNames
                               'BEELINE-data/inputs/scRNA-Seq/mESC/GeneOrdering.csv',
                               'BEELINE-data/inputs/scRNA-Seq/mESC/PseudoTime.csv'),
                          
                          list('mHSC-E', 'BEELINE-data/inputs/scRNA-Seq/mHSC/ExpressionData.csv',
                               'BEELINE-Networks/Networks/mouse/mHSC-ChIP-seq-network.csv',
                               'BEELINE-Networks/mouse-tfs.csv', 'NULL', #RealNames
                               'BEELINE-data/inputs/scRNA-Seq/mHSC/GeneOrdering.csv',
                               'BEELINE-data/inputs/scRNA-Seq/mHSC/PseudoTime.csv'))
    
    beeline.datasets = rbindlist(beeline.datasets) ; 
    colnames(beeline.datasets) <-c('d.name' ,'d.exp','d.gs','d.tf','d.realnames','d.order','d.pseudo')
    
    source.path = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/data"
    dest.path = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/CICT Beeline config/data/L1"
   
    prepDataForBeeline.v1(source.path, dest.path,beeline.datasets,topvaryingGN=1000,randomGN=0,
                          includeTFs = T,ftype = 'csv',genesInCols=F)
  }
  
  #============= OLD editions ====================
  # Prepare beeline data for kuzmanovski pipeline
 if(F){ 
  
    studyDataset = 'hesc'
    d.exp.baseurl = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/Comparison/vkuzmanovski-rn-approach/datasets/For beeline l1/"
    d.exp.baseurl = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/CICT Beeline config/data/"
    #"/BEELINE-data/inputs/scRNA-Seq/hESC/ExpressionData.csv"
    #"/BEELINE-data/inputs/scRNA-Seq/mHSC-E/GeneOrdering.csv"
    #Standardize all db names #"/BEELINE-Networks/Networks/human/hESC-ChIP-seq-network.csv"
    prepBeelineForKuzma<-function(d.exp.baseurl,studyDataset, topvaryingGN=1000,randomGN=1000, includeTFs = T, groundTruthSuffix ="-ChIP-seq-network.csv")
    {
      dbpath = paste0(d.exp.baseurl,"/",studyDataset)
      d.exp = read.csv(paste0(dbpath,"/ExpressionData.csv"),header = T)
      
      all.dt = d.exp %>% rename(gene = X) # %>% select(any_of('vID',hseq.cellset1$V1))
      rownames(all.dt)<- all.dt$gene 
      #hsec1 =  cbind(X = d.exp.sbst1$X, exp = rowMeans(d.exp.sbst1[,2:ncol(d.exp.sbst1)], dims = 1))
      
      geneorder = read.csv(paste0(dbpath,"/GeneOrdering.csv"),header = T)%>% rename(gene = X)
      geneorder.n = geneorder %>%  arrange(desc(Variance)) %>% filter(VGAMpValue<0.01) %>% slice_head(n=topvaryingGN) #
      
      d.gtchip = read.csv(paste0(dbpath,"/",studyDataset,groundTruthSuffix),header = T)
      d.gtchip$edgetyptruth = 'chipseq'
      colnames(d.gtchip) <-c('src','trgt','edgetyptruth') 
      
      tbl.goldStandard = d.gtchip
      
      #TFs are all the  genes controling others in chipseq
      tf.all = unique(tbl.goldStandard$src) # when Chip_seq , gives number of TFs
      #all TFs whose variance had p-value at most 0.01
      tf.set1= geneorder %>% anti_join(geneorder.n,by="gene") %>%  filter(gene %in% tf.all ) %>% arrange(desc(Variance))  #&  VGAMpValue>=0.01
      
      
      semifinalGeneSet = geneorder %>% filter(gene %in% tf.set1$gene) %>% rbind(geneorder.n) %>% unique()
      
      randomgenes = sample_n(setDT(geneorder)[! gene %in% semifinalGeneSet$gene ],size=randomGN)
      
      finalGeneset = rbind(semifinalGeneSet,randomgenes)
      
      all.dt = all.dt %>% filter(gene %in% finalGeneset$gene) 
      all.dt
    }
    
    
    studyDataset = 'hesc'
    hescTop1000 = prepBeelineForKuzma(d.exp.baseurl,studyDataset,topvaryingGN=1000,randomGN=0, includeTFs = T)
    write.csv(hescTop1000,file = paste0(d.exp.baseurl,studyDataset,"/",studyDataset,'_for_Kuzma.csv'),row.names = F)
    
    studyDataset = 'hHep'
    hHepTop1000 = prepBeelineForKuzma(d.exp.baseurl,studyDataset,topvaryingGN=1000,randomGN=0, includeTFs = T)
    write.csv(hHepTop1000,file = paste0(d.exp.baseurl,studyDataset,"/",studyDataset,'_for_Kuzma.csv'),row.names = F)
    
    studyDataset = 'mESC'
    mESCTop1000 = prepBeelineForKuzma(d.exp.baseurl,studyDataset,topvaryingGN=1000,randomGN=0, includeTFs = T)
    write.csv(mESCTop1000,file = paste0(d.exp.baseurl,studyDataset,"/",studyDataset,'_for_Kuzma.csv'),row.names = F)
    
    studyDataset = 'mDC'
    mDCTop1000 = prepBeelineForKuzma(d.exp.baseurl,studyDataset,topvaryingGN=1000,randomGN=0, includeTFs = T)
    write.csv(mDCTop1000,file = paste0(d.exp.baseurl,studyDataset,"/",studyDataset,'_for_Kuzma.csv'),row.names = F)
    
    # length(unique(c(tbl.goldStandard$src,tbl.goldStandard$trgt)))
    # 
    # Nobservations=ncol(all.dt)-1
    # ngens = nrow(all.dt)
    # 
    # a=sapply(similarityMatrices[[i]], function(x) sum(is.na(x)));a[a>0]
    # 
    # d.psuedo = fread(paste0(dbpath,"/PseudoTime.csv"),header = T);#hist(d.psuedo$PseudoTime)
    # 
    # recipes::discretize(d.psuedo$PseudoTime,cuts=5)
    # d.psuedo = d.psuedo %>%  dplyr::mutate(cellset = infotheo::discretize(d.psuedo$PseudoTime, "equalfreq", 5)$X);
    # hist(d.psuedo$cellset);table(d.psuedo$cellset)
    #hseq.cellset1 = d.psuedo[PseudoTime<=0.02,]
    
    
    studyDataset = 'dream5_4'
    prepDataForBeeline<-function(d.exp.baseurl,studyDataset, topvaryingGN=1000,randomGN=1000, includeTFs = T)
    {
      
      dbpath = paste0(d.exp.baseurl,"/",studyDataset)
      d.exp = read.csv(paste0(dbpath,"/DREAM_DS_net4_ss1.csv"),header = T)
      
      all.dt = d.exp %>% rownames_to_column('gene') #rename(gene = X) # %>% select(any_of('vID',hseq.cellset1$V1))
      rownames(all.dt)<- all.dt$gene 
      #hsec1 =  cbind(X = d.exp.sbst1$X, exp = rowMeans(d.exp.sbst1[,2:ncol(d.exp.sbst1)], dims = 1))
      
      library(Seurat)
      library(cowplot)
      #pbmc.data <- Read10X(data.dir = "data/pbmc3k_filtered_gene_bc_matrices/hg19/")
      
      all.dt.s <- CreateSeuratObject(counts = all.dt, min.cells = 3, min.features  = 1, project = "dream5", assay = "RNA")
      all.dt.n <- NormalizeData(object = all.dt.s, normalization.method = "LogNormalize", scale.factor = 10000)
      all.dt.f <- FindVariableFeatures(object = all.dt.n, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125,
                                       x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
      
      geneorder=  HVFInfo(object = all.dt.f) %>%  mutate(VGAMpValue = 0.001 , Variance = variance)  %>% 
        rownames_to_column('gene') %>%  select(gene,Variance,VGAMpValue)
      
      tbl.gs = read.csv(paste0(dbpath,"/DREAM_GS_net4_ss1.csv"),header = T)
      colnames(tbl.gs) <-c('src','trgt','edgetruth') 
      
      tfs = tbl.gs %>% filter(edgetruth ==1) %>% pull(src) %>% unique()
      
      
      finalgeneorder = geneorder %>% arrange(desc(Variance)) %>%
        filter(row_number() <=topvaryingGN | gene %in% tfs) %>% unique()
      write.csv(finalgeneorder,paste0(dbpath,"/GeneOrdering.csv"),row.names = F)
      
      finalExpData = all.dt %>% filter(gene %in% finalgeneorder$gene)
      write.csv(finalExpData,paste0(dbpath,"/ExpressionData.csv"),row.names = F)
      
      
    }
  
  }