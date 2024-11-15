#!/usr/bin/Rscript
set.seed(12)   
args = commandArgs(trailingOnly=TRUE)
library(stringr)
library(RColorBrewer)
library(foreach)
#library(ggsci)
library(ggplot2)
library(gridExtra)
#library(viridis)
#library(gplots)
library(data.table)
library(dplyr)
library(beepr)
library(parallel)
Sjscorethreshold = 10
scorediffthreshold= 10
c = 0.001

metadatafile = args[1]
mutmetadata  =data.frame(fread(metadatafile))
#if(sum(mutmetadata$relevance == "Ambiguous")>0){
  
  
  PANCANCERNFE2L2 =readRDS("sup_files/mc3.v0.2.8.PUBLIC_NRF2.maf.RDS")
  npvalue = 1
  resultadosnumericPANCANCER =   data.frame(fread("output/Scores_Ambiguous.tsv"))
  
  
  resultadosnumericPANCANCER = unique(resultadosnumericPANCANCER)
  resultadosnumericPANCANCER= resultadosnumericPANCANCER[!duplicated(resultadosnumericPANCANCER$TCGA_ID),]
  colnames(resultadosnumericPANCANCER)[3] <- "Mutation labels"
  resultadosnumericPANCANCERforplot = resultadosnumericPANCANCER[order(resultadosnumericPANCANCER$score, resultadosnumericPANCANCER$`Mutation labels`, decreasing = T),]
  resultadosnumericPANCANCERforplot$`Mutation labels` = as.character(resultadosnumericPANCANCERforplot$`Mutation labels`)
  resultadosnumericPANCANCERforplot$`Mutation labels`[resultadosnumericPANCANCERforplot$`Mutation labels`=="NFE2L2_Exon2_Mut"] <- "NRF2_Exon2_Mut"
  resultadosnumericPANCANCERforplot$`Mutation labels`[resultadosnumericPANCANCERforplot$`Mutation labels`=="NFE2L2_Other_Exon_Mut"] <- "NRF2_Other_Exon_Mut"
  resultadosnumericPANCANCERforplot$`Mutation labels`[resultadosnumericPANCANCERforplot$`Mutation labels`=="Exon2_Skipping"] <- "NRF2_Exon2_Skipping"
  old = c("None" ,"NRF2_Exon2_Mut",  "NRF2_Other_Exon_Mut" ,"KEAP1" , "CUL3"  ,"CNA" ,"NRF2_Exon2_Skipping", "Multihit")
  old = c("None" ,"NRF2_Exon2_Hotspot_GOF",  "NRF2_Exon2_Mut" , "NRF2_Other_Exon_Mut" ,"KEAP1" ,"KEAP1 OL" ,"KEAP1 Ambiguous",  "CUL3"  ,"CNA" ,"NRF2_Exon2_Skipping", "Multihit")
  new = c("None" ,"NRF2 Exon 2 Mut",  "NRF2 Other Exon Mut" ,"KEAP1 Mut" ,"KEAP1_OL" ,"KEAP1_Ambiguous", "CUL3 Mut"  ,"NRF2 CNA" ,"NRF2 Exon 2 Skipping", "Multihit")
  new = c("None" ,"NRF2 Exon 2 Hotspot GOF", "NRF2 Exon2 Mut" , "NRF2 Other Exon Mut" ,"KEAP1 Mut" , "KEAP1_OL" ,"KEAP1_Ambiguous", "CUL3 Mut"  ,"NRF2 CNA" ,"NRF2 Exon 2 Skipping", "Multihit")
  new = c("None" ,"NRF2 E2 HGOF, "NRF2 E2" , "NRF2 oE" ,"KEAP1 Mut" , "KEAP1 OLO" ,"KEAP1 Amb", "CUL3"  ,"NRF2 CNA" ,"NRF2 E2 Skip", "Multihit")
  
  
 ######







######  
  
  
  
  
  
  
  for(cont in 1:length(old)){
    resultadosnumericPANCANCERforplot$`Mutation labels`[resultadosnumericPANCANCERforplot$`Mutation labels` == old[cont]] = new[cont]
    
    
  }
  #resultadosnumericPANCANCERforplot$`Mutation labels` =  factor(resultadosnumericPANCANCERforplot$`Mutation labels`, levels = c("None" ,"NRF2_Exon2_Mut",  "NRF2_Other_Exon_Mut" ,"KEAP1" , "CUL3"  ,"CNA" ,"NRF2_Exon2_Skipping", "Multihit"))
  
  resultadosnumericPANCANCERforplot$`Mutation labels` =  factor(resultadosnumericPANCANCERforplot$`Mutation labels`, levels = new[new %in% resultadosnumericPANCANCERforplot$`Mutation labels`])
  
  
  resultadosnumericPANCANCERforplot$is_mut_predicted = resultadosnumericPANCANCERforplot$score >scorediffthreshold
  resultadosnumericPANCANCERforplot$is_mut_observed = resultadosnumericPANCANCERforplot$`Mutation labels` != "None"
  
  resultadosnumericPANCANCERforplot$is_mut_predicted = factor(resultadosnumericPANCANCERforplot$is_mut_predicted, levels = c(FALSE, TRUE))
  
  
  resultadosnumericPANCANCERforplot$is_mut_observed = factor(resultadosnumericPANCANCERforplot$is_mut_observed, levels = c(FALSE, TRUE))
  
  
  #colnames(resultstest) <- c("Observed Not Mutated", "Observed Mutated")
  resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous` = as.character(resultadosnumericPANCANCERforplot$`Mutation labels`)
  
  #resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous`[rowSums(resultadosnumericPANCANCERforplot[,c(7,9:12)])>1]<- "True Multihit"
  #
  #resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous`[resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous` == "Multihit"] <- colnames(resultadosnumericPANCANCERforplot[resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous` == "Multihit",c(7,9:12)] )[which(resultadosnumericPANCANCERforplot[resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous` == "Multihit",c(7,9:12)] == TRUE,arr.ind = T)[,2]]
  #resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous`[resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous` == "True Multihit"] <- "Multihit"
  
  levelsambiguous = c("NRF2 Exon2 Mut","NRF2 Exon2 Hotspot GOF ", "NRF2 Other Exon Mut" ,"CUL3 Mut", "NRF2 CNA" ,"NRF2 Exon 2 Skipping", "Multihit")
  levelsambiguous = levels(resultadosnumericPANCANCERforplot$`Mutation labels`) [levels(resultadosnumericPANCANCERforplot$`Mutation labels`) %in% unique(resultadosnumericPANCANCERforplot$`Mutation labels`)]
  
  resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous` = factor(resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous`, levels = levelsambiguous)
  resultadosnumericPANCANCERforplot$`Mutation labels` <- resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous`
  resultstest = table(resultadosnumericPANCANCERforplot$is_mut_predicted,resultadosnumericPANCANCERforplot$`Mutation labels`)
  rownames(resultstest) <- c("Predicted Not Mutated", "Predicted Mutated")
  
  ggbx =ggplot(resultadosnumericPANCANCERforplot, aes(x = ((1:length(`Mutation labels`))-1) , y =  score, color = `Mutation labels`, fill = `Mutation labels`)) + geom_bar(stat="identity",position="dodge") +xlab("Sample") + ylab("Score") +ggtitle(paste0("TCGA dataset mutation prediction score SJ classifier"))+
    annotation_custom( tableGrob(resultstest), xmin=length(resultadosnumericPANCANCERforplot$score)/2-1, xmax=length(resultadosnumericPANCANCERforplot$score)[1]/2, ymin=max(resultadosnumericPANCANCERforplot$score)/2-1, ymax=max(resultadosnumericPANCANCERforplot$score/2))
  #ggbx = ggbx +scale_fill_nejm()
  
  
  myColors <- brewer.pal(7,"Set1")
  myColors <- brewer.pal(7,"Set1")[c(2,4,5,6,7)]
  myColors <- brewer.pal(7,"Set1")[c(1:7)]
  myColors <-myColors[-3] 
  myColors[1] <-"hotpink"
  #myColors = c("Darkgray", myColors)
  myColors = c(myColors[1:2], "darkseagreen",myColors[3:length(myColors)])
  ggbx = ggbx + scale_fill_manual(values = myColors, drop = FALSE) + scale_color_manual(values = myColors, drop = FALSE)+ theme_minimal()  + geom_hline(yintercept=10, linetype="dashed", color = "red")+
    theme(axis.text.x=element_blank(), #remove x axis labels
          axis.ticks.x=element_blank()  #remove y axis ticks
    )
  
  #ggbx = ggbx + scale_fill_nejm(drop = FALSE) + scale_color_nejm(drop = FALSE)
  ggsave(
    plot = ggbx,
    file = "output/figures/Barplot_Scores_Ambiguous.pdf",
    height = 10,
    width = 50 ,
    units = "cm",
    device =
      "pdf"
  )
  
  
  pdf("output/figures/Boxplot_Scores_Ambiguous.pdf", width = 30)
  
  p = ggplot(resultadosnumericPANCANCERforplot, aes(x = `Mutation labels` , y =  score, fill = `Mutation labels`)) + geom_boxplot( color = "black")+ scale_fill_manual(values = myColors, drop = FALSE) + scale_color_manual(values = myColors, drop = FALSE)+ theme_minimal()  + geom_hline(yintercept=10, linetype="dashed", color = "red")
  print(p)
  dev.off()
  

  
  #stat="identity",position="dodge", width=0.02) +xlab("Sample") + ylab("Score") +ggtitle(paste0("TCGA dataset mutation prediction score SJ classifier"))+
  #  annotation_custom( tableGrob(resultstest), xPANCANCERNFE2L2min=length(resultadosnumericPANCANCERforplot$score)/2-1, xmax=length(resultadosnumericPANCANCERforplot$score)[1]/2, ymin=max(resultadosnumericPANCANCERforplot$score)/2-1, ymax=max(resultadosnumericPANCANCERforplot$score/2))
  #
  
  
  
  
  
  #ggbx = ggbx +scale_fill_nejm()
  
  resultadosnumericPANCANCERforplot$topcancers = resultadosnumericPANCANCERforplot$cancer
  resultadosnumericPANCANCERforplotmut = resultadosnumericPANCANCERforplot[resultadosnumericPANCANCERforplot$`Mutation labels` != "None",]
  topcancers  =names(table(resultadosnumericPANCANCERforplotmut$topcancers)[order(table(resultadosnumericPANCANCERforplotmut$topcancers),decreasing = T)[1:10]])
  resultadosnumericPANCANCERforplotmut$topcancers[!resultadosnumericPANCANCERforplotmut$topcancers %in% topcancers] <- "Others"
  #ggbx = ggbx +scale_fill_nejm()
  resultadosnumericPANCANCERforplotmut$topcancers =factor(resultadosnumericPANCANCERforplotmut$topcancers ,levels = c( topcancers, "Others"))
  
  colnames(resultadosnumericPANCANCERforplotmut)[colnames(resultadosnumericPANCANCERforplotmut) == "Cancer_Type"] <- "Cancer Type"
  
  colnames(resultadosnumericPANCANCERforplotmut)[colnames(resultadosnumericPANCANCERforplotmut) == "topcancers"] <- "Top Cancers"
  ggbx =ggplot(resultadosnumericPANCANCERforplotmut, aes(x = ((1:length(`Mutation labels`))-1) , y =  score, color = `Top Cancers` , fill = `Top Cancers`)) + geom_bar(stat="identity",position="dodge", width=0.6) +xlab("Sample") + ylab("Score") +ggtitle(paste0("TCGA dataset mutation prediction score SJ classifier"))+
    annotation_custom( tableGrob(resultstest), xmin=length(resultadosnumericPANCANCERforplot$score)/2-1, xmax=length(resultadosnumericPANCANCERforplot$score)[1]/2, ymin=max(resultadosnumericPANCANCERforplot$score)/2-1, ymax=max(resultadosnumericPANCANCERforplot$score/2))
  
  
  myColors <-c("#ebac23","#b80058","#008cf9","#006e00","#878500","#d163e6","#b24502","#ff9287","#5954d6","#00c6f8","#bdbdbd")
  ggbx = ggbx + scale_fill_manual(values = myColors, drop = FALSE) + scale_color_manual(values = myColors, drop = FALSE)+ theme_minimal()  + geom_hline(yintercept=10, linetype="dashed", color = "red")+
    theme(axis.text.x=element_blank(), #remove x axis labels
          axis.ticks.x=element_blank(),  #remove y axis ticks
          plot.title=element_text(size=20),
          axis.title.x = element_text(size = 16),     
          axis.title.y = element_text(size = 16),  
          legend.title = element_text(size=14), #change legend title font size
          legend.text = element_text(size=12)) #change legend text font size
  
  
  
  #ggbx = ggbx + scale_fill_nejm(drop = FALSE) + scale_color_nejm(drop = FALSE)
  pdf("output/figures/Barplot_Scores_TopCancers_Mutated_Ambiguous.pdf", width = 20)
  print(ggbx)
  dev.off()
  
  
  
  ######################################
  
  cancers= unique(resultadosnumericPANCANCERforplotmut$cancer)
  cancers = cancers[order(cancers)]
  pdf("output/figures/Barplot_Scores_per_Cancer_Mutated_Ambiguous.pdf", width = 20)
  
  resultadosnumericPANCANCERforplotmut
  for(cont in 1:length(cancers)){
    resultadosnumericPANCANCERforplotcancer = resultadosnumericPANCANCERforplotmut[resultadosnumericPANCANCERforplotmut$cancer ==cancers[cont], ]
    resultstest = table(resultadosnumericPANCANCERforplotcancer$is_mut_predicted,resultadosnumericPANCANCERforplotcancer$`Mutation labels`)
    rownames(resultstest) <- c("Predicted Not Mutated", "Predicted Mutated")
    
    ggbx =ggplot(resultadosnumericPANCANCERforplotcancer, aes(x = ((1:length(`Mutation labels`))-1) , y =  score, color = `Mutation labels`, fill = `Mutation labels`)) + geom_bar(stat="identity",position="dodge") +xlab("Sample") + ylab("Score") +ggtitle(paste0("TCGA dataset mutation prediction score SJ classifier: ", cancers[cont]))+
      annotation_custom( tableGrob(resultstest), xmin=length(resultadosnumericPANCANCERforplotcancer$score)/2-1, xmax=length(resultadosnumericPANCANCERforplotcancer$score)[1]/2, ymin=max(resultadosnumericPANCANCERforplotcancer$score)/2-1, ymax=max(resultadosnumericPANCANCERforplotcancer$score/2))
    #ggbx = ggbx +scale_fill_nejm()
    
    
    myColors <- brewer.pal(7,"Set1")
    myColors <- brewer.pal(7,"Set1")[c(2,4,5,6,7)]
    myColors <- brewer.pal(7,"Set1")
    myColors <-myColors[-3]   
    myColors[1] <-"hotpink"
    #myColors = c("Darkgray", myColors)
    myColors = c(myColors[1:2], "darkseagreen",myColors[3:length(myColors)])

    ggbx = ggbx + scale_fill_manual(values = myColors, drop = FALSE) + scale_color_manual(values = myColors, drop = FALSE)+ theme_minimal()  + geom_hline(yintercept=10, linetype="dashed", color = "red")+
      theme(axis.text.x=element_blank(), #remove x axis labels
            axis.ticks.x=element_blank()  #remove y axis ticks
      )
    
    
    print(ggbx)
    
    
    
  }
  dev.off()
  
  
  
  
  
  resultadosnumericPANCANCERforplotmut$mutornobyclassifier = resultadosnumericPANCANCERforplotmut$score > 10  
  chisq.test(table(resultadosnumericPANCANCERforplotmut$`Cancer Type`, resultadosnumericPANCANCERforplotmut$mutornobyclassifier))
  
  table(resultadosnumericPANCANCERforplotmut$`Mutation labels`,resultadosnumericPANCANCERforplotmut$mutornobyclassifier,resultadosnumericPANCANCERforplotmut$`Top Cancers`)
  resultadosnumericPANCANCERforplotmut = resultadosnumericPANCANCERforplotmut[order(resultadosnumericPANCANCERforplotmut$`Top Cancers`),]
  my_xlab <- paste("(N=",table(resultadosnumericPANCANCERforplotmut$`Top Cancers`),")",sep="")
  ggbx =ggplot(resultadosnumericPANCANCERforplotmut, aes(x =  `Top Cancers` , y =  score, fill = `Top Cancers`)) + geom_boxplot(varwidth = TRUE) + xlab("Cancer Type") + ylab("Score") +ggtitle(paste0("TCGA dataset mutation prediction score SJ classifier")) +
    scale_x_discrete(labels=my_xlab)
  
  myColors <- brewer.pal(7,"Set1")
  myColors <- brewer.pal(9,"Set1")
  myColors = c( myColors, "Cyan","Black")
  
  myColors <- brewer.pal(11,"Paired")
  myColors <- brewer.pal(11,"Set3")
  myColors <-c("#ebac23","#b80058","#008cf9","#006e00","#878500","#d163e6","#b24502","#ff9287","#5954d6","#00c6f8","#bdbdbd")
  ggbx = ggbx + scale_fill_manual(values = myColors, drop = FALSE) + scale_color_manual(values = myColors, drop = FALSE)+ theme_minimal()  + geom_hline(yintercept=10, linetype="dashed", color = "red") + 
    #    theme( axis.title.x = element_blank() ,axis.text.x = element_text(size = 20),
    #  axis.title.y = element_text(size = 20),plot.title = element_text(size = 20)) + 
    theme(text = element_text(size = 15)) 
  
  
  
  
  pdf("output/figures/Boxplot_Scores_TopCancers_Mutated_Ambiguous.pdf", width = 10)
  print(ggbx)
  dev.off()
  
  
  
  
  
  #resultspercancernothres = resultspercancer
  
  
  
  scoreandmutlocation = table(resultadosnumericPANCANCERforplot$`Mutation labels`, resultadosnumericPANCANCERforplot$score > scorediffthreshold)
  scoreandmutlocation = data.frame(scoreandmutlocation )
  colnames(scoreandmutlocation ) <- c("Mut_label", "Mut_predicted", "counts")
  scoreandmutlocation= scoreandmutlocation[scoreandmutlocation$Mut_label != "None",]
  
  pdf("output/figures/Accuracy_Barplot_Ambiguous.pdf", width = 10)
  p = ggplot(scoreandmutlocation , aes(x = Mut_label, y = counts, fill = Mut_predicted)) + geom_bar(stat="identity", position=position_dodge())
  print(p)
  dev.off()
  
  
  
  a = table(resultadosnumericPANCANCER$score > 0, resultadosnumericPANCANCER$`Mutation labels`!= "None",resultadosnumericPANCANCER$cancer )
  resultadosnumericPANCANCERbycancer  = resultadosnumericPANCANCER
  
  
  resultadosnumericPANCANCERbycancer$Mutpredicted = as.factor(resultadosnumericPANCANCER$score > 0)
  levels(resultadosnumericPANCANCERbycancer$Mutpredicted) <- c("No_Pred", "Pred")
  resultadosnumericPANCANCERbycancer$Mutobserved =  as.factor(resultadosnumericPANCANCER$`Mutation labels`!= "None")
  levels(resultadosnumericPANCANCERbycancer$Mutobserved) <- c("No_Obs", "Obs")
  aa = data.frame(table(resultadosnumericPANCANCERbycancer$Mutobserved,resultadosnumericPANCANCERbycancer$Mutpredicted,resultadosnumericPANCANCERbycancer$cancer ))
  cancers = unique(aa$Var3)
  TN = c()
  TP = c()
  FP = c()
  FN = c()
  for(cont in 1:length(cancers)){
    
    caso = aa[aa$Var3 == cancers[cont],]
    TN = c(TN, caso$Freq[1])
    FN = c(FN, caso$Freq[2])
    FP = c(FP, caso$Freq[3])
    TP = c(TP, caso$Freq[4])
    
  }
  resultspercancerPANCANCER = data.frame(cancers,
                                         Sensitivity = TP/(TP+FN),
                                         Specificity = TN/(TN+FP),
                                         Accuracy = (TP + TN )/ ( TP + TN + FP + FN),
                                         TP,
                                         TN,
                                         FP,
                                         FN
  )
  
  write.table(resultspercancerPANCANCER, "output/Accuracy_table_Ambiguous.tsv", sep = "\t", quote = F, row.names = F)
#}

