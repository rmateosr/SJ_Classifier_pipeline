#!/usr/bin/Rscript
set.seed(12)   
args = commandArgs(trailingOnly=TRUE)
Sjscorethreshold = 10
scorediffthreshold= 10
c = 0.001
library(stringr)
library(RColorBrewer)
library(foreach)
library(ggplot2)
library(gridExtra)
library(data.table)
library(dplyr)
library(beepr)
library(parallel)
metadatafile = args[1]
pthresholds =  "p0001"
PANCANCERNFE2L2 = readRDS("sup_files/mc3.v0.2.8.PUBLIC_NRF2.maf.RDS")
resultadosnumericPANCANCER =   data.frame(fread("output/Scores_CV.tsv"))


resultadosnumericPANCANCER = unique(resultadosnumericPANCANCER)
resultadosnumericPANCANCER= resultadosnumericPANCANCER[!duplicated(resultadosnumericPANCANCER$TCGA_ID),]
colnames(resultadosnumericPANCANCER)[3] <- "Mutation labels"


####Relabeling of mutations#####

resultadosnumericPANCANCERforplot = resultadosnumericPANCANCER[order(resultadosnumericPANCANCER$score, resultadosnumericPANCANCER$`Mutation labels`, decreasing = T),]
resultadosnumericPANCANCERforplot$`Mutation labels` = as.character(resultadosnumericPANCANCERforplot$`Mutation labels`)
resultadosnumericPANCANCERforplot$`Mutation labels`[resultadosnumericPANCANCERforplot$`Mutation labels`=="NFE2L2_Exon2_Mut"] <- "NRF2_Exon2_Mut"
resultadosnumericPANCANCERforplot$`Mutation labels`[resultadosnumericPANCANCERforplot$`Mutation labels`=="NFE2L2_Other_Exon_Mut"] <- "NRF2_Other_Exon_Mut"
resultadosnumericPANCANCERforplot$`Mutation labels`[resultadosnumericPANCANCERforplot$`Mutation labels`=="Exon2_Skipping"] <- "NRF2_Exon2_Skipping"

old = c("None" ,"NRF2_Exon2_Mut",  "NRF2_Other_Exon_Mut" ,"KEAP1" , "CUL3"  ,"CNA" ,"NRF2_Exon2_Skipping", "Multihit")
old = c("None" ,"NRF2_Exon2_GOF_Mut",  "NRF2_Exon2_Mut" , "NRF2_Other_Exon_Mut" ,"KEAP1_Mut" ,"KEAP1 OL" ,"KEAP1 Ambiguous",  "CUL3"  ,"CNA" ,"NRF2_Exon2_Skipping", "Multihit")
new = c("None" ,"NRF2 Exon 2 Mut",  "NRF2 Other Exon Mut" ,"KEAP1 Mut" ,"KEAP1_OL" ,"KEAP1_Ambiguous", "CUL3 Mut"  ,"NRF2 CNA" ,"NRF2 Exon 2 Skipping", "Multihit")
new = c("None" ,"NRF2 Exon2 GOF Mut", "NRF2 Exon2 Mut" , "NRF2 Other Exon Mut" ,"KEAP1 Mut" , "KEAP1_OL" ,"KEAP1_Ambiguous", "CUL3 Mut"  ,"NRF2 CNA" ,"NRF2 Exon 2 Skipping", "Multihit")

new = c("None" ,"NRF2 Exon 2 Hotspot", "NRF2 Exon2 Mut" , "NRF2 Other Exon Mut" ,"KEAP1" , "KEAP1_OL" ,"KEAP1_Ambiguous", "CUL3 Mut"  ,"NRF2 CNA" ,"NRF2 Exon 2 Skip", "Multihit")



#Label multihit only with KEAP1 and NRF2 GOF labels
#currently only KEAP1 cases but might need to include GOF in the future
mutmetadata  =data.frame(fread(metadatafile))
mutmetadata[mutmetadata$TCGA_ID %in% resultadosnumericPANCANCERforplot$TCGA_ID[resultadosnumericPANCANCERforplot$`Mutation labels` == "Multihit"],]
resultadosnumericPANCANCERforplot =inner_join(resultadosnumericPANCANCERforplot, mutmetadata)

if(length(unique(resultadosnumericPANCANCERforplot$`Mutation labels`)) <=3){
  resultadosnumericPANCANCERforplot$`Mutation labels`[resultadosnumericPANCANCERforplot$`Mutation labels` == "Multihit" & resultadosnumericPANCANCERforplot$isKEAP1OL == "KEAP1_OL"] <- "KEAP1_OL"
  resultadosnumericPANCANCERforplot$`Mutation labels`[resultadosnumericPANCANCERforplot$`Mutation labels` == "Multihit" & resultadosnumericPANCANCERforplot$ismutexon2 ==  "NRF2_Exon2_Hotspot_GOF"] <- "NRF2_Exon2_Hotspot_GOF"
  resultadosnumericPANCANCERforplot$`Mutation labels`[resultadosnumericPANCANCERforplot$isKEAP1OL == "KEAP1_OL" & resultadosnumericPANCANCERforplot$ismutexon2 ==  "NRF2_Exon2_Hotspot_GOF"] <- "Multihit"
  
}



for(cont in 1:length(old)){
  resultadosnumericPANCANCERforplot$`Mutation labels`[resultadosnumericPANCANCERforplot$`Mutation labels` == old[cont]] = new[cont]
}

resultadosnumericPANCANCERforplot$`Mutation labels` =  factor(resultadosnumericPANCANCERforplot$`Mutation labels`, levels =  new[new %in% resultadosnumericPANCANCERforplot$`Mutation labels`])
#######






resultadosnumericPANCANCERforplot$is_mut_predicted = resultadosnumericPANCANCERforplot$score >scorediffthreshold
resultadosnumericPANCANCERforplot$is_mut_observed = resultadosnumericPANCANCERforplot$`Mutation labels` != "None"

resultadosnumericPANCANCERforplot$is_mut_predicted = factor(resultadosnumericPANCANCERforplot$is_mut_predicted, levels = c(FALSE, TRUE))


resultadosnumericPANCANCERforplot$is_mut_observed = factor(resultadosnumericPANCANCERforplot$is_mut_observed, levels = c(FALSE, TRUE))

resultstest = table(resultadosnumericPANCANCERforplot$is_mut_predicted,resultadosnumericPANCANCERforplot$is_mut_observed)
rownames(resultstest) <- c("Predicted Not Mutated", "Predicted Mutated")
colnames(resultstest) <- c("Observed Not Mutated", "Observed Mutated")

ggbx =ggplot(resultadosnumericPANCANCERforplot, aes(x = ((1:length(`Mutation labels`))-1) , y =  score, color = `Mutation labels`, fill = `Mutation labels`)) + geom_bar(stat="identity",position="dodge", width=0.02) +xlab("Sample") + ylab("Score") +ggtitle(paste0("TCGA dataset mutation prediction score SJ classifier"))+
  annotation_custom( tableGrob(resultstest), xmin=length(resultadosnumericPANCANCERforplot$score)/2-1, xmax=length(resultadosnumericPANCANCERforplot$score)[1]/2, ymin=max(resultadosnumericPANCANCERforplot$score)/2-1, ymax=max(resultadosnumericPANCANCERforplot$score/2))


myColors <- brewer.pal(7,"Set1")
myColors = c("Darkgray", myColors)
myColors = myColors[c(1,2,4,7)]
myColors[4] = "#FFD700"
ggbx = ggbx + scale_fill_manual(values = myColors, drop = FALSE) + scale_color_manual(values = myColors, drop = FALSE)+ theme_minimal()  + geom_hline(yintercept=10, linetype="dashed", color = "red")+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()  #remove y axis ticks
  )

#pdf(paste0("C:/Users/Raul/Documents/NRF2/SJclassifierversions/mutationburden2000_preambiguous//ggbox" , "_betabinomial_recount3_SJ",Sjscorethreshold,"score",scorediffthreshold,"NRF2_nonormalthres",nonormalthres,"newSRApipeline_",pthresholds[npvalue],scorediffthreshold,"_otherchrremoved_mutburden2000.pdf"), width = 20)
#ggbx
#dev.off()
ggsave(
  plot = ggbx,
  file = paste0("output/figures/Barplot_Scores_CV.pdf"),
  height = 10,
  width = 50 ,
  units = "cm",
  device =
    "pdf"
)

#ggplot(resultadosnumericPANCANCERforplot, aes(x = `Mutation labels` , y =  score, fill = `Mutation labels`)) + geom_violin( color = "black")+ scale_fill_manual(values = myColors, drop = FALSE) + scale_color_manual(values = myColors, drop = FALSE)+ theme_minimal()  + geom_hline(yintercept=10, linetype="dashed", color = "red")

pdf(paste0("output/figures/Boxplot_Scores_CV.pdf"), width = 15)

ggplot(resultadosnumericPANCANCERforplot, aes(x = `Mutation labels` , y =  score, fill = `Mutation labels`)) + geom_boxplot( color = "black")+ scale_fill_manual(values = myColors, drop = FALSE) + scale_color_manual(values = myColors, drop = FALSE)+ theme_minimal()  + geom_hline(yintercept=10, linetype="dashed", color = "red") + theme(text = element_text(size = 25))  
dev.off()

pdf(paste0("output/figures/Boxplot_Scores_CV_facet.pdf"), width = 12, height =10)

ggplot(resultadosnumericPANCANCERforplot[resultadosnumericPANCANCERforplot$`Mutation labels` != "None",]) + geom_boxplot(  aes(x = `Mutation labels` , y =  score, fill  =`Mutation labels` ),position = position_dodge2(preserve = "single"))+ scale_fill_manual(values = myColors, drop = FALSE) + scale_color_manual(values = myColors, drop = FALSE)+ theme_minimal()  +
 geom_hline(yintercept=10, linetype="dashed", color = "red") +
    facet_wrap(vars(cancer), nrow = 7)
dev.off()


resultadosnumericPANCANCERforplot2 = resultadosnumericPANCANCERforplot[resultadosnumericPANCANCERforplot$`Mutation labels` != "None",]


myColors <- brewer.pal(7,"Set1")
myColors = myColors[new[-1] %in% unique((resultadosnumericPANCANCERforplot2$`Mutation labels` ))]
myColors <- brewer.pal(7,"Set1")
myColors = c("Darkgray", myColors)
myColors = myColors[c(2,4,7)]
myColors[3] = "#FFD700"

if(length(unique(resultadosnumericPANCANCERforplot2$`Mutation labels`)) == length(c("NRF2 Exon2 Hotspot GOF", "KEAP1_OL"))){
  resultadosnumericPANCANCERforplot2$`Mutation labels` = factor(resultadosnumericPANCANCERforplot2$`Mutation labels`, levels = c("NRF2 Exon 2 Hotspot GOF", "KEAP1_OL"))
  myColors <- brewer.pal(7,"Set1")
  myColors = c(myColors)[c(1, 3)]
}

resultadosnumericPANCANCERforplot2$`Mutation labels` = factor(resultadosnumericPANCANCERforplot2$`Mutation labels`, levels =new[new %in% unique((resultadosnumericPANCANCERforplot2$`Mutation labels` ))])


ggbx1 =ggplot(resultadosnumericPANCANCERforplot2, aes(x = ((1:length(`Mutation labels`))-1) , y =  score, color = `Mutation labels`, fill = `Mutation labels`)) + geom_bar(stat="identity",position="dodge") +xlab("Sample") + ylab("Score") +ggtitle(paste0("TCGA dataset mutation prediction score SJ classifier"))+
  annotation_custom( tableGrob(resultstest,theme=ttheme_default(base_size = 20) ), xmin=length(resultadosnumericPANCANCERforplot2$score)/2-1, xmax=length(resultadosnumericPANCANCERforplot2$score)[1]/2, ymin=max(resultadosnumericPANCANCERforplot2$score)/2-1, ymax=max(resultadosnumericPANCANCERforplot2$score/2))
#ggbx = ggbx +scale_fill_nejm()


ggbx1 = ggbx1 + scale_fill_manual(values = myColors, drop = FALSE) + scale_color_manual(values = myColors, drop = FALSE)+ theme_minimal()  + geom_hline(yintercept=10, linetype="dashed", color = "red")+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),text = element_text(size = 20)) 

#ggbx = ggbx + scale_fill_nejm(drop = FALSE) + scale_color_nejm(drop = FALSE)


#ggbx1
pdf("output/figures/Barplot_Scores_CV_mutated.pdf", width = 30)
ggbx1
dev.off()









#stat="identity",position="dodge", width=0.02) +xlab("Sample") + ylab("Score") +ggtitle(paste0("TCGA dataset mutation prediction score SJ classifier"))+
#  annotation_custom( tableGrob(resultstest), xPANCANCERNFE2L2min=length(resultadosnumericPANCANCERforplot$score)/2-1, xmax=length(resultadosnumericPANCANCERforplot$score)[1]/2, ymin=max(resultadosnumericPANCANCERforplot$score)/2-1, ymax=max(resultadosnumericPANCANCERforplot$score/2))
#


cancers = unique(resultadosnumericPANCANCER$cancer)
cancers = cancers[order(cancers)]

resultadosnumericPANCANCERNFE2L2 = inner_join(resultadosnumericPANCANCER, PANCANCERNFE2L2)
resultadosnumericPANCANCERNFE2L2$Protein_position = as.numeric(resultadosnumericPANCANCERNFE2L2$Protein_position)
resultadosnumericPANCANCERNFE2L2only=resultadosnumericPANCANCERNFE2L2[resultadosnumericPANCANCERNFE2L2$NRF2,]
#TCGAandCancer = dataforbetabin[c("TCGA_ID_shorten","Cancer_Type")]
#colnames(TCGAandCancer) = c("TCGA_ID_shorten","cancer")
#resultspercancer = inner_join(resultadosnumericPANCANCERNFE2L2only, TCGAandCancer)
resultspercancer = resultadosnumericPANCANCERNFE2L2only


resultspercancer$mut_predicted = resultspercancer$score > scorediffthreshold
topfivecancer = names(table(resultspercancer$cancer))[order(table(resultspercancer$cancer),decreasing = T)[1:5]]

resultspercancer$main_cancer_types = resultspercancer$cancer
resultspercancer$main_cancer_types[!resultspercancer$main_cancer_types  %in% topfivecancer] <- "Others"

cancerlabels = (unique(resultspercancer$main_cancer_types))
cancerlabels = cancerlabels[order(cancerlabels)]
cancerlabels=cancerlabels[cancerlabels != "Others"]
cancerlabels = c(cancerlabels, "Others")
resultspercancer$main_cancer_types = factor(resultspercancer$main_cancer_types, levels = cancerlabels)
resultspercancer$Variant_Classification = factor(resultspercancer$Variant_Classification, 
                                                 levels = c("Frame Shift InDel",
                                                            "In Frame InDel",
                                                            "Missense Mutation",
                                                            "Nonsense Mutation",
                                                            "Splice Site"))
colnames(resultspercancer)[c(7:12)] <- c(  "Start Position",
                                           "End Position",
                                           "Protein Position", 
                                           "Variant Classification", 
                                           "Mutation Predicted", 
                                           "Cancer Types")

#resultspercancer = resultspercancer[resultspercancer$`Protein Position` %in%  c(79,
#                                                                                29,
#                                                                                82,
#                                                                                80,
#                                                                                81,
#                                                                                26,
#                                                                                30,
#                                                                                31,
#                                                                                34),]
#



resultspercancer = resultspercancer [resultspercancer$`Mutation labels` == "NRF2_Exon2_Hotspot_GOF",]
colnames(resultspercancer)[4] <- "Score"
resultspercancer$`Variant Classification` = factor(resultspercancer$`Variant Classification`, levels = unique(resultspercancer$`Variant Classification`))
####indels removed####3

resultspercancernoindel = resultspercancer[grep( "Mutation",resultspercancer$`Variant Classification`),]
pdf(paste0("output/figures/Hotspot_GOF_NRF2_noINDEL.pdf"))

xaxis = c( min(resultspercancer$`Protein Position`,na.rm = T),  max(resultspercancer$`Protein Position`,na.rm = T))
yaxis = c( min(resultspercancer$Score,na.rm = T),  max(resultspercancer$Score,na.rm = T))

p = ggplot(resultspercancernoindel, aes(x = `Protein Position`, y =Score, group= `Cancer Types`)) + geom_point(aes(color = `Cancer Types`, shape =`Variant Classification` )) + #xlim(xaxis) + ylim(yaxis) +
  geom_hline(yintercept=scorediffthreshold, linetype="dashed", color = "red")+ scale_color_brewer(palette="Set1") + theme_minimal()
print(p)

for(cont in 1:length(cancers)){
  resultspercancerLUAD = resultspercancer[resultspercancer$cancer ==cancers[cont],]
  if(dim(resultspercancerLUAD)[1]>0){
    p = ggplot(resultspercancerLUAD, aes(x = `Protein Position`, y = Score, group= `Cancer Types`)) + geom_point(aes(color =  "Black", shape =`Variant Classification` ))+ xlim(xaxis) + ylim(yaxis) + 
      geom_hline(yintercept=scorediffthreshold, linetype="dashed", color = "red") + theme_minimal()
    print(p)
  }
  
}

dev.off()

####

pdf(paste0("output/figures/Hotspot_GOF_NRF2_INDELincluded.pdf"))
xaxis = c( min(resultspercancer$`Protein Position`,na.rm = T),  max(resultspercancer$`Protein Position`,na.rm = T))
yaxis = c( min(resultspercancer$Score,na.rm = T),  max(resultspercancer$Score,na.rm = T))

p = ggplot(resultspercancer, aes(x = `Protein Position`, y =Score, group= `Cancer Types`)) + geom_point(aes(color = `Cancer Types`, shape =`Variant Classification` )) + #xlim(xaxis) + ylim(yaxis) +
  geom_hline(yintercept=scorediffthreshold, linetype="dashed", color = "red")+ scale_color_brewer(palette="Set1") + theme_minimal()
print(p)

for(cont in 1:length(cancers)){
  resultspercancerLUAD = resultspercancer[resultspercancer$cancer ==cancers[cont],]
  if(dim(resultspercancerLUAD)[1]>0){
    p = ggplot(resultspercancerLUAD, aes(x = `Protein Position`, y = Score, group= `Cancer Types`)) + geom_point(aes(color =  "Black", shape =`Variant Classification` ))+ xlim(xaxis) + ylim(yaxis) + 
      geom_hline(yintercept=scorediffthreshold, linetype="dashed", color = "red") + theme_minimal()
    print(p)
  }
  
}

dev.off()

#Not sure if necessary
write.table(resultspercancer, "output/NRF2_mutations.tsv", sep = "\t", quote = F, row.names = F) 


resultspercancernothres = resultspercancer



scoreandmutlocation = table(resultadosnumericPANCANCERforplot$`Mutation labels`, resultadosnumericPANCANCERforplot$score > scorediffthreshold)
scoreandmutlocation = data.frame(scoreandmutlocation )
colnames(scoreandmutlocation ) <- c("Mut_label", "Mut_predicted", "counts")
scoreandmutlocation= scoreandmutlocation[scoreandmutlocation$Mut_label != "None",]
#ggplot(scoreandmutlocation , aes(x = Mut_label, y = counts, fill = Mut_predicted)) + geom_bar(stat="identity", position=position_dodge())



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

write.table(resultspercancerPANCANCER, "output/Accuracy_table.tsv", sep = "\t", quote = F, row.names = F)




scoreandmutlocation = table(resultadosnumericPANCANCERforplot$`Mutation labels`, resultadosnumericPANCANCERforplot$score > scorediffthreshold)
scoreandmutlocation = data.frame(scoreandmutlocation )
colnames(scoreandmutlocation ) <- c("Mutation Labels", "Mutation Predicted", "Counts")
scoreandmutlocation= scoreandmutlocation[scoreandmutlocation$`Mutation Labels` != "None",]
pdf(paste0("output/figures/Accuracy_Barplot.pdf"), width  =12)

ggplot(scoreandmutlocation , aes(x = `Mutation Labels`, y = Counts, fill = `Mutation Predicted`)) + geom_bar(stat="identity", position=position_dodge()) + theme_minimal()


dev.off()

