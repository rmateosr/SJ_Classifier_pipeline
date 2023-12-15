#!/bin/bash
module load R/4.2
mkdir tmp
mkdir output
mkdir output/figures

# Specify the path of the required file (modify it accordingly)
input_file_folder="/rshare1/ZETTAI_path_WA_slash_home_KARA//home/rmateosr/SJproject/individualSJselectionprewilcoxonnormalized/"
metadata_path="mutmetadata_less_StrictKEAP1OL_rankedbyme_ambiguousintegrated.tsv"
annotated_folder="/rshare1/ZETTAI_path_WA_slash_home_KARA/home/rmateosr/SJproject/annotated/"
not_annotated_folder="/rshare1/ZETTAI_path_WA_slash_home_KARA/home/rmateosr/SJproject/notannotatedmorethan10/"

####   PIPELINE STARTS   #####
i=0
nsamples=$(ls -1 "$input_file_folder" | wc -l)
while [ $i -ne $nsamples ]
do
        i=$(($i+1))
        echo "$i"
        other_inputs=("$i")
        qsub -N Step7 -cwd -j yes -l s_vmem=64G script/Step7_Wilcoxon.sh "$input_file_folder" "$metadata_path" "$i"
done

qsub -N Step8 -cwd -j yes -l s_vmem=64G -hold_jid Step7 script/Step8_GatherSJs.sh

nsamples=$(ls -1 "$input_file_folder" | wc -l)
i=0
while [ $i -ne $nsamples ]
do
        i=$(($i+1))
        echo "$i"
        other_inputs=("$i")
        qsub -N Step9 -cwd -j yes -l s_vmem=64G -hold_jid Step8 script/Step9_Get_Normal.sh "$annotated_folder" "$i"
        qsub -N Step9 -cwd -j yes -l s_vmem=64G -hold_jid Step8 script/Step9_Get_Altered.sh "$not_annotated_folder" "$i"
done

qsub -N Step10 -cwd -j yes -l s_vmem=64G -hold_jid Step9 script/Step10_Matrix_Input.sh

qsub -N Step11_CV -cwd -j yes -l s_vmem=64G -hold_jid Step10 script/Step11_CV.sh $metadata_path
qsub -N Step11_Full -cwd -j yes -l s_vmem=64G -hold_jid Step10 script/Step11_Full_Model.sh $metadata_path

qsub -N Step12_Ambiguous -cwd -j yes -l s_vmem=64G -hold_jid Step11_Full script/Step12_Ambiguous.sh $metadata_path

qsub -N Figures_CV -cwd -j yes -l s_vmem=64G -hold_jid Step11_CV script/Figures_CV.sh $metadata_path
qsub -N Figures_Ambiguous -cwd -j yes -l s_vmem=64G -hold_jid Step12_Ambiguous script/Figures_Ambiguous.sh $metadata_path