#!/bin/sh
i=0
while [ $i -ne 33 ]
do
        i=$(($i+1))
        echo "$i"
	qsub -cwd -j yes -l s_vmem=128G pipelinemutburden_includingambiguous2023_gathernormalcountsindividualGATHEREDinner.sh ${i}
done
