#!/bin/zsh


subjectsList=('004*ed.set' '005*ed.set' '006*ed.set' 'dic*ed.set' 'rou*ed.set')
N_SUBJECTS=${#subjectsList[@]}
numberFiles=$N_SUBJECTS

echo processing

for (( i=0;i<$ELEMENTS;i++)); do
    echo ${subjectsList[${i}]}
done 


logFileProgress='/users/aarslan/log_progression_jobBatchFeature.sh.1.txt'
logFile='/users/aarslan/log_list_files_jobRunClassif.sh.1.txt'
filesDirectory='/gpfs/data/tserre/data/neural_readout/results/'
eventField_1='correct'
eventField_2='mask'

normalization='minMax'
bin='6'
classifType='generalization'
regressorsEventField='type'
eventCode_1='1'
eventCode_2='0'
classifieur='fastLinearSVMAdrien'
fileSuffix=${normalization}'_bin'${bin}'_'${classifType}'_reg-'${regressorsEventField}'_sel-'${eventField_1}${eventCode_1}'-'${eventField_2}${eventCode_2}'_'${classifieur}'.mat'


##echo $indFeature 
for indSubject in ${subjectsList}
do
	fileName=${subjectsList[${indSubject}]}
	if ! [  -e ${fileName}.txt ]
	then
		echo this file ${fileName} is being processed > ${filesDirectory}/${fileName}.txt
		echo  now processing ....... ${fileName} 
		nohup matlab -nodesktop -r  "sillyscript() , runClassif(${indSubject},${indFeature})"
		echo ${fileName} .... 100% processed and saved	
		echo ${fileName} >> ${logFile}
		number_files_processed=$(cat ${logFile} | wc -l)
		percentage_progression=$((number_files_processed*100/numberFiles))
		echo percentage progression is ${percentage_progression} % and number files processed are $((number_files_processed)) / ${numberFiles} > ${logFileProgress}
	fi

done
