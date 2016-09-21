#!/bin/bash
file=$1
dir=tmp
outDir=tmpOut

#STARTING SPILL READ AT TIME (1us): 0x01EE8310 RUN: 878 EVENT: 400
#Board header: on FMC-IO 9, trig_count in mem= 400, sk_status = 1
#Event header for event 0 with (200ns) timestamp 0x51FABA20  global tts(us)  0x00AA8338 and CKOV= 0

#EVENT: run spill  readtime nfmcio
#eventtime0 globaltts0 eventtime1 globaltts1
#sk0.0 sk0.1 sk1.0 sk1.1 etc


echo ${file}
if [  -e "${dir}" ];then 
	rm ${dir} -Rf
fi

mkdir ${dir}
if [ -e "${outDir}" ];then 
	rm ${outDir} -Rf
fi
mkdir ${outDir}


awk -v dir=${dir} -f f.awk $file


# take the list of spills and create the subdirs
run_fills=`grep STARTING $file | awk '{print $7, $9}'`
# take the list of boards and create the subdirs
echo "[INFO] Now reordering files" 


IFS=$'\n'
for run_fill in ${run_fills}
do
	let spillID=$spillID+1
	run=`echo $run_fill | awk '{print $2}'`
	spill=`echo ${run_fill} | awk '{print $1}'`

	dir=tmp/$spillID-$spill/
	events=`ls -1 ${dir}/*  | sort -n | uniq | grep dat | sed 's|.dat||'`
	boards=`ls -1 tmp/*/ |grep -v tmp | sort -n | uniq | sed '/^$/ d'`
	echo "Spill: $spill" 
	for event in ${events}
	do
		for board in ${boards}
		do
			f=${dir}/$board/$event.dat
			if [ ! -e "$f" ];then continue; fi
			cat $f >> ${outDir}/$spillID-$event.dat
		done
		cat ${outDir}/${spillID}-${event}.dat >> ${outDir}/${spillID}.dat
	done
	
done


# IFS=$'\n'

# 	if [ ! -d "${dir}/${run}/${board}" ];then
# 		mkdir ${dir}/${run}/${boardID}-${board} -p
# 	fi

# done




