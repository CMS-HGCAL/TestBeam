#!/bin/bash
file=$1
tmpdir=/dev/shm/$USER/
splitDir=$tmpdir/split
mergeDir=$tmpdir/merge
finalDir=tmpOut

#STARTING SPILL READ AT TIME (1us): 0x01EE8310 RUN: 878 EVENT: 400
#Board header: on FMC-IO 9, trig_count in mem= 400, sk_status = 1
#Event header for event 0 with (200ns) timestamp 0x51FABA20  global tts(us)  0x00AA8338 and CKOV= 0

#EVENT: run spill  readtime nfmcio
#eventtime0 globaltts0 eventtime1 globaltts1
#sk0.0 sk0.1 sk1.0 sk1.1 etc


echo "[`basename $0`] Processing file: ${file}"
if [  -e "${tmpdir}" ];then 
	rm ${tmpdir}/* -Rf
fi
mkdir -p $tmpdir/{split,merge}/


awk -v baseDir=${splitDir} -f scripts/rearrangeTxtFile.awk $file


IFS=$'\n'
for run in $splitDir/RUN_*
do
	spills=`ls --color=none -d $splitDir/RUN_000930/SPILL_* | sed 's|.*SPILL_||' | sort -n`
	run=`basename $run | sed 's|RUN_||'`
	mkdir -p $mergeDir/RUN_${run}

	for spill in $spills
	do
		echo $spill
		#tmp/RUN_000930/SPILL_01/SPILL_01-EVENT_000400-BOARD_23.txt
		#boards=`ls ${dir}/RUN_${run}/SPILL_${spill}/*.txt | sed 's|.*-BOARD_||;s|.txt||' | sort | uniq`
		
		events=`ls ${splitDir}/RUN_${run}/SPILL_${spill}/*.txt | sed 's|.*EVENT_\([0-9]*\)-BOARD_.*|\1|;s|.txt||' | sort | uniq`
		for event in $events
		do
#		echo $spill " -- " $event
			paste $splitDir/RUN_${run}/SPILL_${spill}/SPILL_${spill}-EVENT_${event}-BOARD_* | sed -e 's|[[:space:]]RUN.*BOARD{,1}|\tBOARD|g' > $mergeDir/RUN_${run}/SPILL_${spill}-EVENT_${event}.txt
		done
	done
	cat $mergeDir/RUN_${run}/SPILL_*-EVENT_*.txt > $finalDir/RUN_$run.txt
	
done



exit 0
cat $dir/RUN_*/SPILL_*-EVENT_*-BOARD_*.txt > $outDir/RUN_
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




