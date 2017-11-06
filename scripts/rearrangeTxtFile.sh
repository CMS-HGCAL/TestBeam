#!/bin/bash
file=$1
tmpdir=/dev/shm/$USER
splitDir=$tmpdir/split
mergeDir=$tmpdir/merge
finalDir=$2

splitDir=$splitDir/`basename $file .txt`
mergeDir=$mergeDir/`basename $file .txt`
finalFile=$finalDir/`basename $file` || exit

#STARTING SPILL READ AT TIME (1us): 0x01EE8310 RUN: 878 EVENT: 400
#Board header: on FMC-IO 9, trig_count in mem= 400, sk_status = 1
#Event header for event 0 with (200ns) timestamp 0x51FABA20  global tts(us)  0x00AA8338 and CKOV= 0

#EVENT: run spill  readtime nfmcio
#eventtime0 globaltts0 eventtime1 globaltts1
#sk0.0 sk0.1 sk1.0 sk1.1 etc


echo "[`basename $0`] Processing file: ${file}"
if [ -e "${finalFile}" ];then
	echo "[`basename $0`] Final file ${finalFile} exists"
	exit 0
fi

if [ ! -e "${tmpdir}" ];then 
	mkdir -p $tmpdir/{split,merge}/ || exit 1
fi


awk -v baseDir=${splitDir} -f scripts/rearrangeTxtFile.awk $file || exit 1

IFS=$'\n'
for run in $splitDir/RUN_*
do
	
	run=`basename $run | sed 's|RUN_||'`
	spills=`ls --color=none -d $splitDir/RUN_${run}/SPILL_* | sed 's|.*SPILL_||' | sort -n`
	mkdir -p $mergeDir/RUN_${run}

	for spill in $spills
	do
		echo "[`basename $0`] Merging spillID: $spill"
		#tmp/RUN_000930/SPILL_01/SPILL_01-EVENT_000400-BOARD_23.txt
		#boards=`ls ${dir}/RUN_${run}/SPILL_${spill}/*.txt | sed 's|.*-BOARD_||;s|.txt||' | sort | uniq`
		
		events=`ls ${splitDir}/RUN_${run}/SPILL_${spill}/*.txt | sed 's|.*EVENT_\([0-9]*\)-BOARD_.*|\1|;s|.txt||' | sort | uniq`
		for event in $events
		do
			(paste $splitDir/RUN_${run}/SPILL_${spill}/SPILL_${spill}-EVENT_${event}-BOARD_* | sed -e 's|[[:space:]]RUN.*BOARD{,1}|\tBOARD|g' > $mergeDir/RUN_${run}/SPILL_${spill}-EVENT_${event}.txt) &
		done
		wait
	done
	cat $mergeDir/RUN_${run}/SPILL_*-EVENT_*.txt > $finalFile || exit 1
	
done

rm {$splitDir,$mergeDir}/ -Rf


exit 0
