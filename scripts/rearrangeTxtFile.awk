# directory structure:
# RUN/EVENT/BOARD
BEGIN{


}

#STARTING SPILL READ AT TIME (1us): 0x01EE8310 RUN: 878 EVENT: 400
(/STARTING SPILL/){
	spillID+=1
	spillTime=$7
	run=$9
#	print "[DEBUG SPILL]", spillID, spillTime, run
	dir=sprintf("%s/RUN_%06d", baseDir, run)
	system("mkdir -p "dir)
}

#Board header: on FMC-IO 9, trig_count in mem= 396, sk_status = 1
(/Board header/){
	boardID=$5
	gsub(",", "", boardID)

	triggerCount=$9
	gsub(",", "", triggerCount)

	skiStatus=$12
#	print "[DEBUG]", boardID, triggerCount, skiStatus
}

#Event header for event 0 with (200ns) timestamp 0xBEC20B15  global tts(us)  0x0072D911 and CKOV= 0
#DANGER: stale event warning!
(/Event header/){
	if(match($0, "DANGER")){
		eventID=$9
		eventTime=$13
		triggerTime=$16
	}else{
		eventID=$5
		eventTime=$9
		triggerTime=$12
	}
	
	#print "[DEBUG EVENT HEADER]", eventID, eventTime, triggerTime
	file=sprintf("SPILL_%02d-EVENT_%06d-BOARD_%02d.txt", spillID, eventID, boardID)
	printf("RUN=%06d\tSPILL=%02d\tEVENT=%06d\tGLOBALTIME=%s\tBOARD=%02d\n", run, spillID, eventID, spillTime, boardID) >> dir"/"file
	print "T "triggerTime, eventTime >> dir"/"file
}

(NF==4){
	sk0=$3
	sk1=$4
	print sk0, sk1 >> dir"/"file

}

END{
#	print spillID

}
