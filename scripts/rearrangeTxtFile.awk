BEGIN{


}

#STARTING SPILL READ AT TIME (1us): 0x01EE8310 RUN: 878 EVENT: 400
(/STARTING SPILL/){
	spillID+=1
	spillTime=$7
	run=$9
#	print "[DEBUG]" spillID, spillTime, run

}

(/Board header/){

	boardId=$5
	gsub(",", "", boardId)
	print "[DEBUG]", boardId
}




END{
	print spillID

}
