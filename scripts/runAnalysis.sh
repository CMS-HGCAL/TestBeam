#!/bin/bash
usage(){
	echo "Usage: `basename $0` 
    	   --mode <beam|pedestal>
" 
}

dir=/tmp/shervin/

#------------------------------ parsing

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt -u -o hf: -l help,mode:,beam: -- "$@")
then
    # something went wrong, getopt will put out an error message for us
    exit 1
fi


set -- $options
#echo $options

while [ $# -gt 0 ]
do
    case $1 in
        -h|--help) usage; exit 0;;
		-f) FILE=$2; shift;;
        --mode) echo "[OPTION] MODE=$2"; MODE=$2; shift;;
        (--) shift; break;;
        (-*) usage; echo "$0: error - unrecognized option $1" 1>&2; usage >> /dev/stderr; exit 1;;
        (*) break;;
    esac
    shift
done


if [ -z "${FILE}" ];then
	echo "[ERROR] -f not specified at the command line: log file option is mandatory" >> /dev/stderr
	usage >> /dev/stderr
	exit 1
fi


if [ -z "${MODE}" ];then
	echo "[ERROR] --mode not specified at the command line: MODE option is mandatory" >> /dev/stderr
	usage >> /dev/stderr
	exit 1
fi

grep -i ${MODE} ${FILE}  | grep -v '#'

runs=`grep -i ${MODE} ${FILE} | grep -v '#' | cut -f 1`

for run in $runs
do
	echo $run
	find $dir -name "HGC_Output_${run}.txt"
done

