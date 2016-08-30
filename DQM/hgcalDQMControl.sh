#!/bin/bash

DQMPROCESSINFO=$(ps aux | grep -v grep | grep runhgcalDQM.py)
DQMPROCESSID=$(echo $DQMPROCESSINFO | awk '{print $2}')

if [ -z $1 ]; then
    echo "Error: hgcalDQMControl.sh called without arguments. You can either give the command:"
    echo "./hgcalDQMControl.sh start"
    echo "to start the DQM processing, or"
    echo "./hgcalDQMControl.sh stop"
    echo "to stop the DQM processing."
    exit
fi

if [ $1 = "start" ]; then

    if [ -z $DQMPROCESSID ]; then
        echo "Starting hgcalDQM"
        cat runhgcalDQM_output.txt >> runhgcalDQM_output_old.txt
        python runhgcalDQM.py &> runhgcalDQM_output.txt &
    else
        echo "hgcalDQM already running; process ID is $DQMPROCESSID"
    fi
    
elif [ $1 = "stop" ]; then

    if [ -z $DQMPROCESSID ]; then
        echo "no process named hgcalDQM found to be running"
    else
        kill $DQMPROCESSID
        echo "Stopped process"
    fi

else
    
    echo "Invalid command: $1"
    echo "Error: hgcalDQMControl.sh called with meaningless arguments. You can either give the command :"
    echo "./hgcalDQMControl.sh start"
    echo "to start the DQM processing, or"
    echo "./hgcalDQMControl.sh stop"
    echo "to stop the DQM processing."

fi
