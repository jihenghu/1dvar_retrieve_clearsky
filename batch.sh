#!/bin/bash
first=$1
second=$2
while [ "$first" -le "$second" ]
do
    echo + [thread-main] $first  start ++++++++++++++++++++++++++++++++++
    ./main_clear_retrieve_landonly.exe $first	
    echo + [thread-main] $first finished
    let first=`date -d "-1 days ago ${first}" +%Y%m%d`
done
echo + [thread-main] all days done
