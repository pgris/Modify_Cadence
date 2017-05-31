#!/bin/bash

thedir=../Ana_Cadence/OpSimLogs/WFD

for val in `ls -tr ${thedir}`
do 
fieldid=`echo ${val} | cut -d '_' -f3 | cut -d '.' -f1`
echo ${fieldid} ${val}
done