#!/bin/bash

fichmin=$1
fichmax=$2

for fich in $(eval echo {$fichmin..$fichmax})
do
./scripts/fake_rolling_${fich}.sh &
done