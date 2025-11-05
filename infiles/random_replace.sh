#!/bin/bash

num_in=1

for i in {11..11..1}
do

	input="test${num_in}.txt"
	inputname="test_${i}.txt"
	replace="1234"
	toreplace=`expr $replace + ${i} \* 2`
	echo $inputname
	echo $toreplace
	cp $input $inputname
#	root -l 'simc_replay_randomseed.cc\("$inputname},"test2.txt","${replace}","${toreplace}"\)'
#	root -l `printf "simc_replay_randomseed.cc(\"${inputname}\",\"test2.txt\",\"${replace}\",\"${toreplace}\")"`
	root -l `printf "simc_replay_randomseed.cc(\"${inputname}\",\"test2.txt\",\"1234\",\"5678\")"`
#	cp test2.txt $inputname

done

