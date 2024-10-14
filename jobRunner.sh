#! /usr/bin/bash
  
sbatch bcSub.sh

# wait until job submits outputs
while [ ! -f ./*.out ];
do
	echo "waiting"
	sleep 5
done

# allow time for all to come through
wait 

mv *.out outputs
cd outputs

for file in *.out; do
    mv -- "$file" "${file%.out}.txt"
done

# find average time for each iterayion
cd ../programs
python ./timer.py
