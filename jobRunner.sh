#! /usr/bin/bash
  
sbatch bcSub.sh

# wait until job submits outputs
while [ ! -f ./*.out ];
do
	echo "waiting"
	sleep 5
done

# allow time for all to come through
sleep 1 

mv *.out outputs
cd outputs

for file in *.out; do
    mv -- "$file" "${file%.out}.txt"
done

echo "Ignore first warning message. Here are the times of the run, they have also been collated in summary.txt in outputs."

# find average time for each iterayion
cd ../programs

sleep 2

python ./timer.py

