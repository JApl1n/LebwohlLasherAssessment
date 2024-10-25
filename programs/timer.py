import numpy as np
import os

path = "./../outputs"

os.chdir(path)

times = np.empty(0)
progName = ""
L = ""
steps = ""

# iterate through all file
for file in os.listdir():
    # Check whether file is in text format or not 
    if file.endswith(".txt"):
        if file.startswith("slurm"):
            f = open(f"{path}/{file}", "r")
            split = f.read().split()
            print(split)
            times = np.append(times, float(split[10]))

progName = str(split[0][11:-1])
L = str(split[2][:-1])
steps = str(split[4][:-1])

output = str("Average time taken for "+progName+" to complete "+steps+" steps on a grid of size "+L+"x"+L+" for "+str(len(times))+" iterations was "+str(times.mean())+" seconds. \n")
print(output)
f = open("summary.txt", "a")
f.write(output)
f.close()
