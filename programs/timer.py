import numpy as np
import os

path = "./../outputs"

os.chdir(path)

times = np.empty(0)

# iterate through all file
for file in os.listdir():
    # Check whether file is in text format or not 
    if file.endswith(".txt"):
        f = open(f"{path}/{file}", "r")
        split = f.read().split()
        times = np.append(times, float(split[10]))

print("Average time taken for ", str(len(times)), " runs was ", str(times.mean()), " seconds.")
