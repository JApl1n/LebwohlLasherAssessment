import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from os.path import isfile, join

# This takes the set of validations outputs from each method then compares how
# they evolve over time. The only two that are the same are numpy and the original.
# The reason for this is different amounts of random features are used in eachmethod.

path = "../validations"
fileNames = [f for f in listdir(path) if isfile(join(path, f))]

colNames = ["MC step", "Ratio", "Energy", "Order"]
fileNames

i = 0
order = np.empty((len(fileNames),51))
for f in fileNames:
    df = pd.read_csv(join(path, f), skiprows=8, delimiter=r"\s+", names=colNames, header=0, index_col=colNames[0])
    order[i][:] = df[colNames[-1]]
    i+=1
x = df.index

for i in range(len(fileNames)):
    plt.plot(x,order[i], label=fileNames[i])
plt.legend()
plt.title("Order over time using same starting seed.")
plt.xlabel("step")
plt.ylabel("order")
