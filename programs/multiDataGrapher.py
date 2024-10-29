import pandas as pd
import matplotlib.pyplot as plt

from os import listdir
from os.path import isfile, join

path = "../summaries"
fileNames = [f for f in listdir(path) if isfile(join(path, f))]

for fileName in fileNames:
    if fileName[-6:-4] == "60": # replace with 20, 40 or 60 as with size of grid used in simulations
        colNames = ["MC step", "Ratio", "Energy", "Order"]
        df = pd.read_csv(join(path, fileName), skiprows=8, delimiter=r"\s+", names=colNames, header=0, index_col=colNames[0])

        # Get metadata

        f = open(join(path, fileName), "r")
        metaData = f.read()[108:226].split()

        L = metaData[4][:3]
        # If not 3 digit number reduce length by 1
        if L[-1] == "x":
            L=L[:2]
        # If not 2 digit number reduce length by 1
        if L[-1] == "x":
            L=L[:1]

        numSteps = metaData[10]
        T = metaData[14]
        time = metaData[-1]

        f.close()

        index = df.index

        y = 3 #1 is ratio, 2 energy, 3 order
        yAxis = colNames[y]
        yData = df[yAxis]

        plt.plot(index, yData, label=str(fileName))
        plt.grid(True)
        plt.xlabel("MC Step")
        plt.ylabel(yAxis)
        plt.legend()
        plt.title(str(yAxis) + " over time for simulation with N="+str(fileName[-6:-4])+" and T*=1.0") # chnage title here
        # plt.title(str(yAxis) + " over time for simulation")

        
plt.savefig("order60.png")
