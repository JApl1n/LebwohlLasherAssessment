import pandas as pd
import matplotlib.pyplot as plt

fileName = "LL-Output-Sun-13-Oct-2024-at-04-28-41PM.txt"

colNames = ["MC step", "Ratio", "Energy", "Order"]
df = pd.read_csv(fileName, skiprows=8, delimiter=r"\s+", names=colNames, header=0, index_col=colNames[0]
        
index = df.index

y = 3
yAxis = colNames[y]
yData = df[yAxis]

plt.plot(MCStep, yData)
plt.grid(True)
plt.xlabel("MC Step")
plt.ylabel(yAxis)

plt.title(str(yAxis) + " over time for simulation")

# plt.savefig("initial"+str(yAxis)+"Fig.png"))
