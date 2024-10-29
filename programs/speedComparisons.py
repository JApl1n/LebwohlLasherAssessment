import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# This is genuinely the way that works for my workflow and theres not enough time to start from the beginning to align everything

fileName = "../outputs/summary.txt"

df = pd.read_csv(fileName, delimiter=" ",header=0)

names = df["4"]
steps = np.array(df["7"].astype("int64"))

length = np.array(df["14"].str[:3])
for i in range(len(length)):
    l = length[i]
    length[i] = int(l.split("x")[0])

time = np.array(df["19"].astype('float64'))

# Because the order of some of the tests is out of order I have manually selected regions
# for outputting. It works.

# Effect of grid size on time taken (length vs time)
# Number of iterations 100 and other factors constant
xOriginalLength = length[6:19]
yOriginalLength = time[6:19]
xNumpyLength = length[24:37]
yNumpyLength = time[24:37]
xNumbaLength = length[41:50]
yNumbaLength = time[41:50]
xCythonLength = length[54:63]
yCythonLength = time[54:63]
xMpiLength4 = length[67:76]
yMpiLength4 = time[67:76]
xMpiLength8 = length[76:85]
yMpiLength8 = time[76:85]
xMpiLength2 = length[95:104]
yMpiLength2 = time[95:104]


# Effect of number of time steps on time taken (iterations vs time)
# Size of grid 25x25 and other factors constant
xOriginalIterations = np.concatenate((steps[9:10],steps[20:24]))
yOriginalIterations = np.concatenate((time[9:10],time[20:24]))
xNumpyIterations = np.concatenate((steps[27:28],steps[37:41]))
yNumpyIterations = np.concatenate((time[27:28],time[37:41]))
xNumbaIterations = steps[50:54]
yNumbaIterations = time[50:54]
xCythonIterations = steps[63:67]
yCythonIterations = time[63:67]
xMpiIterations4 = steps[85:90]
yMpiIterations4 = time[85:90]
xMpiIterations8 = steps[90:95]
yMpiIterations8 = time[90:95]
xMpiIterations2 = steps[104:109]
yMpiIterations2 = time[104:109]

oCol = "blue"
npCol = "red"
nbCol = "orange"
cCol = "green"
mCol2 = "lightgrey"
mCol4 = "grey"
mCol8 = "black" 

# Plot different methods for grid size variation
plt.scatter(xOriginalLength[:11], yOriginalLength[:11], marker=".", color=oCol, label = "Original")
plt.scatter(xNumpyLength, yNumpyLength, marker=".", color=npCol, label = "Numpy")
plt.scatter(xNumbaLength, yNumbaLength, marker=".", color=nbCol, label = "Numba")
plt.scatter(xCythonLength, yCythonLength, marker=".", color=cCol, label = "Cython")
plt.scatter(xMpiLength2, yMpiLength2, marker=".", color=mCol2, label = "MPI (2 nodes)")
plt.scatter(xMpiLength4, yMpiLength4, marker=".", color=mCol4, label = "MPI (4 nodes)")
plt.scatter(xMpiLength8, yMpiLength8, marker=".", color=mCol8, label = "MPI (8 nodes)")


plt.xlabel("Grid length")
plt.ylabel("Average run time (s)")
plt.grid()
plt.legend()
plt.title("Effect of grid size on average run time for each method for 100 time steps")
plt.show()
plt.savefig("gridSizeEffect.png")
plt.close()


# Plot different methods for num iterations variations
plt.scatter(xOriginalIterations, yOriginalIterations, marker=".", color=oCol, label = "Original")
plt.scatter(xNumpyIterations, yNumpyIterations, marker=".", color=npCol, label = "Numpy")
plt.scatter(xNumbaIterations, yNumbaIterations, marker=".", color=nbCol, label = "Numba")
plt.scatter(xCythonIterations, yCythonIterations, marker=".", color=cCol, label = "Cython")
plt.scatter(xMpiIterations2, yMpiIterations2, marker=".", color=mCol2, label = "MPI (2 nodes)")
plt.scatter(xMpiIterations4, yMpiIterations4, marker=".", color=mCol4, label = "MPI (4 nodes)")
plt.scatter(xMpiIterations8, yMpiIterations8, marker=".", color=mCol8, label = "MPI (8 nodes)")

# plt.plot(x, y, 'o', t, fit(t), '-')
plt.xlabel("Grid length")
plt.ylabel("Average run time (s)")
plt.grid()
plt.legend()
plt.title("Effect of number of iterations on average run time for each method for grid length 25")
plt.show()
plt.savefig("iterationsSizeEffect.png")
plt.close()

# Now below compare numba method to 8 nodes MPI

xMpiBigger = length[109:114]
yMpiBigger = time[109:114]
xNumbaBigger = length[114:119]
yNumbaBigger = time[114:119]

def curveVals(x,y):

    a_guess = 1
    b_guess = 0.1
    c_guess = 1
    
    # Fit the function a * np.exp(b * t) + c to x and y
    popt, pcov = curve_fit(
        lambda t, a, b, c: a * np.exp(b * t) + c,
        x, y, p0=(a_guess, b_guess, c_guess)
    )
    
    # The optimised values of the parameters are
    return popt

# Get curve_fit then plot with scatter of actual speeds for grid sizes for both methods
[a,b,c] = curveVals(xMpiBigger, yMpiBigger)
x_fitted = np.linspace(np.min(xMpiBigger), np.max(xMpiBigger), 100)
y_fitted = a * np.exp(b * x_fitted) + c

plt.scatter(xMpiBigger, yMpiBigger, marker=".", color="blue", label = "MPI (8 nodes)")
plt.plot(x_fitted, y_fitted, '-', label=(f"Ae^(Bx)+c where A={round(a,2)}, B={round(b,4)}, c={round(c,2)}"))

[a,b,c] = curveVals(xNumbaBigger, yNumbaBigger)
x_fitted = np.linspace(np.min(xNumbaBigger), np.max(xNumbaBigger), 100)
y_fitted = a * np.exp(b * x_fitted) + c

plt.scatter(xNumbaBigger, yNumbaBigger, marker=".", color="red", label = "Numba")
plt.plot(x_fitted, y_fitted, '-', label=(f"Ae^(Bx)+c where A={round(a,2)}, B={round(b,4)}, c={round(c,2)}"))

plt.xlabel("Grid length")
plt.ylabel("Average run time (s)")
plt.grid()
plt.legend()
plt.title("Effect of grid size on average run time MPI vs. Numba")
plt.savefig("mpiVsNumba.png")
