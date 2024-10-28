"""
Basic Python Lebwohl-Lasher code.  Based on the paper 
P.A. Lebwohl and G. Lasher, Phys. Rev. A, 6, 426-429 (1972).
This version in 2D.

Run at the command line by typing:

python LebwohlLasher.py <ITERATIONS> <SIZE> <TEMPERATURE> <PLOTFLAG>

where:
  ITERATIONS = number of Monte Carlo steps, where 1MCS is when each cell
      has attempted a change once on average (i.e. SIZE*SIZE attempts)
  SIZE = side length of square lattice
  TEMPERATURE = reduced temperature in range 0.0 - 2.0.
  PLOTFLAG = 0 for no plot, 1 for energy plot and 2 for angle plot.
  
The initial configuration is set at random. The boundaries
are periodic throughout the simulation.  During the
time-stepping, an array containing two domains is used; these
domains alternate between old data and new data.

SH 16-Oct-23
"""

import sys
import time
import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpi4py import MPI

#=======================================================================
def initdat(nmax):
    """
    Arguments:
      nmax (int) = size of lattice to create (nmax,nmax).
    Description:
      Function to create and initialise the main data array that holds
      the lattice.  Will return a square lattice (size nmax x nmax)
	  initialised with random orientations in the range [0,2pi].
	Returns:
	  arr (float(nmax,nmax)) = array to hold lattice.
    """
    arr = np.random.random_sample((nmax,nmax))*2.0*np.pi
    return arr
#=======================================================================
def plotdat(arr,pflag,nmax,figN):
    """
    Arguments:
	  arr (float(nmax,nmax)) = array that contains lattice data;
	  pflag (int) = parameter to control plotting;
      nmax (int) = side length of square lattice.
    Description:
      Function to make a pretty plot of the data array.  Makes use of the
      quiver plot style in matplotlib.  Use pflag to control style:
        pflag = 0 for no plot (for scripted operation);
        pflag = 1 for energy plot;
        pflag = 2 for angles plot;
        pflag = 3 for black plot.
	  The angles plot uses a cyclic color map representing the range from
	  0 to pi.  The energy plot is normalised to the energy range of the
	  current frame.
	Returns:
      NULL
    """
    if pflag==0:
        return
    u = np.cos(arr)
    v = np.sin(arr)
    x = np.arange(nmax)
    y = np.arange(nmax)
    cols = np.zeros((nmax,nmax))
    if pflag==1: # colour the arrows according to energy
        mpl.rc('image', cmap='rainbow')
        for i in range(nmax):
            for j in range(nmax):
                cols[i,j] = one_energy(arr,i,j,nmax)
        norm = plt.Normalize(cols.min(), cols.max())
    elif pflag==2: # colour the arrows according to angle
        mpl.rc('image', cmap='hsv')
        cols = arr%np.pi
        norm = plt.Normalize(vmin=0, vmax=np.pi)
    else:
        mpl.rc('image', cmap='gist_gray')
        cols = np.zeros_like(arr)
        norm = plt.Normalize(vmin=0, vmax=1)

    quiveropts = dict(headlength=0,pivot='middle',headwidth=1,scale=1.1*nmax)
    fig, ax = plt.subplots()
    q = ax.quiver(x, y, u, v, cols,norm=norm, **quiveropts)
    ax.set_aspect('equal')
    plt.show()
    plt.savefig("outputfig"+str(figN)+".png")

#=======================================================================
def savedat(arr,nsteps,Ts,runtime,ratio,energy,order,nmax):
    """
    Arguments:
	  arr (float(nmax,nmax)) = array that contains lattice data;
	  nsteps (int) = number of Monte Carlo steps (MCS) performed;
	  Ts (float) = reduced temperature (range 0 to 2);
	  ratio (float(nsteps)) = array of acceptance ratios per MCS;
	  energy (float(nsteps)) = array of reduced energies per MCS;
	  order (float(nsteps)) = array of order parameters per MCS;
      nmax (int) = side length of square lattice to simulated.
    Description:
      Function to save the energy, order and acceptance ratio
      per Monte Carlo step to text file.  Also saves run data in the
      header.  Filenames are generated automatically based on
      date and time at beginning of execution.
	Returns:
	  NULL
    """
    # Create filename based on current date and time.
    current_datetime = datetime.datetime.now().strftime("%a-%d-%b-%Y-at-%I-%M-%S%p")
    filename = "LL-Output-{:s}.txt".format(current_datetime)
    FileOut = open(filename,"w")
    # Write a header with run parameters
    print("#=====================================================",file=FileOut)
    print("# File created:        {:s}".format(current_datetime),file=FileOut)
    print("# Size of lattice:     {:d}x{:d}".format(nmax,nmax),file=FileOut)
    print("# Number of MC steps:  {:d}".format(nsteps),file=FileOut)
    print("# Reduced temperature: {:5.3f}".format(Ts),file=FileOut)
    print("# Run time (s):        {:8.6f}".format(runtime),file=FileOut)
    print("#=====================================================",file=FileOut)
    print("# MC step:  Ratio:     Energy:   Order:",file=FileOut)
    print("#=====================================================",file=FileOut)
    # Write the columns of data
    for i in range(nsteps+1):
        print("   {:05d}    {:6.4f} {:12.4f}  {:6.4f} ".format(i,ratio[i],energy[i],order[i]),file=FileOut)
    FileOut.close()
#=======================================================================
def one_energy(arr,ix,iy,nmax,leftColumn,rightColumn):
    """
    Arguments:
	  arr (float(nmax,nmax)) = array that contains lattice data;
	  ix (int) = x lattice coordinate of cell;
	  iy (int) = y lattice coordinate of cell;
      nmax (int) = side length of square lattice.
    Description:
      Function that computes the energy of a single cell of the
      lattice taking into account periodic boundaries.  Working with
      reduced energy (U/epsilon), equivalent to setting epsilon=1 in
      equation (1) in the project notes.
	Returns:
	  en (float) = reduced energy of cell.
    """
    en = 0.0
    ixp = (ix+1)%nmax # These are the coordinates
    ixm = (ix-1)%nmax # of the neighbours
    iyp = (iy+1)%nmax # with wraparound
    iym = (iy-1)%nmax #
#
# Add together the 4 neighbour contributions
# to the energy
#
    mid = arr[ix,iy]

    if (ix==arr.shape[0]-1): #if in final column of chunk use extra column given to it
        ang = mid-rightColumn[iy] 
    else: 
        ang = mid-arr[ixp,iy]

    en += 0.5*(1.0 - 3.0*np.cos(ang)**2)
    
    if (ix==0): #If in first column of chunk use left column given to it
        ang = mid-leftColumn[iy]
    else:
        ang = mid-arr[ixm,iy]

    en += 0.5*(1.0 - 3.0*np.cos(ang)**2)
    ang = mid-arr[ix,iyp]
    en += 0.5*(1.0 - 3.0*np.cos(ang)**2)
    ang = mid-arr[ix,iym]
    en += 0.5*(1.0 - 3.0*np.cos(ang)**2)
    return en
#=======================================================================
def all_energy(arr,nmax):
    """
    Arguments:
	  arr (float(nmax,nmax)) = array that contains lattice data;
      nmax (int) = side length of square lattice.
    Description:
      Function to compute the energy of the entire lattice. Output
      is in reduced units (U/epsilon).
	Returns:
	  enall (float) = reduced energy of lattice.
    """
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    numCols = nmax//size
    startCol = rank*numCols
    if rank != (size-1):
        endCol = (rank+1)*numCols
    else:
        endCol = nmax

    if rank == 0:
        leftColumn = arr[:,-1].copy() #column to left of beginning is last column due to boundary conditions
        rightColumn = arr[:,endCol+1].copy()
    elif rank == size-1:
        rightColumn = arr[:,0].copy() #column to right of final chunk is first column
        leftColumn = arr[:,startCol-1].copy()
    else:
        leftColumn = arr[:,startCol-1].copy() #send neighbour columns for one_energy
        rightColumn = arr[:,endCol+1].copy()


    enLocal = 0.0
    for i in range(startCol, endCol):
        for j in range(nmax):
            enLocal += one_energy(arr,i,j,nmax,leftColumn,rightColumn)
    enAll = comm.reduce(enLocal, op=MPI.SUM, root=0)
    
    if rank == 0:
        return enAll
    else:
        return None
#=======================================================================
def get_order(arr,nmax):
    """
    Arguments:
	  arr (float(nmax,nmax)) = array that contains lattice data;
      nmax (int) = side length of square lattice.
    Description:
      Function to calculate the order parameter of a lattice
      using the Q tensor approach, as in equation (3) of the
      project notes.  Function returns S_lattice = max(eigenvalues(Q_ab)).
	Returns:
	  max(eigenvalues(Qab)) (float) = order parameter for lattice.
    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    numCols = nmax//size
    startCol = rank*numCols
    if rank != (size-1):
        endCol = (rank+1)*numCols
    else:
        endCol = nmax

    QabLocal = np.zeros((3,3))
    delta = np.eye(3,3)
    #
    # Generate a 3D unit vector for each cell (i,j) and
    # put it in a (3,i,j) array.
    #
    lab = np.vstack((np.cos(arr),np.sin(arr),np.zeros_like(arr))).reshape(3,nmax,nmax
            )
    for a in range(3):
        for b in range(3):
            for i in range(startCol, endCol):
                for j in range(nmax):
                    QabLocal[a,b] += 3*lab[a,i,j]*lab[b,i,j] - delta[a,b]
    
    QabGlobal = comm.reduce(QabLocal, op=MPI.SUM, root=0)

    if rank == 0:
        QabGlobal /= (2*nmax*nmax)
        eigenvalues, _ = np.linalg.eig(QabGlobal)
        return eigenvalues.max()
    else:
        return None
#=======================================================================
def MC_step(arr,Ts,nmax):
    """
    Arguments:
	  arr (float(nmax,nmax)) = array that contains lattice data;
	  Ts (float) = reduced temperature (range 0 to 2);
      nmax (int) = side length of square lattice.
    Description:
      Function to perform one MC step, which consists of an average
      of 1 attempted change per lattice site.  Working with reduced
      temperature Ts = kT/epsilon.  Function returns the acceptance
      ratio for information.  This is the fraction of attempted changes
      that are successful.  Generally aim to keep this around 0.5 for
      efficient simulation.
	Returns:
	  accept/(nmax**2) (float) = acceptance ratio for current MCS.
    """
    #
    # Pre-compute some random numbers.  This is faster than
    # using lots of individual calls.  "scale" sets the width
    # of the distribution for the angle changes - increases
    # with temperaturie.

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    numCols = nmax//size
    startCol = rank*numCols
    if rank != (size-1):
        endCol = (rank+1)*numCols
    else:
        endCol = nmax

    leftNeighbour = (rank-1)%size
    rightNeighbour = (rank+1)%size

    localArr = arr[startCol:endCol,:].copy()
     
    if rank == 0:
        leftColumn = arr[-1,:].copy() #column to left of beginning is last column due to boundary conditions
        rightColumn = arr[endCol+1,:].copy() 
    elif rank == size-1:
        rightColumn = arr[0,:].copy() #column to right of final chunk is first column
        leftColumn = arr[startCol-1,:].copy()
    else:
        leftColumn = arr[startCol-1,:].copy() #send neighbour columns for one_energy
        rightColumn = arr[endCol+1,:].copy()
    
    
    scale = 0.1 + Ts
    xran = np.random.randint(0, high=localArr.shape[0], size=(localArr.shape[0],nmax))
    aran = np.random.normal(scale=scale, size=(localArr.shape[0],nmax))
  
    # im defining this function here so it has access to data within MC_Step
    def update_positions(offset):
        accept=0
        for i in range(offset,localArr.shape[0], 2):
            for j in range(nmax):
                ix = xran[i,j]
                ang = aran[i,j]

                en0 = one_energy(localArr,ix,j,nmax,leftColumn,rightColumn)
                localArr[ix,j] += ang
                en1 = one_energy(localArr,ix,j,nmax,leftColumn,rightColumn)

                if en1 <= en0 or np.exp(-(en1-en0)/Ts) >= np.random.uniform(0,1):
                    accept += 1
                else:
                    localArr[ix,j] -= ang
        return accept

    # Look at squares with even index
    acceptEven = update_positions(0)
    sync_boundaries(localArr, leftNeighbour, rightNeighbour, comm, nmax)
    # Look at squares with odd indices
    acceptOdd = update_positions(1)
    sync_boundaries(localArr, leftNeighbour, rightNeighbour, comm, nmax)

    totalAccept = comm.reduce(acceptEven + acceptOdd, op=MPI.SUM, root=0)
    if rank == 0:
        return totalAccept / (nmax*nmax)
    else:
        return None
#=======================================================================
def sync_boundaries(localArr, leftNeighbour, rightNeighbour, comm, nmax):

    # create contiguous buffers for sending and receiving boundary columns
    leftColSend = np.ascontiguousarray(localArr[0,:])
    rightColSend = np.ascontiguousarray(localArr[-1,:])
    leftColRecv = np.empty(nmax, dtype=localArr.dtype)
    rightColRecv = np.empty(nmax, dtype=localArr.dtype)

    # exchange with left neighbour
    comm.Sendrecv(sendbuf=rightColSend, dest=rightNeighbour, sendtag=1,
                  recvbuf=leftColRecv, source=leftNeighbour, recvtag=1)

    # exchange with right neighbour
    comm.Sendrecv(sendbuf=leftColSend, dest=leftNeighbour, sendtag=0,
                  recvbuf=rightColRecv, source=rightNeighbour, recvtag=0)

    localArr[0,:] = leftColRecv
    localArr[-1,:] = rightColRecv

#=======================================================================
def main(program, nsteps, nmax, temp, pflag):
    """
    Arguments:
	  program (string) = the name of the program;
	  nsteps (int) = number of Monte Carlo steps (MCS) to perform;
      nmax (int) = side length of square lattice to simulate;
	  temp (float) = reduced temperature (range 0 to 2);
	  pflag (int) = a flag to control plotting.
    Description:
      This is the main function running the Lebwohl-Lasher simulation.
    Returns:
      NULL
    """
    figN = int(0)
    # Create and initialise lattice
    lattice = initdat(nmax)
    # Plot initial frame of lattice
    #plotdat(lattice,pflag,nmax,figN)
    #figN += 1
    # Create arrays to store energy, acceptance ratio and order parameter
    energy = np.zeros(nsteps+1,dtype=np.dtype)
    ratio = np.zeros(nsteps+1,dtype=np.dtype)
    order = np.zeros(nsteps+1,dtype=np.dtype)
    # Set initial values in arrays
    energy[0] = all_energy(lattice,nmax)
    ratio[0] = 0.5 # ideal value
    order[0] = get_order(lattice,nmax)

    # Begin doing and timing some MC steps.
    initial = time.time()
    for it in range(1,nsteps+1):
        ratio[it] = MC_step(lattice,temp,nmax)
        energy[it] = all_energy(lattice,nmax)
        order[it] = get_order(lattice,nmax)
    final = time.time()
    runtime = final-initial
    
    # Final outputs
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank==0:
        print("{}: Size: {:d}, Steps: {:d}, T*: {:5.3f}: Order: {:5.3f}, Time: {:8.6f} s".format(program, nmax,nsteps,temp,order[nsteps-1],runtime))
        # Plot final frame of lattice and generate output file
        savedat(lattice,nsteps,temp,runtime,ratio,energy,order,nmax)
        #plotdat(lattice,pflag,nmax,figN)


#=======================================================================
# Main part of program, getting command line arguments and calling
# main simulation function.
#
if __name__ == '__main__':
    if int(len(sys.argv)) == 5:
        PROGNAME = sys.argv[0]
        ITERATIONS = int(sys.argv[1])
        SIZE = int(sys.argv[2])
        TEMPERATURE = float(sys.argv[3])
        PLOTFLAG = int(sys.argv[4])
        main(PROGNAME, ITERATIONS, SIZE, TEMPERATURE, PLOTFLAG)
    else:
        print("Usage: python {} <ITERATIONS> <SIZE> <TEMPERATURE> <PLOTFLAG>".format(sys.argv[0]))
#=======================================================================
