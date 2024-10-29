# LebwohlLasherAssessment

Project files of assessment 1 of Core Computing 2024 on speeding up a program following the Lebwohl Lasher model.

testFig1Init and 2 and 3 are the results of running the input values 1000 20 1 and then 1,2 or 3. These are different visual outputs with very insignificant time differences however their plotting can be skipped entirely for speed comparisons.

testFig1Fin and 2 and 3 show the final time shot of the output after as shown in the command line 1000 time steps. 

Within 'programs' is all the programs required to run the different methods, each name linking to a module used to alter the original only using that module. This allows for a direct comparison of each method.

-------

original, numpy and numba all are executed as follows:

python filename.py I N T flag

where I is the number of time steps to run through, N the length of one of the sides of the square grid, T the initial temperature and flag the flag number for graphing, leave as 0 for no graph output.

For running the cython version, it is already compiled with the cythonSetup file, then run using:

python cythonRunner.py I N T flag

If you wish to change something in a function in cythonLebwohlLasher.pyx, do so then do:

python setupiCython.py build_ext --inplace

to compile the c file then you can run the cython runner. And lastly for the MPI method use:

mpiexec -n 4 python mpiLebwohlLasher.py I N T flag

where you can replace 4 with the number of nodes you'd like to use.

-------

timer.py can be used while in the programs directory to average times taken of any outputs in the outputs folder. This means if you want to only time the last batch of outputs you must remove any other output times before generating the batch for timing. This is usually alright as the individual outputs dont have much information anyway.

To automatiocally do all of this jobSub.sh will use the bcSub.sh file to run the desired file in bcSub.sh a desired amount of times and with a desired amount of nodes as chosen there, then collate the resulting times for processing by the timer automatically to get an average run time.

Once you have a collection of time for various inputs in summary.txt, you can go to the programs folder and output a graphed version of this using speedComparison.py with 'python speedComparison'. I have manual selections in there from the data I collected and you would have to manually select indices for your own which may be easier to just empty summary.txt. I will also not I added a string of numbers spaced by spaces to use to access the columns for the dataframe so they would also need to be added to the top.

-------
