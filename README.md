# synple

An Easy-to-Use Python Wrapper for the Spectral Synthesis Code Synspec


** install **

Starting in the main synple directory

1- Compile synspec and rotin

 Make sure you have a fortran compiler -- gfortran, is expected and
 if you are using a different compiler you will need to modify the
 makefile accordingly.

 cd synspec

 make

 cd ..
 
2- Download the line list data

 Make sure you have wget and a working internet connection

 cd linelists

 make

 cd ..
 
3- Add the synple directory to your $PYTHONPATH enviromental variable (optional)

4- Make sure you have a working python install!

 
** test **

 Copy one of the model atmospheres in the 'models' folder (e.g. to your working
 directory) and test that the code works. For example if you use the MARCS model
 sun.mod and want to compute the solar spectrum between 6160. and 6164. AA 
 
   synple.py sun.mod 6160. 6164.
 
 or, similarly, from a python interpreter (e.g. ipython)
 
   from synple import syn

   x, y, z = syn('ksun.mod', (6160.,6164.))

   #and to plot the continuum normalized spectrum

   %pylab

   plot(x,y/z)

