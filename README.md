Quantum Espresso 6.1 with the Environ 0.2 module modified to allow for the calculation of equilibrium semiconductor-solution interfaces. 


The example folder provides a series of inputs for a silicon slab in contact with water. 
It also contains scripts allowing for the rapid creation of several charges at once and the easy compilation of q-v.dat files.

Installation instructions: 

Preliminary steps:
     -1. Configure QE following the standard procedure (running './configure' from QE-main-dir should be enough in most cases).
      0. Compile QE without the Environ module (running 'make pw' from QE-main-dir).
If there are problems with the preliminary steps, look up for solutions on the PW forum or Quantum-Espresso webpage.
For PW and NEB:
      1. In QE-main-dir/ run the script addsonpatch.sh
        ./install/addsonpatch.sh Environ Environ/src Modules -patch
      2. In QE-main-dir/PW/src dir run  
        ../../Environ/patches/environpatch.sh -patch
      3. In QE-main-dir/install/ dir run 
        ./makedeps.sh -addson Modules Modules
        ./makedeps.sh PW/src
      4. In QE-main-dir re-compile pw.x and neb.x
                     make pw
                     make neb
      5. Run pw.x (or neb.x) with -environ flag, e.g. 
                     pw.x -environ < pwinput > pwoutput


Uninstallation   
It is possible to revert the Quantum-ESPRESSO package to a pristine (Environ-free) state by following these steps:
     1. In QE main dir run the script addsonpatch.sh
                      ./install/addsonpatch.sh Environ Environ/src Modules -revert
     2. In PW/src dir run 
               ../../Environ/patches/environpatch.sh -revert
     3. Be sure to remove objects, modules and executables
                        make clean

