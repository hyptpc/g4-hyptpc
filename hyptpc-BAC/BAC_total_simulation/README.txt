2023.01.17

This is a simulation code for BAC.

Build this simulation code from the command line
--------------------------------
* Create a build directory as a subdirectory of this code

    mkdir build

* Run cmake from within the build directory, pointing to Gmsh's source directory, then run
  "make"

    cd build
    cmake ..
    make


How to start?
--------------------------------
* Execute the simulation with visualization, e.g.
    ./main

* Execute the simulation with root file that contains simualtion result, e.g.
    ./main run.mac {name of root file}.root {geometry version} {x position} {y position} {momentum of a particle}

* The explanation of the version of the geometry is in 'main.cc'.

* You can change the type of particle in 'include/BACPrimaryGeneratorAction.hh'. ('pion' or 'kaon')
    
    


