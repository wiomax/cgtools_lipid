cgtools
=======
This software simulates lipid bilayer with a systematically solvent-free coarse-grained model. 

General information
-------------------

Portal: http://www.wiomax.com/cgtools
E-MAIL: info@wiomax.com

License description
-------------------

See the [Creative Commons Non-Commercial License 3.0] (https://creativecommons.org/licenses/by-nc/3.0/us/) for more details.

Please acknowledge the author(s) if you use this code in any way.

REFERENCES:
-----------

[1] Zun-Jing Wang, Markus Deserno. Systematic implicit solvent coarse-graining 
of bilayer membranes: lipid and phase transferability of the force field. 
New Journal of Physics, 2010, 12(9): 095004. 

[2] Zun-Jing Wang, Markus Deserno. A systematically coarse-grained solvent-free 
model for quantitative phospholipid bilayer simulations. 
Journal of Physical Chemistry B, 2010, 114(34): 11207-11220. 

PROGRAM :
---------

cgtools is now for building and simulating a membrane only. Its updated version 
(coming soon) will be a builder for a system with both membrane and peptide.

Initial positions of lipid molecules can be set with two ways:
1) Randomly fill up the whole simulation box 
2) Read from a configuration file 

There are two main outputs to the program :
      - ESPResSo - one can feed the membrane in the molecular dynamics (MD)
                software package. cgtools takes care of defining bonds, interactions,
		positions,  etc. All the relevant ESPResSo parameters may be set in the
		cgtools parameter file which will be read into the MD package. Then,
		ESPResSo runs a simulation on the membrane(s).

      - VMD - Visual Molecular Dynamics is a GUI that lets one visualize
        molecules in cgtools is able to write .pdb & .psf files that are
        read by VMD. 

cgtools is able to combine these two characteristics by applying MD
simulations on a membrane and output it in VMD as a function of simulation
time.


GETTING STARTED :
-----------------
System requirements :
Script language :   tcl 8.4 or higher
OS              :   Linux, Unix, Mac OS X (tested)
Additional soft :   ESPResSo <www.espresso.mpg.de/>
		    VMD <www.ks.uiuc.edu/Research/vmd/>

Set path for executable ESPResSo file and scripts:
For example:
If installing ESPResSo in the HOME directory, and the executable file "Espresso_bin" 
is in the folder of "$HOME/Espresso/obj-Core-pc-linux", one should add to '~/.bashrc':
export PATH=$HOME/Espresso/obj-Core-pc-linux/:$PATH
export ESPRESSO_SCRIPTS=~/Espresso/scripts"

INSTALLATION :
--------------
Written in Tcl, thus nothing need to be compiled.
Untar the archive in any directory. Set the CGTOOLS_DIR environment variable
to your shell, setting it to its installation path. For example, under bash:
If extracting the archive in the HOME directory, one should add to '~/.bashrc':
export CGTOOLS_DIR=~/cgtools/"
where the executable cgtoolsmain.tcl is located inside $CGTOOLS_DIR.
Don't forget to add this folder to your $PATH as well, and source your bash.

 * In order for cgtools to be used whithin ESPResSo: 
   Open your espresso configuration file '~/.espressorc' (create a new one if
   necessary) and add the following two lines:
   
   lappend auto_path "$env(CGTOOLS_DIR)/srcmembrane/"
   source $env(ESPRESSO_SCRIPTS)/parallel_tempering.tcl
   
   Also, for the force field to work properly, you should copy the
   'myconfig.h' file included in the home directory of cgtools to the home
   directory of ESPResSo. This will activate all (or almost all!) the
   necessary features of ESPResSo to run the simulations properly.


 * In order to use the 'loadseries' command in VMD, you must source a file to
   your ~/.vmdrc. To do so, copy the file 'vmd_plg.tcl' located in the
   'vmd_plugin/' directory wherever you want. Then append (or create) the init
   file ~/.vmdrc by using
   
   source [VMD_PLUGIN_DIRECTORY]/vmd_plg.tcl
   
   where you should replace [...] by the directory you chose for the
   'vmd_plg.tcl' file.


QUICK START :
-------------
cgtools will *not* run unless a configuration file is specified in argument.
The main script should be parsed inside Espresso : 

   > Espresso cgtoolsmain.tcl <CONFIG_FILE> 
	      [-new]
	      [-replica [-connect HOST]]
	      [-annealing]
	      [-annealfast]
	      [-hybrid]
	      [-ffs]

   where the brackets represent optional choices. 

   '-new' starts a new computation by deleting the existing folder if there is
   any, without this choice, the simulation will resume from the last checkpoint.

   '-replica' starts a parallel tempering simulation host. The different
   simulation temperatures are parameterized in the CONFIG_FILE. Adding the
   '-connect HOST' argument is required for slaves to connect to the host.

   '-annealing' starts an annealing simulation.

   '-annealfast' starts a fast annealing process, afterwards equilibrium the system 
   at the end temperature.

   '-hybrid' starts a hybrid simulation of Monte Carlo and Molecular dynamics.
   
   '-ffs' starts a forward-flux-sampling Molecular dynamics simulation.

   Note that the optional arguments can be alternated, and they do not require
   brackets. 

At the end, the script will have created a new directory with stored configuration 
files inside.


The 'configs/' folder contains several configuration files that are ready to
be ran.

PARAMETER FILE :
----------------
We refer the user to the 'configs/examples' folder which contains examples. 

INITIAL CONFIGURATION FILEs :
-----------------------------
We refer the user to the 'readfiles/' folder which contains working examples. 

FORCE FIELD :
-------------
Parameters for bonded interacting functions: set in 'overlapcoffs/' folder.
Nonbonded tabulated force field: set in 'forcetables/' folder.
For detailed information, read the online supporting information of reference [1].


STRUCTURE OF THE CODE :
-----------------------

The main script 'cgtoolsmain.tcl' is at the root of the archive. It mainly
reads the configuration file, and starts relevant routines for the output.

The 'configs/examples' folder contains example scripts that can be used to build
membranes. 

The 'srcmembrane/' folder contains all the source code for membrane simulations.
 
 * 'espresso/' - all routines that are related to the integration of the
   script in Espresso. From here will be ran the MD integrator to run a
   simulation.

 * 'generation/' - will build a membrane, generate topology and configurations. 

 * 'mmsg/' - a standard output wrapper. Easier to determine to which part of
   the code a given message belongs to. This was written by Ira Cooke. 

 * 'utils/' - contains lots of useful routines. Mainly math-geometry functions
   that are used to build a membrane, read and write files.

 * 'analysis/' - contains analysis routines for bilayer membrane.

TESTS & EXAMPLES:
-----------------

 * Compute RDFs, density profile of a POPC bilayer membrane, in the folder of "cgtools", 
   run a command line:  
	   mpirun -np 4 Espresso_bin cgtoolsmain.tcl -n 4 configs/examples/popcbilayer.tcl
   a output folder "popcbilayer_288" will be produced in the folder of cgtools.
   
 * Simulate a self-assembly process of a lipid system, in the folder of "cgtools", 
   run a command line:  
	   mpirun -np 4 Espresso_bin cgtoolsmain.tcl -n 4 configs/examples/random.tcl
   a output folder "random_288" will be produced in the folder of cgtools.

 * Simulate a DOPC bilayer membrane, in the folder of "cgtools", 
   run a command line:  
	   mpirun -np 4 Espresso_bin cgtoolsmain.tcl -n 4 configs/examples/dopcbilayer.tcl
   a output folder "dopcbilayer_288" will be produced in the folder of cgtools.

 * Simulate a DPPC bilayer membrane, in the folder of "cgtools", 
   run a command line:  
	   mpirun -np 4 Espresso_bin cgtoolsmain.tcl -n 4 configs/examples/dppcbilayer.tcl
   a output folder "dppcbilayer_288" will be produced in the folder of cgtools.

