CGTools
=======
This software simulates lipid bilayer with a systematically solvent-free coarse-grained model. 

License description
-------------------

See the [Creative Commons Non-Commercial License 3.0] (https://creativecommons.org/licenses/by-nc/3.0/us/) for more details.

Please acknowledge the author(s) if you use this code in any way.

Quick start
-----------
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

The 'configs/' folder contains several configuration files that are ready to be ran.

References
----------
[1] Zun-Jing Wang, Markus Deserno. Systematic implicit solvent coarse-graining 
of bilayer membranes: lipid and phase transferability of the force field. 
New Journal of Physics, 2010, 12(9): 095004. 

[2] Zun-Jing Wang, Markus Deserno. A systematically coarse-grained solvent-free 
model for quantitative phospholipid bilayer simulations. 
Journal of Physical Chemistry B, 2010, 114(34): 11207-11220. 

General information
-------------------
Portal: http://www.wiomax.com/cgtools
E-MAIL: info@wiomax.com



