# Zun-Jing Wang 2008 Sep. 23 - Oct. 9
# 
# Simulate selfassembly of a POPC lipid system
# Parameter file for simulating single lipid molecule with only bonded interactions
# at fixed box size and temperature.
#
::mmsg::send [namespace current] "loading parameter file .. " nonewline
flush stdout

# Specify the name of the job <ident> the location of forcetables and
# the output and scripts directories, dir and filename of the initial positions
set currentrundir "runcode"
set numberoflipid 288
set ident "random_"
append ident $numberoflipid

set tabledir "./forcetables"
set overlapdir "./overlapcoffs"
set outputdir "./$ident"
set topofile "$ident.top"
set readpdbdir "./readfiles"
set readpdbname "singlemol.init.crd"
set rdfaadir "./rdfaa"
set rdfcgdir "./$ident/rdfcg"
set rdfcgintermoldir "./$ident/rdfcgintermol"
# Specify how we want to use vmd for visualization: allowable values
# are "interactive" "offline" "none".  interactive rarely works
set use_vmd "offline"

# --- Specify a global key for molecule types -----#
# In this line we specify the details of all molecule types in the
# system.  In this example molecule type 0 is assigned to be a lipid.
# The atom types in this lipid are respectively 0:CH, 1:PH, 2:GL, 3:E1, 4:E2, 5:AS, 6:AD, 7:AE. 
# BOND:
# The bonding interactions are respectively 0:CHPH, 1:PHGL, 2:GLE1, 3:GLE2, 4:E1AS, 5:E2AS, 6:ASAS, 7:ASAD, 8:ASAE. 
#
# BEND:
# Choos
# The bending interactions are respectively 0:CHPHGL, 1:PHGLE1,  2:PHGLE2, 3:GLE1AS, 4:GLE2AS, 5:E1ASAS, 6:E2ASAS, 7:ASASAS, 
#			       		    8:ASASAD, 9:ASADAS, 10:ASASAE. 
# Or
# The cosangl interactions are respectively 0:CHPHGL, 1:PHGLE1,  2:PHGLE2, 3:GLE1AS, 4:GLE2AS, 5:E1ASAS, 6:E2ASAS, 7:ASASAS, 
#			       		    8:ASASAD, 9:ASADAS, 10:ASASAE. 
#
# DIHE:
# The dihedral interactions are respectively 0:CHPHGLE2, 1:E2PHGLE1, 2:PHGLE1AS, 3:PHGLE2AS,  4:GLE1ASAS,  5:GLE2ASAS, 
#					     6:E1ASASAS, 7:E2ASASAD, 8:ASASASAS, 9:ASASADAS, 10:ASASASAE, 11:ADASASAE. 
#			.......
#		     }
# list of beads:{beadid}
set beadlist [list 0 1 2 3 5 5 5 5 7 4 5 5 6 5 5 7 ] 
# list of bonds:{bondid {bead1# bead2# in the bead lists}}
set bondlist [list { 0 { 0 1 } } { 1 { 1 2 } } { 2 { 2 3 } }  { 4 { 3 4 } } { 6 { 4 5 } } { 6 { 5 6 }} { 6 { 6 7 } } { 8 { 7 8 } } { 3 { 2 9 } } { 5 { 9 10 } } { 6 { 10 11 } } { 7 { 11 12 } } { 7 { 12 13 } } { 6 { 13 14 } } { 8 { 14 15 } } { 9 { 3 9 } } { 10 { 0 2 } } { 11 { 1 3 } } { 13 { 2 4 } } { 15 { 3 5 } } { 17 { 4 6 } } { 17 { 5 7 } } { 20 { 6 8 } } { 12 { 1 9 } } { 14 { 2 10 } } { 16 { 9 11 } } { 18 { 10 12 } } { 19 { 11 13 } } { 18 { 12 14 } } { 20 { 13 15 } } ]
# list of angls:{bondid {bead1# bead2# bead3# in the bead lists}}
set angllist [list ]
# list of dihes:{bondid {bead1# bead2# bead3# bead4# in the bead lists}} 
set dihelist [list ]
lappend popclist $beadlist
lappend popclist $bondlist
lappend popclist $angllist
lappend popclist $dihelist
unset beadlist
unset bondlist
unset angllist
unset dihelist

# list of beads names for each of types:{beadid beadname beadmass beadcharge}
set beadtypelist [list { 0 CH 87.166 1.1} { 1 PH  94.970 -1.2} { 2 GL 41.073 0.36} { 3 E1 58.036 -0.13} { 4 E2 58.036 -0.13 } { 5 AS 42.081 0.0} { 6 AD 26.038 0.0} { 7 AE 29.062 0.0} ]
# list of bond names for each of types:{bondid {bondnames}} 
set bondtypelist [list { 0 CH PH } { 1 PH GL } { 2 GL E1 } { 3 GL E2 } { 4 E1 AS } { 5 E2 AS } { 6 AS AS } { 7 AS AD } { 8 AS AE } { 9 E1 GL E2 } { 10 CH PH GL } { 11 PH GL E1 } { 12 PH GL E2 } { 13 GL E1 AS } { 14 GL E2 AS } { 15  E1 AS AS } { 16 E2 AS AS } { 17 AS AS AS } { 18 AS AS AD } { 19 AS AD AS } {20 AS AS AE } ]
# list of angl names for each of types:{bondid {anglnames}}  
set angltypelist [list ]
# list of dihe names for each of types:{bondid {dihenames}} 
set dihetypelist [list ]
# list of nonbonded interaction names for each of types:{{typeid1 typeid2} {nonbnames}} 
set nonbtypelist [list { { 0 0 } { CH CH } } { { 0 1 } { CH PH } } { { 0 2 } { CH GL } } { { 0 3 } { CH ES } } { { 0 4 } { CH ES } } { { 0 5 } { CH AS } } { { 0 6 } { CH AD } } { { 0 7 } { CH AE } } { { 1 1 } { PH PH } } { { 1 2 } { PH GL } } { { 1 3 } { PH ES } } { { 1 4 } { PH ES } } { { 1 5 } { PH AS } } { { 1 6 } { PH AD } } { { 1 7 } { PH AE } } { { 2 2 } { GL GL } } { { 2 3 } { GL ES } } { { 2 4 } { GL ES } } { { 2 5 } { GL AS } } { { 2 6 } { GL AD } }  { { 2 7 } { GL AE } } { { 3 3 } { ES ES } } { { 3 4 } { ES ES } } { { 3 5 } { ES AS } } { { 3 6 } { ES AD } } { { 3 7 } { ES AE } } { { 4 4 } { ES ES } } { { 4 5 } { ES AS } } { { 4 6 } { ES AD } } { { 4 7 } { ES AE } } { { 5 5 } { AS AS } } { { 5 6 } { AS AD } } { { 5 7 } { AS AE } } { { 6 6 } { AD AD } } { { 6 7 } { AD AE } } { { 7 7 } { AE AE } } ]
lappend popctypelist $beadtypelist
lappend popctypelist $bondtypelist
lappend popctypelist $angltypelist
lappend popctypelist $dihetypelist
lappend popctypelist $nonbtypelist
unset beadtypelist
unset bondtypelist
unset angltypelist
unset dihetypelist
unset nonbtypelist

# list of charmm beads names for each of types:{beadid cahrmmbeadname}
set charmmbeadlist [list { 0 CH } { 1 PH } { 2 GL } { 3 ES1 } { 4 AS11 } { 5 AS12 } { 6 AS13 } { 7 AS14 } { 8 AE15 } { 9 ES2 } { 10 AS21 } { 11 AS22 } { 12 AD23 } { 13 AS24 } { 14 AS25 } { 15 AE26 } ]
lappend popccharmmbeadlist $charmmbeadlist
unset charmmbeadlist

#moltypeid
lappend molpopclist "0" 
#molspec
lappend molpopclist "POPC" 
#
lappend molpopclist $popclist 
lappend molpopclist $popctypelist 
lappend molpopclist $popccharmmbeadlist 
unset popclist
unset popctypelist
unset popccharmmbeadlist

# moltypelists is used for setting global variable moltypeskey in generataion namespace
lappend moltypelists $molpopclist
unset molpopclist

#Notice rename E1 E2 to ES, total 28 types of nonbonded interactions

# --- Specify the system geometry and composition ----#
# Set the geometry 
# geometry structure:    geometry  "geometry characteristic" }
set readcrdfilename "\"random_" 
append readcrdfilename $numberoflipid
append readcrdfilename "/singlemol"
append readcrdfilename ".init.crd\""
set howtoset "\"random\""
set aspace " "
set geometryconponent "geometry" 
set geometry [set geometryconponent][set aspace][set howtoset][set aspace][set readcrdfilename] 
unset readcrdfilename 
unset howtoset
unset aspace
unset geometryconponent

# In this line we specify 360 of type 1 (ie lipid) are to be used
# n_molslist structure:    n_molslist  {  { 0 1 } }
set componet_n_molslist 0
lappend componet_n_molslist $numberoflipid
lappend vector_n_molslist $componet_n_molslist
set n_molslist "n_molslist"
lappend n_molslist $vector_n_molslist
unset componet_n_molslist
unset vector_n_molslist

#set n_molslist { n_molslist {  { 0 288 } } }

# Now bundle the above info into a list
lappend lipidspec $geometry
lappend lipidspec $n_molslist

# Now group the lipidspec with other specs into a list of such
# systems (we can have multiple systems if we like each with different
# composition of molecule types
lappend system_specs $lipidspec
unset lipidspec
unset geometry
unset n_molslist

# Set the box size: for pressure profile need a bigger z, i
# thus we set z=2*boxz, accordingly, we divide 2 in rdf computation.
set lengthx [expr sqrt(68.293696*$numberoflipid*0.5)]
lappend setbox_l $lengthx
lappend setbox_l $lengthx
lappend setbox_l $lengthx
unset lengthx

#don't compute line tension
set linetension 0

# Warmup parameters
#----------------------------------------------------------#
set warm_time_step 0.1

set free_warmsteps 0 
set free_warmtimes 1 

# ------ Integration parameters -----------------#
set main_time_step 0.1
set verlet_skin 2.0 
# -------Constant Temperature-----------------#
set langevin_gamma 0.2
set systemtemp 1.0
# -------Constant Pressure-----------------#
#set npt "on"
#set p_ext 0.000
#set piston_mass 0.0005
#set gamma_0 0.2
#set gamma_v 0.00004

# -------DPD-----------------#
set thermo "DPD"
set dpd_gamma 1.0 

# -------- Set the espresso integration steps and times 
# The number of steps to integrate with each call to integrate
set int_steps   100
# The number of times to call integrate
set int_n_times 1000

# --------- Set the replica exchange temperatures:  310 K --- 1.0 
set replica_temps { 0.9 0.95 1.0 1.05 }
# number of times to integrate between each replica MC step --- replace int_steps
set replica_timestep 100
# number of replica exchange rounds --- replace int_n_times
set replica_rounds 1000

# -------- Frequency of backup and analysis
# backup frequency 
set write_frequency 10 
# analysis frequency 
set analysis_write_frequency 100


# MPI distribution
# MPI distribution of bilayer
if { [lindex geometry 1] != "random"} {
        if {[setmd n_nodes]==2} {
            setmd node_grid 2 1 1
        } elseif {[setmd n_nodes]==4} {
            setmd node_grid 2 2 1
        } elseif {[setmd n_nodes]==6} {
            setmd node_grid 3 2 1
        } elseif {[setmd n_nodes]==8} {
            setmd node_grid 4 2 1
        } elseif {[setmd n_nodes]==16} {
            setmd node_grid 4 4 1
        } elseif {[setmd n_nodes]==32} {
            setmd node_grid 8 4 1
        } elseif {[setmd n_nodes]==64} {
            setmd node_grid 8 8 1
        }
}

# MPI distribution of random 
if { [lindex geometry 1] == "random"} {
        if {[setmd n_nodes]==8} {
            setmd node_grid 2 2 2
        } elseif {[setmd n_nodes]==64} {
            setmd node_grid 4 4 4
        }
}
# Potentials and Forces
set partbondnamelist [lindex [lindex $moltypelists 0] 3]
set beadnamelist [lindex $partbondnamelist 0 ]
set bondnamelist [lindex $partbondnamelist 1 ]
set anglnamelist [lindex $partbondnamelist 2 ]
set dihenamelist [lindex $partbondnamelist 3 ]
set nonbnamelist [lindex $partbondnamelist 4 ]
unset partbondnamelist

# BOND Potentials
#----------------------------------------------------------#
set nbond [llength  $bondnamelist]
for { set p 0 } { $p < $nbond } { incr p } {
	set bondname [lindex $bondnamelist $p]
	set bondfeature [llength $bondname]
	if {$bondfeature == 3} {
		set bondid [lindex $bondname 0] 
		set bn0 [lindex $bondname 1]
		set bn1 [lindex $bondname 2]
		set overlapnamenow overlap_bond.$bn0$bn1\.coff
		lappend overlapnames $overlapnamenow 
		lappend bonded_parms [list $bondid overlapped bond $outputdir/$overlapdir/$overlapnamenow]
	}
	if {$bondfeature == 4} {
		set bondid [lindex $bondname 0] 
		set bn0 [lindex $bondname 1]
		set bn1 [lindex $bondname 2]
		set bn2 [lindex $bondname 3]
		set overlapnamenow overlap_bend.$bn0$bn1$bn2\.coff
		lappend overlapnames $overlapnamenow 
		lappend bonded_parms [list $bondid overlapped bond $outputdir/$overlapdir/$overlapnamenow]
	}
}
unset nbond bondname bondid bn0 bn1 overlapnamenow
unset bondnamelist
#----------------------------------------------------------#

# ANGL Potentials
#----------------------------------------------------------#
set nangl [llength  $anglnamelist]
for { set p 0 } { $p < $nangl } { incr p } {
	set anglname [lindex $anglnamelist $p]
	set bondid [lindex $anglname 0] 
	set an0 [lindex $anglname 1]
	set an1 [lindex $anglname 2]
	set an2 [lindex $anglname 3]
	set overlapnamenow overlap_cosangl.$an0$an1$an2\.coff
	lappend overlapnames $overlapnamenow 
	lappend bonded_parms [list $bondid overlapped angle $outputdir/$overlapdir/$overlapnamenow]
}

if {$nangl > 0} {
	unset nangl anglname bondid an0 an1 an2 overlapnamenow
	unset anglnamelist
}
#----------------------------------------------------------#

# DIHE Potentials
#----------------------------------------------------------#
set ndihe [llength  $dihenamelist]
for { set p 0 } { $p < $ndihe } { incr p } {
	set dihename [lindex $dihenamelist $p]
	set bondid [lindex $dihename 0] 
	set dn0 [lindex $dihename 1]
	set dn1 [lindex $dihename 2]
	set dn2 [lindex $dihename 3]
	set dn3 [lindex $dihename 4]
	set overlapnamenow overlap_dihe.$dn0$dn1$dn2$dn3\.coff
	lappend overlapnames $overlapnamenow 
	lappend bonded_parms [list $bondid overlapped dihedral $outputdir/$overlapdir/$overlapnamenow]
}
if {$ndihe > 0} {
        unset ndihe dihename bondid dn0 dn1 dn2 dn3 overlapnamenow
        unset dihenamelist
}
#----------------------------------------------------------#


# Non Bonded Potentials 
#----------------------------------------------------------#
# Define the interactions between lipid beads
# NO Non-bonded interactions for this one molecular system

set nnonb [llength $nonbnamelist]
for { set i 0 } { $i < $nnonb } { incr i } {
	set nonbinfo [lindex $nonbnamelist $i]
	set btype [lindex $nonbinfo 0]
	set typei [lindex $btype 0]
	set typej [lindex $btype 1]
	set tabname [lindex $nonbinfo 1]
	set tabnamei [lindex $tabname 0]
	set tabnamej [lindex $tabname 1]
	set tablenamenow pot_force_nonbond.$tabnamei$tabnamej\.tab
	lappend tablenames $tablenamenow 
	lappend nb_interactions [list $typei $typej tabulated $outputdir/$tabledir/$tablenamenow]

	set rdfnamenow rdf.$tabnamei$tabnamej\.tab
#	lappend rdfnames $rdfnamenow
	lappend rdfaalist [list $typei $typej $rdfaadir/$rdfnamenow] 
	lappend rdfcglist [list $typei $typej $rdfcgdir/$rdfnamenow] 
	lappend rdfcgintermollist [list $typei $typej $rdfcgintermoldir/$rdfnamenow]
}
unset nnonb nonbinfo btype typei typej tabname tabnamei tabnamej tablenamenow 
unset beadnamelist  
unset nonbnamelist  

# set rdfoutputlist : merge the rdfs which have the smae name
set rdfoutputlist [list {0} {1} {2} {3 4} {5} {6} {7} {8} {9} {10 11} {12} {13} {14} {15} {16 17} {18} {19} {20} {21 22 26} {23 27} {24 28} {25 29} {30} {31} {32} {33} {34} {35}]

#::mmsg::send [namespace current] "nrdf = [llength $rdfoutputlist]"
#exit 1

# Analysis Parameters
#----------------------------------------------------------# 

# These are are parameters that will be passed to the setup_analysis
# command when it is called in the main.tcl script.

#For mode analysis, mgrid: # of grid, stray_cut_off: cutoff of below of up membrane
#Only work for a flat membrane
#Need check if there is implicit x_mgrid y_mgrid in code to permit box_x != box_y
set mgrid 8 
set stray_cut_off 30. 

# Use these flags to specify which observables to calculate during the
# simulation.  Values are calculated after every call to the espresso
# integrate command and written to files like
# time_vs_parametername. See the analysis package for more
# details
# Notice: fluctuations (mode analysis) is needed first for analysis_flags
#lappend analysis_flags orient_order
#lappend analysis_flags flipflop
#lappend analysis_flags stray
#lappend analysis_flags pressure 
#lappend analysis_flags stress_tensor 
lappend analysis_flags boxl
lappend analysis_flags temperature
lappend analysis_flags energy
#lappend analysis_flags fluctuations
#set profile_beadtypes [list  0 1 2 3 4 5 6 7 ]
#lappend analysis_flags "density_profile -beadtypes \{ $profile_beadtypes \} -nogrid"
#lappend analysis_flags stress_profile 
#set rminnow "0." 
#set rmaxnow "20.0"
#set nbinnow "200"
#lappend analysis_flags "rdf \{ $rdfcglist \} \{ $rdfoutputlist \} -rmin $rminnow -rmax $rmaxnow -nbin $nbinnow" 

::mmsg::send [namespace current] "Input paramter: done"
