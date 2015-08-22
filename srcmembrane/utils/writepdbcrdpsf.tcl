# cgtools::utils -- 
#
# Write charmm crd pdb files and xplor psf files
# Zunjing Wang 
# Author: Zun-Jing Wang
# June-August 2008

namespace eval ::cgtools::utils {
    namespace export writecrd_charmm  
    namespace export writepdb_charmm 
    namespace export writepsf_xplor
}

# ::cgtools::utils::writepdb_charmm --
# This will write a PDB file that in the standard CHARMM format
proc ::cgtools::utils::writecrd_charmm { file topology args} {

    set options {
        {periodbox.arg  "0"    "whether use periodic boundary condition"}
        {computecomz.arg  "0"    "whether compute center of mass of sytem"}
    }
    set usage "Usage: writecrd_charmm file topology \[periodbox] "
    # Strip off optional arguments and store in params
    array set params [::cmdline::getoptions args $options $usage]

    set useperiod $params(periodbox)
    set adjustcomz $params(computecomz)

    if {  $useperiod > 0 } {
	#mmsg::send [namespace current] "using period boundary condition"
	set inst_boxl [setmd box_l]
    }
    if {  $adjustcomz > 0 } {
    	set comz_mem [compute_membrane_comz  $topology]
    }

    set f [open $file "w"]
    puts $f "*  AUTHOR: ZUNJING WANG"
    puts $f "*  CRD FILE:   CHARMM-FORMAT COORDINATE OF CG LIPID (IMPLICIT SOLVENT MODEL)"
    puts $f "*  T [setmd temp]"
    puts $f "*  boxl [setmd box_l]"
    puts $f "*"
    puts $f [format "%5d" [setmd n_part]]

    set imol 1
    foreach mol $topology {
  	set moltype [lindex $mol 0]
    	set typeinfo [matchtype $moltype]
  	set molname [lindex $typeinfo 1]
    	set partbondlists [lindex $typeinfo 2]
    	set partbondtypelists [lindex $typeinfo 3]
    	set partbondcharmmlists [lindex $typeinfo 4]
    
    	#set and write particle list
    	set beadlists [lindex $partbondlists 0]
    	set nbeads [llength $beadlists]
    
    	#set beadtypelists [lindex $partbondtypelists 0]
        #set itype_begin [lindex [lindex $beadtypelists 0] 0]
    	set beadcharmmlists [lindex $partbondcharmmlists 0]

    	if {  $useperiod > 0 } {
		#mmsg::send [namespace current] "using period boundary condition"
		set av_posmol { 0.0 0.0 0.0 }
		for { set b 0 } { $b < $nbeads } {incr b } {
	                #current particle
        	        set partnum [lindex $mol [ expr $b + 1] ]
	                set posvec [part $partnum print pos]
			lset av_posmol 0 [expr [lindex $av_posmol 0] + [lindex $posvec 0] ]
			lset av_posmol 1 [expr [lindex $av_posmol 1] + [lindex $posvec 1] ]
			lset av_posmol 2 [expr [lindex $av_posmol 2] + [lindex $posvec 2] ]
	        }
		set av_posmol_x [expr [lindex $av_posmol 0]/($nbeads*1.0)]
        	set av_posmol_y [expr [lindex $av_posmol 1]/($nbeads*1.0)]
        	set av_posmol_z [expr [lindex $av_posmol 2]/($nbeads*1.0)]
		set posmol_mov { 0.0 0.0 0.0 }
		#x direction
		while { [expr $av_posmol_x + [lindex $posmol_mov 0] ] < 0. } { 
			lset posmol_mov 0 [expr [lindex $posmol_mov 0] + [lindex $inst_boxl 0] ]
		}
		while { [expr $av_posmol_x + [lindex $posmol_mov 0] ] > [lindex $inst_boxl 0] } { 
			lset posmol_mov 0 [expr [lindex $posmol_mov 0] - [lindex $inst_boxl 0] ]
		}
		#y direction
		while { [expr $av_posmol_y + [lindex $posmol_mov 1] ] < 0. } { 
			lset posmol_mov 1 [expr [lindex $posmol_mov 1] + [lindex $inst_boxl 1] ]
		}
		while { [expr $av_posmol_y + [lindex $posmol_mov 1] ] > [lindex $inst_boxl 1] } { 
			lset posmol_mov 1 [expr [lindex $posmol_mov 1] - [lindex $inst_boxl 1] ]
		}
		#z direction
		while { [expr $av_posmol_z + [lindex $posmol_mov 2] ] < 0. } { 
			lset posmol_mov 2 [expr [lindex $posmol_mov 2] + [lindex $inst_boxl 2] ]
		}
		while { [expr $av_posmol_z + [lindex $posmol_mov 2] ] > [lindex $inst_boxl 2] } { 
			lset posmol_mov 2 [expr [lindex $posmol_mov 2] - [lindex $inst_boxl 2] ]
		}
   	}


	for { set b 0 } { $b < $nbeads } {incr b } {
        	#current particle
        	set partnum [lindex $mol [ expr $b + 1] ]
        	set partnumcharmm [expr $partnum + 1]

        	set partname [lindex [lindex $beadcharmmlists $b] 1]

		set posvec [part $partnum print pos]
		set posx [lindex $posvec 0]
		set posy [lindex $posvec 1]
		set posz [lindex $posvec 2]
    		if {  $adjustcomz > 0 } {
			set posz [expr $posz - $comz_mem + 0.5*[lindex $inst_boxl 2]]
	    	}
    		if {  $useperiod > 0 } {	
			set posx [expr $posx + [lindex $posmol_mov 0]]
			set posy [expr $posy + [lindex $posmol_mov 1]]
			set posz [expr $posz + [lindex $posmol_mov 2]]
		}
		set seqname "L$imol"
		set linecrd [format "%5d %4d %4s %-4s%10.5f%10.5f%10.5f %-4s %-6d 0.00000" $partnumcharmm $imol $molname $partname $posx $posy $posz $seqname $imol]
		puts $f $linecrd

	}
	incr imol		
    }
    close $f
}

# ::cgtools::utils::writepdb_charmm --
# This will write a PDB file that in the standard CHARMM format
proc ::cgtools::utils::writepdb_charmm { file topology args} {

    set options {
        {periodbox.arg  "0"    "whether use periodic boundary condition"}
        {computecomz.arg  "0"    "whether compute center of mass of sytem"}
    }
    set usage "Usage: writecrd_charmm file topology \[periodbox] "
    # Strip off optional arguments and store in params
    array set params [::cmdline::getoptions args $options $usage]

    set useperiod $params(periodbox)
    set adjustcomz $params(computecomz)

    if {  $useperiod > 0 } {
        #mmsg::send [namespace current] "using period boundary condition"
        set inst_boxl [setmd box_l]
    }
    if {  $adjustcomz > 0 } {
    	set comz_mem [compute_membrane_comz  $topology]
    }

    set f [open $file "w"]
    puts $f "REMARK  AUTHOR: ZUNJING WANG"
    puts $f "REMARK  PDB FILE:   CHARMM-FORMAT COORDINATE OF CG LIPID (IMPLICIT SOLVENT MODEL)"
    puts $f "REMARK  T [setmd temp]"
    puts $f "REMARK  boxl [setmd box_l]"

    set imol 1
    foreach mol $topology {
  	set moltype [lindex $mol 0]
    	set typeinfo [matchtype $moltype]
  	set molname [lindex $typeinfo 1]
    	set partbondlists [lindex $typeinfo 2]
    	set partbondtypelists [lindex $typeinfo 3]
    	set partbondcharmmlists [lindex $typeinfo 4]
    
    	#set and write particle list
    	set beadlists [lindex $partbondlists 0]
    	set nbeads [llength $beadlists]
    
    	#set beadtypelists [lindex $partbondtypelists 0]
        #set itype_begin [lindex [lindex $beadtypelists 0] 0]
    	set beadcharmmlists [lindex $partbondcharmmlists 0]

    	if {  $useperiod > 0 } {
		#mmsg::send [namespace current] "using period boundary condition"
		set av_posmol { 0.0 0.0 0.0 }
		for { set b 0 } { $b < $nbeads } {incr b } {
	                #current particle
        	        set partnum [lindex $mol [ expr $b + 1] ]
	                set posvec [part $partnum print pos]
			lset av_posmol 0 [expr [lindex $av_posmol 0] + [lindex $posvec 0] ]
			lset av_posmol 1 [expr [lindex $av_posmol 1] + [lindex $posvec 1] ]
			lset av_posmol 2 [expr [lindex $av_posmol 2] + [lindex $posvec 2] ]
	        }
		set av_posmol_x [expr [lindex $av_posmol 0]/($nbeads*1.0)]
        	set av_posmol_y [expr [lindex $av_posmol 1]/($nbeads*1.0)]
        	set av_posmol_z [expr [lindex $av_posmol 2]/($nbeads*1.0)]
		set posmol_mov { 0.0 0.0 0.0 }
		#x direction
		while { [expr $av_posmol_x + [lindex $posmol_mov 0] ] < 0. } { 
			lset posmol_mov 0 [expr [lindex $posmol_mov 0] + [lindex $inst_boxl 0] ]
		}
		while { [expr $av_posmol_x + [lindex $posmol_mov 0] ] > [lindex $inst_boxl 0] } { 
			lset posmol_mov 0 [expr [lindex $posmol_mov 0] - [lindex $inst_boxl 0] ]
		}
		#y direction
		while { [expr $av_posmol_y + [lindex $posmol_mov 1] ] < 0. } { 
			lset posmol_mov 1 [expr [lindex $posmol_mov 1] + [lindex $inst_boxl 1] ]
		}
		while { [expr $av_posmol_y + [lindex $posmol_mov 1] ] > [lindex $inst_boxl 1] } { 
			lset posmol_mov 1 [expr [lindex $posmol_mov 1] - [lindex $inst_boxl 1] ]
		}
		#z direction
		while { [expr $av_posmol_z + [lindex $posmol_mov 2] ] < 0. } { 
			lset posmol_mov 2 [expr [lindex $posmol_mov 2] + [lindex $inst_boxl 2] ]
		}
		while { [expr $av_posmol_z + [lindex $posmol_mov 2] ] > [lindex $inst_boxl 2] } { 
			lset posmol_mov 2 [expr [lindex $posmol_mov 2] - [lindex $inst_boxl 2] ]
		}
   	}

	for { set b 0 } { $b < $nbeads } {incr b } {
        	#current particle
        	set partnum [lindex $mol [ expr $b + 1] ]
        	set partnumcharmm [expr $partnum + 1]

        	set partname [lindex [lindex $beadcharmmlists $b] 1]

		set posvec [part $partnum print pos]
		set posx [lindex $posvec 0]
		set posy [lindex $posvec 1]
		set posz [lindex $posvec 2]
    		if {  $adjustcomz > 0 } {
			set posz [expr $posz - $comz_mem + 0.5*[lindex $inst_boxl 2]]
	    	}
    		if {  $useperiod > 0 } {	
			set posx [expr $posx + [lindex $posmol_mov 0]]
			set posy [expr $posy + [lindex $posmol_mov 1]]
			set posz [expr $posz + [lindex $posmol_mov 2]]
		}
		set seqname "L$imol"
		
		set linepdb [format "ATOM  %5d %4s %4s %4d    %8.3f%8.3f%8.3f  1.00  0.00     %4s" $partnumcharmm $partname $molname $imol $posx $posy $posz $seqname]
		puts $f $linepdb

	}
	incr imol		
    }
    puts $f "TER"
    puts $f "END"
    close $f
}


# ::cgtools::utils::writepsf_xplor --
# This will write a PSF file that in the standard XPLOR format
proc ::cgtools::utils::writepsf_xplor { file topology } {

    # set particles, bonds, angles, dihedrals
    set partlinelist 0
    unset partlinelist
    set bondlinelist 0
    unset bondlinelist
    set angllinelist 0
    unset angllinelist
    set dihelinelist 0
    unset dihelinelist
    set imol 1
    foreach mol $topology {
  	set moltype [lindex $mol 0]
    	set typeinfo [matchtype $moltype]
  	set molname [lindex $typeinfo 1]
    	set partbondlists [lindex $typeinfo 2]
    	set partbondtypelists [lindex $typeinfo 3]
    	set partbondcharmmlists [lindex $typeinfo 4]
    
    	#set particle list
    	set beadlists [lindex $partbondlists 0]
    	set nbeads [llength $beadlists]
    	set beadtypelists [lindex $partbondtypelists 0]
        set itype_begin [lindex [lindex $beadtypelists 0] 0]
    	set beadcharmmlists [lindex $partbondcharmmlists 0]
	for { set b 0 } { $b < $nbeads } {incr b } {
        	#current positions of particles
        	set partnum [lindex $mol [ expr $b + 1] ]
        	set partnumcharmm [expr $partnum + 1]
		set parttype [lindex $beadlists $b]
                set parttypeinfo [lindex $beadtypelists [expr $parttype - $itype_begin]]
        	set parttypename [lindex $parttypeinfo 1]
        	set partmass [lindex $parttypeinfo 2]
        	set partcharge [lindex $parttypeinfo 3]
        	set partname [lindex [lindex $beadcharmmlists $b] 1]
		set seqname "L$imol"
		set linepart [format "%8d %-5s %-4d %4s %-4s %-4s%10.5f%14.4f           0   0.00000     -0.301140E-02" $partnumcharmm $seqname $imol $molname $partname $parttypename $partcharge $partmass ]
		#puts "$partnumcharmm $seqname $imol $molname $partname $parttypename $partcharge $partmass"
		lappend partlinelist $linepart
	}

    	#set bond list
    	set bondlists [lindex $partbondlists 1]
    	set nbonds [llength $bondlists]
    	for { set b 0 } { $b < $nbonds } {incr b } {
	        set curbond [lindex $bondlists $b ]
        	# index of the particles inside mol 
	        set partlists_inmol [lindex $curbond 1]
        	set part1_inmol [lindex $partlists_inmol 0]
   	     	set part2_inmol [lindex $partlists_inmol 1]
        	# partnumbers [1:npart] to link the bond 
  	      	set partnum1 [expr [lindex $mol [ expr $part1_inmol + 1]] + 1 ]
       	 	set partnum2 [expr [lindex $mol [ expr $part2_inmol + 1]] + 1 ]

		set linebond [format "%8d%8d" $partnum1 $partnum2]
		lappend bondlinelist $linebond
    	}

	
    	#set angle list
    	set angllists [lindex $partbondlists 2]
    	set nangls [llength $angllists]
    	for { set b 0 } { $b < $nangls } {incr b } {
        	set curbond [lindex $angllists $b ]
        	# index of the particles inside mol 
        	set partlists_inmol [lindex $curbond 1]
        	set part1_inmol [lindex $partlists_inmol 0]
        	set part2_inmol [lindex $partlists_inmol 1]
        	set part3_inmol [lindex $partlists_inmol 2]
        	# partnumbers [1:npart] to link the angle 
        	set partnum1 [expr [lindex $mol [ expr $part1_inmol + 1] ] + 1]
        	set partnum2 [expr [lindex $mol [ expr $part2_inmol + 1] ] + 1]
        	set partnum3 [expr [lindex $mol [ expr $part3_inmol + 1] ] + 1]
		
		set linebond [format "%8d%8d%8d" $partnum1 $partnum2 $partnum3]
		lappend angllinelist $linebond
    	}

    	#set dihedral list
	set dihelists [lindex $partbondlists 3]
    	set ndihes [llength $dihelists]
    	for { set b 0 } { $b < $ndihes } {incr b } {
        	set curbond [lindex $dihelists $b ]
       	 	# index of the particles inside mol 
	        set partlists_inmol [lindex $curbond 1]
        	set part1_inmol [lindex $partlists_inmol 0]
	        set part2_inmol [lindex $partlists_inmol 1]
        	set part3_inmol [lindex $partlists_inmol 2]
	        set part4_inmol [lindex $partlists_inmol 3]
        	# partnumbers [1:npart] to link the dihedral
	        set partnum1 [expr [lindex $mol [ expr $part1_inmol + 1] ] + 1]
        	set partnum2 [expr [lindex $mol [ expr $part2_inmol + 1] ] + 1]
	        set partnum3 [expr [lindex $mol [ expr $part3_inmol + 1] ] + 1]
        	set partnum4 [expr [lindex $mol [ expr $part4_inmol + 1] ] + 1]

		set linedihe [format "%8d%8d%8d%8d" $partnum1 $partnum2 $partnum3 $partnum4]
		lappend dihelinelist $linedihe
    	}
	incr imol		
    }

    
    # write titles, particles, bonds, dihedrals
    set f [open $file "w"]
    # write title 
    puts $f "PSF"
    puts $f ""
    set ntitle 1
    puts $f [format "%8d !NTITLE" $ntitle]
    puts $f "* PSF OF CG LIPID (IMPLICIT SOLVENT MODEL) GENERATED by ZUNJING WANG"
    puts $f ""
    # write particles 
    puts $f [format "%8d !NATOM" [setmd n_part]]
    foreach partline $partlinelist {
	puts $f $partline
    }
    puts $f ""

    # write bonds 
    set totbond [llength $bondlinelist]
    puts $f [format "%8d !NBOND: bonds" $totbond]
    set ibond 0
    foreach bondline $bondlinelist {
	puts $f $bondline nonewline
	incr ibond
	if { [expr $ibond % 4] == 0 || $ibond == $totbond} {
		puts $f ""
	}
    }
    puts $f ""

    # write angles 
    if { $nangls > 0} {
    	set totangl [llength $angllinelist]
    	puts $f [format "%8d !NTHETA: angles" $totangl]
    	set iangl 0
    	foreach anglline $angllinelist {
		puts $f $anglline nonewline
		incr iangl
		if { [expr $iangl % 3] == 0 || $iangl == $totangl} {
			puts $f ""
		}
    	}
    	puts $f ""
    }

    # write dihedrals 
    if { $ndihes > 0} {
    	set totdihe [llength $dihelinelist]
    	puts $f [format "%8d !NPHI: dihedrals" $totdihe]
    	set idihe 0
    	foreach diheline $dihelinelist {
		puts $f $diheline nonewline
		incr idihe
		if { [expr $idihe % 2] == 0 || $idihe == $totdihe} {
			puts $f ""
		}
    	}
    }
    close $f

}

