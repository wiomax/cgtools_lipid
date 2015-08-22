# Routine for reading a pdb file 
#
# Author: Zun-Jing Wang
# Sep.29 2008 Done 

namespace eval ::cgtools::utils {
    namespace export readpdb
    namespace export readcrd
}

#::cgtools::generation::readpdb
# Reading a pdb file line by line
proc ::cgtools::utils::readpdb { readfile } {

    set linelist 0
    unset linelist

    # open readfile
    set errcode [ catch { set readchannel [open $readfile "r"] } ]
    if { $errcode } {
        ::mmsg::err [namespace current]  "couldn't open file $readfile"
    } else {
        ::mmsg::send [namespace current] "reading positions of particles from the pdb file $readfile "

        # read lines in pdb file 
        while {[gets $readchannel pdbline] >= 0} {
                if { [lindex $pdbline 0] == "ATOM" } {
                        lappend linelist $pdbline
                }
        }
        close $readchannel
    }
    return $linelist
}

#::cgtools::generation::readcrd
# Reading a crd file line by line and trasfer them pdb lines
proc ::cgtools::utils::readcrd { readfile } {

    set linelist 0
    unset linelist

    set flagsetboxl 0
    
    # open readfile
    set errcode [ catch { set readchannel [open $readfile "r"] } ]
    if { $errcode } {
        ::mmsg::err [namespace current]  "couldn't open file $readfile"
    } else {
        ::mmsg::send [namespace current] "reading positions of particles from the crd file $readfile "

        # read lines in crd file 
        while {[gets $readchannel crdline] >= 0} {
                if { [lindex $crdline 0] != "*" } {
                	if { [llength $crdline] != 1 } {
				scan $crdline "%5d%5d %4s %4s%f%f%f %4s %6d 0.00000" partnumcharmm imol molname partname posx posy posz seqname imol_1
				set pdbline [format "ATOM  %5d %4s %4s %4d    %8.3f%8.3f%8.3f  1.00  0.00     %4s" $partnumcharmm $partname $molname $imol $posx $posy $posz $seqname]
                        	lappend linelist $pdbline
    				#puts "length of crdline is: [llength $crdline]"
    				#puts "$pdbline"
			}
                } else {
                	if { [lindex $crdline 1] == "boxl" } {
				setmd box_l [lindex $crdline 2] [lindex $crdline 3] [lindex $crdline 4]
    				set flagsetboxl 1
			}
		}


        }
        close $readchannel
    }

    if { $flagsetboxl == 0 } {
        ::mmsg::err [namespace current]  "The crd file lack the box_l parameter list"
    }

    return $linelist
}


