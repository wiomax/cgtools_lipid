# ::cgtools::analysis::analyze_fluctuations --
# 
# Obtained from MBtools, modified by Zun-Jing Wang Apr. 2009
# Routine for calculating the power spectrum of fluctuations for a
# flat bilayer sheet.  Uses the "modes_2d" routine in espresso to
# calculate the height function and perform the fft.  
# 
# It also calculate the thickness spectrum.
# Basic output:
# \# q^2 <|h(q)|^2> <|t(q)|^2>
# where h(q) correspond to the height fluctuation,
# and t(q) corresponds to the thickness.
#
#
#
# A summary of what happens is as follows
#
# (1) do_analysis calls analyze_fluctuations which in turn calls
#  <analyze modes2d> in espresso.

# (2) The resulting raw mode output (raw output from fftw) is sorted
# into increasing mode (q) order and some correction factors are
# applied. This is done via the function <arrange_modes>

# (3) The result is stored in av_pow 

# (4) When print_averages is called, the values stored in av_pow are
# printed by the (perhaps poorly named) power_analysis
# function. Generally for a mode analysis we would only call
# print_averages right at the end so that the entire simulation is
# used as an average.  For this reason this analysis is best done on
# stored configuration files rather than during simulation.
#
# Calculating kappa
#
# In order to calculate the bending stiffness from the results of this
# analysis one needs to carefully choose the fitting range so that
# only the q^4 regime is used.  An example fitting gnuplot script is
# as follows:
#
# 
# f(x) = a+b*x**2
#
# fit [0.02:0.2] f(x) 'powav.dat' u 1:($2*$1**2) via a,b
#
# set log
#
# plot 'powav.dat' u 1:2,\
#     f(x)/x**2,\
#     a/x**2
#
# print "kappa = ", 1.0/a
# 
#
# See the references: 3beads model by M. Deserno
#
#

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::fluctuations {
    variable av_pow
    variable av_pow_i 0

    variable av_pow_cntfile


    variable rawmodes
    variable verbose

    namespace export setup_fluctuations
    namespace export analyze_fluctuations
    namespace export printav_fluctuations
    namespace export resetav_fluctuations
}

proc ::cgtools::analysis::fluctuations::resetav_fluctuations { } {
    # Do nothing because we want to average over the whole simulation
}


proc ::cgtools::analysis::fluctuations::printav_fluctuations { } {    
    variable av_pow_i
    if { $av_pow_i > 0 } {
	# Do a power analysis on the intermediate results which does its own printing
	power_analysis
    }
}

proc ::cgtools::analysis::fluctuations::setup_fluctuations { args } {
    global ::cgtools::analysis::mgrid
    global ::cgtools::analysis::straycutoff
    global ::cgtools::analysis::outputdir
    variable rawmodes
    variable verbose
    variable av_pow_ht
    variable av_pow_th
    variable av_pow_cntfile
    set options {
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_fluctuations verbose "
    array set params [::cmdline::getoptions args $options $usage]

    for { set r 0 } { $r < [expr $mgrid*($mgrid)] } { incr r } {
	lappend av_pow_ht 0.0
	lappend av_pow_th 0.0
    }
    #fold_allcoordinates
    set rawmodes "[analyze modes2d ]"    

    set av_pow_cntfile 0

     # checkwith if the file exists already 
     set newfile 0
     while { $newfile == 0} {
        if { [file exists "$outputdir/powav.dat.$av_pow_cntfile"] } {
                incr av_pow_cntfile
        } else {
                set newfile 1
        }
     }

     #puts "av_pow_cntfile = $av_pow_cntfile"
}

proc ::cgtools::analysis::fluctuations::analyze_fluctuations {  } {
    variable rawmodes
    variable av_pow_i
    variable av_pow_ht
    variable av_pow_th
    global ::cgtools::analysis::mgrid
    global ::cgtools::analysis::straycutoff
    variable verbose
    ::mmsg::send [namespace current] "analyzing fluctuations"
    # Obtain unprocessed results from espresso
    #fold_allcoordinates
    set errnum [ catch { set rawmodes "[analyze modes2d]" } ]
    if { !$errnum } {
	set boxlen [setmd box_l]
	# Assume a square membrane in xy plane
	set Lx [lindex $boxlen 1]
	
	#Calculate the power spectrum for this configuration
	set pow_tmp [arrange_modes $rawmodes $Lx $mgrid ]
	# Add this set of results to the accumulated data
	for { set k 0 } { $k < [llength $av_pow_ht] } {incr k } {
	    set tmpval_ht "[expr [lindex $av_pow_ht $k] + [lindex [lindex $pow_tmp $k] 0]]"
	    lset av_pow_ht $k $tmpval_ht
	    set tmpval_th "[expr [lindex $av_pow_th $k] + [lindex [lindex $pow_tmp $k] 1]]"
	    lset av_pow_th $k $tmpval_th
	}	
	incr av_pow_i
    } else {
	mmsg::warn [namespace current] "fluctuation analysis failed"
    }    
}

# ::cgtools::analysis::arrange_modes --
#
# Arranges rawmodes into a power spectrum
#
proc ::cgtools::analysis::fluctuations::arrange_modes { rawmodes box_l mgrid } {
     # The arrangement is done in two stages 
    set round1 [modes_analysis_half $rawmodes $box_l $mgrid]
    set final [expand_powspec $round1 $box_l $mgrid]
    return $final
}


# ::cgtools::analysis::modes_analysis_half --
#
# Performs magic corrections and sorting.  Bloody painful
#
proc ::cgtools::analysis::fluctuations::modes_analysis_half { rawmodes box_l mgrid } {
    set boxlen [setmd box_l]
    # Assume a square membrane in xy plane
    set Lx [lindex $boxlen 0]

    # This prefactor has the advantage of giving a
    # physical spectrum - it does not depend on gridding.
    # The L^2 factor of the spectrum is now incorporated.
    set norm [expr $Lx/[expr $mgrid*$mgrid*1.0]]

    for { set i  0 } { $i < $mgrid } { incr i } {
	set qi $i
	if { $qi > [expr $mgrid/2 - 1] } { set qi [expr $i - $mgrid] }
	set sincx [sinc [expr $qi/($mgrid*1.0)] ]
	set normx [expr 1.0/$sincx]
	set normx [expr $normx*$normx ]
	
	for { set j  0 } { $j < [expr $mgrid/2 + 1] } { incr j } {  
	    set qj $j
	    if { $qj > [expr $mgrid/2 - 1] } { set qj [expr $j - $mgrid] }
	    set sincy [sinc [expr $qj/($mgrid*1.0)]]
	    set normy [expr 1.0/$sincy]
	    set normy [expr $normy*$normy]
	    
	    # ht: height grid; th: thickness.
	    set ht_re [expr [lindex [lindex $rawmodes 1 $i $j] 0]*$norm]
	    set ht_im [expr [lindex [lindex $rawmodes 1 $i $j] 1]*$norm]
	    set th_re [expr [lindex [lindex $rawmodes 1 $i $j] 2]*$norm]
	    set th_im [expr [lindex [lindex $rawmodes 1 $i $j] 3]*$norm]

	    # Add the Deserno correction factor
	    set pow "[expr ($ht_re*$ht_re + $ht_im*$ht_im)*$normy*$normx] [expr ($th_re*$th_re + $th_im*$th_im)*$normy*$normx]" 
	    
	    lappend powspec $pow
	}
    }
    return $powspec
}

# ::cgtools::analysis::expand_powspec --
#
# More of the nastiness
#
proc ::cgtools::analysis::fluctuations::expand_powspec { rawmodes box_l mgrid } {	
    for { set i  0 } { $i < $mgrid } { incr i } {
	for { set j  0 } { $j < [expr $mgrid] } { incr j } {  
	    set xindex $i
	    set yindex $j
	    if { $j > [expr $mgrid/2] } { 
		if { $i !=  0 } {
		    set xindex [expr $mgrid - $i ]
		    set yindex [expr $mgrid - $j ]
		} else {
		    set yindex [expr $mgrid - $j ]
		}
	    }
	    set index [expr $xindex*($mgrid/2+1) + $yindex]
	    lappend powspec [lindex $rawmodes $index]
	}
    }
    return $powspec
}


# ::cgtools::analysis::sinc --
#
# Calculate the sinc(x) function required in order to apply the magic
# Deserno correction
#
proc ::cgtools::analysis::fluctuations::sinc { x } {	
    set retval 0
    if { $x != 0.0 } {
	set retval [expr sin($x*3.141592654)/($x*3.141592654) ]
    } else {
	set retval 1
    }
    return $retval
}
    
# Basically just creates a file of the squared modulus of the
# wavevector vs the squared power.  This is essentially just a sorting
# and printing function
proc ::cgtools::analysis::fluctuations::power_analysis { } {
    global ::cgtools::analysis::mgrid
    global ::cgtools::analysis::outputdir
    variable av_pow_ht
    variable av_pow_th
    variable av_pow_i
    variable av_pow_cntfile

    set file "$outputdir/powav.dat.$av_pow_cntfile"
    set Lx [lindex [setmd box_l] 0]

    set f_powav [open $file w]
    puts $f_powav "\# q^2 <|h(q)|^2> <|t(q)|^2>"
    set wavevectors [modes_q2 $Lx $mgrid]
    set k 0
    foreach val $av_pow_ht {
	puts $f_powav "[lindex $wavevectors $k] [ expr $val/($av_pow_i*1.0) ] [expr [lindex $av_pow_th $k]/($av_pow_i*1.0)]"
	incr k
    }
    close $f_powav

    incr av_pow_cntfile
    
}

proc modes_q2 { box_l mgrid } {

    for { set i  0 } { $i < $mgrid } { incr i } {
        for { set j  0 } { $j < [expr $mgrid] } { incr j } {
            set qi $i
            set qj $j
            #           puts -nonewline "raw $qi $qj "     
            if { $qi > [expr $mgrid/2 - 1] } { set qi [expr $i - $mgrid] }
            if { $qj > [expr $mgrid/2 - 1] } { set qj [expr $j - $mgrid] }

            set qx [expr $qi*2*3.141592/$box_l]
            set qy [expr $qj*2*3.141592/$box_l]
            #           puts -nonewline "$qi $qj "          
            set q2 [expr $qx*$qx + $qy*$qy]
            lappend wavecs $q2
            #           puts "$q2 "
        }
    }

    return $wavecs
}


### To fold all the coordinates from periodic boundary condition
proc ::cgtools::analysis::fluctuations::fold_allcoordinates { } {
		variable ::cgtools::analysis::topology

		set inst_boxl [setmd box_l]
		set inst_boxx [lindex $inst_boxl 0]
		set inst_boxy [lindex $inst_boxl 1]
		puts "inst_boxx : $inst_boxx"
		puts "inst_boxy : $inst_boxy"

   		foreach mol $topology {
   			set moltype [lindex $mol 0]
		    	set typeinfo [::cgtools::utils::matchtype $moltype ]
  
			set partbondlists [lindex $typeinfo 2]

			set beadlists [lindex $partbondlists 0]
			set nbeads_mol [llength $beadlists]

			set center_x 0.
			set center_y 0.

			for { set b 0 } { $b < $nbeads_mol } {incr b } {
			        set partnum [lindex $mol [ expr $b + 1] ]

        			#current positions of particles
				set curpos [part $partnum print id pos]
			        set posx [lindex $curpos 1]
			        set posy [lindex $curpos 2]

			        #geometry center of molecule
				set center_x [expr $center_x + $posx]
				set center_y [expr $center_y + $posy]
			}
			
			set center_x [expr $center_x/$nbeads_mol]
			set center_y [expr $center_y/$nbeads_mol]
		puts "center_x : $center_x"
		puts "center_y : $center_y"

			#periodic boundary
			set dx 0.
			while {$center_x < 0.} {	
				set center_x [expr $center_x + $inst_boxx ]
				set dx [expr $dx + $inst_boxx ]
			}
			while {$center_x >= $inst_boxx} {	
				set center_x [expr $center_x - $inst_boxx ]
				set dx [expr $dx - $inst_boxx ]
			}
			set dy 0.
			while {$center_y < 0.} {	
				set center_y [expr $center_y + $inst_boxy ]
				set dy [expr $dy + $inst_boxy ]
			}
			while {$center_y >= $inst_boxy} {	
				set center_y [expr $center_y - $inst_boxy ]
				set dy [expr $dy - $inst_boxy ]
			}
		puts "center_x : $center_x"
		puts "center_y : $center_y"

			for { set b 0 } { $b < $nbeads_mol } {incr b } {

			        set partnum [lindex $mol [ expr $b + 1] ]

        			#current positions of particles
				set curpos [part $partnum print id pos]
			        set posx [lindex $curpos 1]
			        set posy [lindex $curpos 2]
			        set posz [lindex $curpos 3]

			        #reset x, y by considering periodic boundary 
				set posx [expr $posx + $dx]
				set posy [expr $posy + $dy]

				#update configuraiton
				part $partnum pos $posx $posy $posz 

			}
		puts "Done Periodic boundary check"

	   	}
}

