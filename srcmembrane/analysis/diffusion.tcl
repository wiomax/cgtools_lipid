# ::cgtools::analysis::analyze_diffusion --
#
# Author: Zun-Jing Wang 
# Oct 14 2009

# When print_averages is called the values stored in time_vs_msd are
# printed.  Generally for a diffusion analysis we would only call
# print_averages right at the end so that the entire simulation is
# used as an average.  For this reason this analysis is best done on
# stored configuration files rather than during simulation.
#

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::diffusion {  
    variable i_conf
    variable moltypenow
    variable ntimesteps
    variable numberconfs  
    variable MSD2dlist
    namespace export printav_diffusion
    namespace export setup_diffusion
    namespace export analyze_diffusion
    namespace export resetav_diffusion
}
proc ::cgtools::analysis::diffusion::resetav_diffusion {} {

    variable i_conf   
    set i_conf 0   
    analyze remove 

}

proc ::cgtools::analysis::diffusion::printav_diffusion { } {
    variable i_conf 
    variable moltypenow
    variable ntimesteps
    variable numberconfs  
    variable MSD2dlist
    global ::cgtools::analysis::outputdir

   puts "moltypenow : $moltypenow"
   puts "ntimesteps: $ntimesteps"
   puts "numberconfs : $numberconfs"
   puts "i_conf: $i_conf"
   puts "numberconfs: $numberconfs"

    if {  $i_conf == [set numberconfs] } {
    	::mmsg::send [namespace current] "analyzing the MSD2d"
    	set MSD2dlist [ analyze MSD2d $moltypenow $ntimesteps $numberconfs ]
    	puts "MSD2dlist: $MSD2dlist"
    	puts " [lindex [lindex $MSD2dlist 0] 0]"
    	puts " [lindex [lindex $MSD2dlist 0] 1]"
    	puts " [lindex [lindex $MSD2dlist 0] 2]"


        ::mmsg::debug [namespace current] "opening $outputdir/time_vs_msd2d"
    	set f_msd2d [open "$outputdir/time_vs_msd2d" w]
	puts $f_msd2d "\# [lindex [lindex $MSD2dlist 0] 0]"
	puts $f_msd2d "\# [lindex [lindex $MSD2dlist 0] 2]"
        set msd2d [lindex [lindex $MSD2dlist 0] 1]

        #set tlist ""
        #set msdlist ""
        foreach value $msd2d {
		set tnow [lindex $value 0]
                set msdnow [lindex $value 1]
	        puts $f_msd2d "$tnow $msdnow"
        }
    	close $f_msd2d


    } else {
	::mmsg::warn [namespace current] "can't print msd2d"
	flush stdout
    }
}

proc ::cgtools::analysis::diffusion::setup_diffusion { args } {
    variable i_conf
    variable moltypenow
    variable ntimesteps
    variable numberconfs  
    variable MSD2dlist
    global ::cgtools::analysis::outputdir

    set options {
	{typem.arg "0"  "type of molecule in computing MSD2d"}
	{ntimesteps.arg "1" "delta_t in computing MSD2d"}
	{numberofconf.arg "10" "number of configurations in computing MSD2d"}
    }
    set usage "Usage: ::cgtools::analysis::diffusion::setup_diffusion \[type_m:n_time_steps:number_of_conf]"
    array set params [::cmdline::getoptions args $options $usage]

    set i_conf 0   
    set moltypenow $params(typem)
    set ntimesteps $params(ntimesteps)
    set numberconfs $params(numberofconf)
    set MSD2dlist ""
    analyze remove 
}

proc ::cgtools::analysis::diffusion::analyze_diffusion { } {
    variable i_conf
    variable moltypenow
    variable numberconfs  
    variable MSD2dlist

    ::mmsg::send [namespace current] "save configurations"
    analyze append
 
    incr i_conf
}
