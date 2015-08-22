# Obtained from MBTOOLS, modified by Zun-Jing Wang

namespace eval ::cgtools::utils {
    namespace export warmup  
    namespace export nptwarmup
    namespace export mcwarmup
}

#----------------------------------------------------------#
# ::cgtools::utils::warmup--
#
# basic warmup routine
#
# Can be used in a variety of ways depending on the parameters and
# options that are set. 
#
# Options:
#
# mindist: Set this option to a minimum particle distance requirement.
# Routine will use the mindist function of expresso (very very slow )
# to determine if this criterion is satisfied and if it is the warmup
# will be terminated.
#
# imd: Use imd to watch the warmup on the fly 
#
# capincr: This is normally calculated from capgoal and startcap.  If
# it is set then it will override the calculated value
#
proc ::cgtools::utils::warmup { steps times topology args } {
    variable warmcfg
    ::mmsg::send [namespace current] "warming up $times times $steps timesteps "
    set options {
	{mindist.arg     0   minimum distance criterion to terminate warmup }
	{cfgs.arg -1   write out a record of the warmup at given intervals}
	{outputdir.arg  "./" place to write warmup configs}
	{vmdflag.arg "offline" vmd settings}
	{startcap.arg     5    initial forcecap}	
	{capgoal.arg      1000 forcecap goal}
    }
    set usage "Usage: warmup steps times \[mindist:imd:startcap:capincr:capgoal:]"
    array set params [::cmdline::getoptions args $options $usage]
    
    #Make a sanity check
    if { $steps == 0 || $times == 0 } {
	::mmsg::warn [namespace current] "warmup steps are zero"
	return
    }

    # Work out the cap increment if it is not set
    set capincr [expr ($params(capgoal) - $params(startcap))/($times*1.0)]


    # Set the initial forcecap
    set cap $params(startcap)


    for { set i 0 } { $i < $times } { incr i } {
	# Check the mindist criterion if necessary
	if { $params(mindist) } {
	    set act_min_dist [analyze mindist]
	    if { $act_min_dist < $params(mindist) } {
		break
	    }
	}

	# Write out configuration files and pdb files
	if { $i%$params(cfgs)==0 && ($params(cfgs) > 0 ) } {
	    #polyBlockWrite "$params(outputdir)/warm.[format %04d $warmcfg].out" {time box_l npt_p_diff } {id pos type v f molecule} 
	    #mmsg::send [namespace current] "wrote file $params(outputdir)/warm.[format %04d $warmcfg].out " 

	    flush stdout

	    if { $params(vmdflag) == "offline" } {
		writepdb_charmm "$params(outputdir)/warm.vmd[format %04d $warmcfg].pdb" $topology 
		#writepdbfoldtopo "$params(outputdir)/warm.vmd[format %04d $warmcfg].pdb"  
	    }
	    incr warmcfg
	}

	# Set the new forcecap into espresso and integrate
	# catch tabulated force cap error in case tabulated option is not turned on
	catch {inter tabforcecap $cap}
	catch {inter ljforcecap $cap}
	integrate $steps
	set cap [expr $cap + $capincr ]
	::mmsg::send [namespace current]  "run $i of $times at time=[setmd time] (cap=$cap) " 

	flush stdout
	
    }
    
    # Turn off all forcecapping
    ::mmsg::send [namespace current] "uncapping forces"
    catch {inter tabforcecap 0}
    catch {inter ljforcecap 0}
 }


#----------------------------------------------------------#
# ::cgtools::utils::nptwarmup--
#
# Routine for decreasing the box mass and friction and allowing box to reach
#
# Can be used in a variety of ways depending on the parameters and
# options that are set. 
#
#
# Note: Assumes a 2d npt ie constant tension
# Options:
#
# 
#
# 
proc ::cgtools::utils::nptwarmup { steps times iparms fparms args } {
    variable warmcfg
    ::mmsg::send [namespace current] "nptwarmup $times times $steps timesteps "
    set options {
	{cfgs.arg -1   write out a record of the warmup at given intervals}
	{outputdir.arg  "./" place to write warmup configs}
	{vmdflag.arg "offline" vmd settings}
	{onedim.arg -1 "only does npt in 1 dimension"}
    }
    set usage "Usage: warmup steps times \[mindist:imd:startcap:capincr:capgoal:]"
    array set params [::cmdline::getoptions args $options $usage]
    
    #Make a sanity check
    if { $steps == 0 || $times == 0 } {
	::mmsg::warn [namespace current] "npt warmup steps are zero"
	return
    }
    
    for { set p 0 } { $p < [llength  $iparms] } { incr p } {
	if {[lindex $iparms $p] == 0} {
	    if  {[lindex $fparms $p] == 0} {
		lappend parmincs 1.0
	    } else {
		mmsg::warn [namespace current] "cannot gradually increment an npt parameter starting from 0"
	    }
	} else {
	    lappend parmincs [expr pow([lindex $fparms $p ] / [lindex $iparms $p],1.0/(1.0*$times-1.0))]
	}
    }
    
    set systemtemp [setmd temperature ]
    
    for { set i 0 } { $i < $times } { incr i } {
	
	# Write out configuration files and pdb files
	if { $i%$params(cfgs)==0 && ($params(cfgs) > 0 ) } {
	    polyBlockWrite "$params(outputdir)/warm.[format %04d $warmcfg].out" {time box_l npt_p_diff } {id pos type v f molecule} 
	    mmsg::send [namespace current] "wrote file $params(outputdir)/warm.[format %04d $warmcfg].out " 
	    
	    flush stdout
	    
	    if { $params(vmdflag) == "offline" } {
		writepdb_charmm "$params(outputdir)/warm.vmd[format %04d $warmcfg].pdb"  $topology
		#writepdbfoldtopo "$params(outputdir)/warm.vmd[format %04d $warmcfg].pdb"  
	    }
	    incr warmcfg
	}
	
	set p_ext [expr [lindex $iparms 0] * pow([lindex $parmincs 0],$i) ]
	set piston_mass [expr [lindex $iparms 1] * pow([lindex $parmincs 1],$i) ]	
	set gamma_0 [expr [lindex $iparms 2] * pow([lindex $parmincs 2],$i)]
	set gamma_v [expr [lindex $iparms 3] * pow([lindex $parmincs 3],$i)]
	
	if {$params(onedim) == 0} {
	    integrate set npt_isotropic $p_ext $piston_mass 1 0 0 
	} elseif {$params(onedim) == 1} {
	    integrate set npt_isotropic $p_ext $piston_mass 0 1 0 
	} elseif {$params(onedim) == 2} {
	    integrate set npt_isotropic $p_ext $piston_mass 0 0 1 
	} else {
	    integrate set npt_isotropic $p_ext $piston_mass 1 1 0 
	}
	thermostat set npt_isotropic $systemtemp $gamma_0 $gamma_v 
	mmsg::send [namespace current] "npt integrator has been set $p_ext $piston_mass $systemtemp $gamma_0 $gamma_v "
	
	integrate $steps
	
	::mmsg::send [namespace current]  "run $i of $times at time=[setmd time] " 
	
	flush stdout
	
    }
    
}

#----------------------------------------------------------#
# ::cgtools::utils::warmup--
#
# warmup routine after MC moves
#
# Options:
#
# mindist: Set this option to a minimum particle distance requirement.
# Routine will use the mindist function of expresso (very very slow )
# to determine if this criterion is satisfied and if it is the warmup
# will be terminated.
#
# imd: Use imd to watch the warmup on the fly 
#
# capincr: This is normally calculated from capgoal and startcap.  If
# it is set then it will override the calculated value
#
proc ::cgtools::utils::mcwarmup { steps times topology args } {
    variable warmcfg
    ::mmsg::send [namespace current] "warming up $times times $steps timesteps "
    set options {
	{mindist.arg     0   minimum distance criterion to terminate warmup }
	{cfgs.arg -1   write out a record of the warmup at given intervals}
	{outputdir.arg  "./" place to write warmup configs}
	{vmdflag.arg "offline" vmd settings}
	{startcap.arg     5    initial forcecap}	
	{capgoal.arg      1000 forcecap goal}
    }
    set usage "Usage: warmup steps times \[mindist:imd:startcap:capincr:capgoal:]"
    array set params [::cmdline::getoptions args $options $usage]
    
    #Make a sanity check
    if { $steps == 0 || $times == 0 } {
	::mmsg::warn [namespace current] "warmup steps are zero"
	return
    }

    # Work out the cap increment if it is not set
    set capincr [expr ($params(capgoal) - $params(startcap))/($times*1.0)]


    # Set the initial forcecap
    set cap $params(startcap)


    for { set i 0 } { $i < $times } { incr i } {
	# Check the mindist criterion if necessary
	if { $params(mindist) } {
	    set act_min_dist [analyze mindist]
	    if { $act_min_dist < $params(mindist) } {
		break
	    }
	}

	# Write out configuration files and pdb files
	if { $i%$params(cfgs)==0 && ($params(cfgs) > 0 ) } {
	    polyBlockWrite "$params(outputdir)/warm.[format %04d $warmcfg].out" {time box_l npt_p_diff } {id pos type v f molecule} 
	    mmsg::send [namespace current] "wrote file $params(outputdir)/warm.[format %04d $warmcfg].out " 

	    flush stdout

	    if { $params(vmdflag) == "offline" } {
		writepdb_charmm "$params(outputdir)/warm.vmd[format %04d $warmcfg].pdb" $topology 
		#writepdbfoldtopo "$params(outputdir)/warm.vmd[format %04d $warmcfg].pdb"  
	    }
	    incr warmcfg
	}

	# Set the new forcecap into espresso and integrate
	# catch tabulated force cap error in case tabulated option is not turned on
	catch {inter tabforcecap $cap}
	catch {inter ljforcecap $cap}
	integrate $steps
	set cap [expr $cap + $capincr ]
	::mmsg::send [namespace current]  "run $i of $times at time=[setmd time] (cap=$cap) " 

	flush stdout
	
    }
    
    # Turn off all forcecapping
    ::mmsg::send [namespace current] "uncapping forces"
    catch {inter tabforcecap 0}
    catch {inter ljforcecap 0}
 }
