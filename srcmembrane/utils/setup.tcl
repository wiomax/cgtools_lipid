# cgtools::utils -- 
# Author: Zun-Jing Wang 
# Aug.-Sep. 2008 done
#
# This package provides a set of routines for basic bookkeeping tasks
# that typically arise in a membrane simulation #
# 

namespace eval ::cgtools::utils {
    namespace export setup_outputdir
    namespace export read_startfile
    namespace export readcheckpoint
    namespace export read_topology
    namespace export set_topology
    namespace export set_bonded_interactions
    namespace export set_nb_interactions
    namespace export init_random
    namespace export initialize_vmd
    namespace export cleanup_tabfiles
    namespace export cleanup_cofffiles
    namespace export cleanup_readpdbfiles
    namespace export uniqkey
    namespace export sleep
}

# ::cgtools::utils::setup_outputdir --
#
# This routine is designed to setup a directory for simulation
# output. It copies forcetables, overlapped coefficient files and the parameter file to the
# directory after creating it if necessary.  
#
# Arguments: outputdir The full pathname of the directory to be setup.
#                        At least the parent of this directory must
#                        exist
#
#            tabdir    The directory containing forcetables
#
#            tabnames  The list of all forcetables see
#                         set_tab_interactions
# 
#            coffdir    The directory containing coefficients of
#                          overlapped potential 
#
#            coffnames  The list of all coefficient file of 
#                          overlapped potential 
# 
#            paramsfile: The name of the parameter file 
#                        to be transferred to the directory
#
proc ::cgtools::utils::setup_outputdir { outputdir args } {

    # ---- Process Arguments ---------------------------# 
    # the ntabs option should not be necessary .. remove in a later
    # version
    set options {
	{paramsfile.arg  ""    "parameter file"}
	{tabdir.arg      ""    "forcetable directory"}
	{tabnames.arg    ""    "list of forcetable names"}
	{coffdir.arg     ""    "potcoff directory"}
	{coffnames.arg   ""    "list of potcoff names"}
	{readpdbdir.arg  ""    "readpdb directory"}
	{readpdbname.arg ""    "list of readpdb file name"}
    }
    set usage "Usage: setup_outputdir \[paramsfile:tabdir:tabnames:coffdir:coffnames:readpdbdir:readpdbname] outputdir "
    # Strip off optional arguments and store in params
    array set params [::cmdline::getoptions args $options $usage]
    
    
    ::mmsg::send [namespace current]  "setting up $outputdir for $params(paramsfile)"
    # If <outputdir> doesn't exist then create it
    set errcode [catch { exec mkdir $outputdir }]
    
    # If <outputdir/$params(tabdir)> doesn't exist then create it
    set errcode [catch { exec mkdir $outputdir/$params(tabdir) }]
    
    set ntabs [llength $params(tabnames)]
    for { set i 0 } { $i < $ntabs } { incr i } {
	set tablename [lindex $params(tabnames) $i ]
        # Copy forcetables to current directory
	#set errcode [ catch { exec cp $params(tabdir)/$tablename [pwd]/ } ]	    
	#if { $errcode } {
	#    ::mmsg::warn [namespace current]  "couldn't transfer forcetable $params(tabdir)/$tablename to [pwd]"
	#} else {
	#    ::mmsg::send [namespace current]  "copied $params(tabdir)/$tablename to [pwd] "
	#}

        # Copy forcetables to output directory
	set errcode [ catch { exec cp $params(tabdir)/$tablename $outputdir/$params(tabdir)/ } ]
	if { $errcode } {
	    ::mmsg::warn [namespace current]  "couldn't transfer forcetable $params(tabdir)/$tablename to $outputdir/$params(tabdir)"
	} else {
	    ::mmsg::send [namespace current]  "copied $params(tabdir)/$tablename to $outputdir/$params(tabdir) "
	}
    }

    # If <outputdir/$params(coffdir)> doesn't exist then create it
    set errcode [catch { exec mkdir $outputdir/$params(coffdir) }]

    set ncoffs [llength $params(coffnames)]
    for { set i 0 } { $i < $ncoffs } { incr i } {
        set coffname [lindex $params(coffnames) $i ]
        # Copy coefficient files to current directory
        #set errcode [ catch { exec cp $params(coffdir)/$coffname [pwd]/ } ]
        #if { $errcode } {
        #    ::mmsg::warn [namespace current]  "couldn't transfer overlapped coefficient $params(coffdir)/$coffname to [pwd]"
        #} else {
        #    ::mmsg::send [namespace current]  "copied $params(coffdir)/$coffname to [pwd] "
        #}

        # Copy coefficient files to output directory
        set errcode [ catch { exec cp $params(coffdir)/$coffname $outputdir/$params(coffdir) } ]
        if { $errcode } {
            ::mmsg::warn [namespace current]  "couldn't transfer overlapped coefficient $params(coffdir)/$coffname to $outputdir/$params(coffdir)"
        } else {
            ::mmsg::send [namespace current]  "copied $params(coffdir)/$coffname to $outputdir/$params(coffdir) "
        }
    }

    
    #Copy the paramsfile to the outputdir
    catch { exec cp $params(paramsfile) $outputdir }
    
    #Copy the readfile to the outputdir and current dir
    catch { exec cp $params(readpdbdir)/$params(readpdbname) $outputdir }
    #catch { exec cp $params(readpdbdir)/$params(readpdbname) [pwd]/ }
    
    # Construct a directory for checkpoint backups inside outputdir
    catch { exec mkdir $outputdir/checkpoint_bak }    
}



# ::cgtools::utils::read_startfile --
#
# Read in configuration information from an existing config file
# 
#
proc ::cgtools::utils::read_startfile { file } {
    ::mmsg::send [namespace current]   "reading config: $file" nonewline 
    flush stdout
    set mark .gz
    set fname $file
    if { [ regexp $mark $fname ] } {
	set f [open "|gzip -cd $file" "r"]
    } else {
	set f [open "$file" "r"]
    }

    while { [ blockfile $f read auto ] != "eof" } {
	::mmsg::send [namespace current]  "." nonewline
	flush stdout
    }
    
    close $f
    ::mmsg::send [namespace current]  " done" 
}


# ::cgtools::utils::readcheckpoint --
#
#   A wrapper routine for reading checkpoints that checks for success
#
proc ::cgtools::utils::readcheckpoint { checkpointdir } {
    set result [ catch { set f [ open "$checkpointdir/checkpoint.latest.out" "r" ]  } ]
    if { $result } {
	::mmsg::warn [namespace current]  "no checkpoint named $checkpointdir/checkpoint.latest.out"
	return 0
    }
    #::mmsg::send [namespace current] "reading Checkpoint $checkpointdir/checkpoint.latest.out"
    checkpoint_read "$checkpointdir/checkpoint.latest.chk"
    #set testresult [ catch { checkpoint_read "$checkpointdir/checkpoint.latest.chk" } ]
    #while { $testresult } {		
    #	::cgtools::utils::sleep 50000
    # 	set testresult [ catch { checkpoint_read "$checkpointdir/checkpoint.latest.chk" } ]
    #}
    close $f
    return 1
}

proc ::cgtools::utils::uniqkey { } {
     set key   [ expr { pow(2,31) + [ clock clicks ] } ]
     set key   [ string range $key end-8 end-3 ]
     set key   [ clock seconds ]$key
     return $key
 }

proc ::cgtools::utils::sleep { ms } {
     set uniq [ uniqkey ]
     set ::__sleep__tmp__$uniq 0
     after $ms set ::__sleep__tmp__$uniq 1
     vwait ::__sleep__tmp__$uniq
     unset ::__sleep__tmp__$uniq
 }

# ::cgtools::utils::read_topology --
#
#   Uses the blockfile read command to read topology information
#   contained in a file "$dir/file" into espresso.  The analyze
#   "topo_part_sync" command is then used to synchronise the molecule
#   information contained in topology and particle molecule id's
#
proc ::cgtools::utils::read_topology { file } {
    set f [open $file r]
    blockfile $f read topology
    analyze set "topo_part_sync"
    close $f
}

# ::cgtools::utils::set_topology --
#
#  Takes a topology list "topo" and sets this into espresso. The
#   analyze "topo_part_sync" command is then used to synchronise the
#   molecule information contained in topology and particle molecule
#   id's
#
proc ::cgtools::utils::set_topology { topo } {
    eval analyze set $topo
    analyze set "topo_part_sync"
}


# ::cgtools::utils::set_bonded_interactions --
#
#  This routine uses the inter command of espresso to set all bonded
#  interactions.
#
proc ::cgtools::utils::set_bonded_interactions { bonded_parms } {
    foreach bondtype $bonded_parms {
        set bondsetsuc [catch {eval [concat inter $bondtype] } ]
        set iloop 0
	while { $bondsetsuc } {
	   mmsg::err [namespace current] "couldn't set interaction: [concat [lindex $bondtype 0]], try again"
	   ::cgtools::utils::sleep 5000
           set bondsetsuc [catch {eval [concat inter $bondtype] } ]
           set iloop [expr $iloop + 1]
           if {$iloop > 10} break
	} 
	mmsg::send [namespace current] "set interaction: $bondtype "
    }
    return

    #foreach bondtype $bonded_parms {
    #	if { [catch {eval [concat inter $bondtype] } ] } {
    #	    mmsg::err [namespace current] "couldn't set interaction: [concat [lindex $bondtype 0]]"
    #	} else {	
    #	    mmsg::send [namespace current] "set interaction: $bondtype "
    #	}
    #}
    #return
}

# ::cgtools::utils::set_nb_interactions --
# 
# Set all the non-bonded interactions apart from tabulated ones, eg
# for constraints for instance.
#
proc ::cgtools::utils::set_nb_interactions { interactionlist } {
    foreach intertype $interactionlist {
        set nobondsetsuc [catch { eval [concat inter  $intertype ] } ] 
        set iloop 0
	while { $nobondsetsuc } {
	    mmsg::err [namespace current] "could not set interaction: $intertype, try again"
            ::cgtools::utils::sleep 5000
            set nobondsetsuc [catch { eval [concat inter  $intertype ] } ]
            set iloop [expr $iloop + 1]
            if {$iloop > 10} break
	}
	mmsg::send [namespace current] "set interaction: $intertype "
    }
    return

    #foreach intertype $interactionlist {
    #	if { [catch { eval [concat inter  $intertype ] } ] } {
    #	    mmsg::err [namespace current] "could not set interaction: $intertype"
    #	}
    #	mmsg::send [namespace current] "set interaction: $intertype "
    #}
    #return
}

# ::cgtools::utils::init_random --
#
# Initialize the random number generators on each processor 
#
proc ::cgtools::utils::init_random { n_procs } { 
    
    set c [expr [clock seconds] - 1068130800]    
    #    set c  1068130800    
    
    ::mmsg::send [namespace current]  "Setting the random seed to clock in seconds: $c "
    set seedbase $c
    for { set i 1 } { $i < $n_procs } { incr i } {
	lappend seedbase [expr $c + 2304*$i]
    }
    eval t_random seed $seedbase
    
    flush stdout
}
    
# ::cgtools::utils::initialize_vmd --
#
# Depending on the value of <flag> initialize vmd to one of two possible states:
#
#  interactive: VMD is started and a connection to espresso
#               established for immediate viewing of the current
#               espresso process. With some luck this might even work
#               sometimes!!!  If VMD doesn't get a proper connection
#               to espresso then it will crash.
#
#  offline: Just constructs the appropriate psf and
#               vmd_animation.script files and writes them to the
#               output directory so that pdb files generated with
#               writepdb can be viewed with vmd -e
#               vmd_animation.script
#
#  default: Any value other than those above for flag will just result
#               in vmd not being initialized.
#
proc ::cgtools::utils::initialize_vmd { flag outputdir ident topology args } {
    # ---- Process Arguments ---------------------------# 
    set options {
	{extracommands.arg       ""     "additional stuff to be written to vmd_animation script"}
    }
    set usage "Usage: initialize_vmd flag $outputdir $ident $topology $args \[extracommands]"
    # Strip off optional arguments and store in params
    array set params [::cmdline::getoptions args $options $usage]

    # set some defaults
    set filename "vmd"  
    set wait "0"
    set start "1" 
    
    ::mmsg::send [namespace current]  "initializing vmd to : " nonewline
    
    switch -regexp $flag  {
	"interactive" {
	    ::mmsg::send [namespace current]  "interactive" 
	    writepsf_xplor "$filename.psf" $topology
	    #writepsf "$filename.psf"
	    writecrd_charmm "$filename.crd" $topology
	    writepdb_charmm "$filename.pdb" $topology
	    #writepdb "$filename.pdb"
	    for {set port 10000} { $port < 65000 } { incr port } {
		catch {imd connect $port} res
		if {$res == ""} break
	    }
	    set HOSTNAME [exec hostname]
	    set vmdout_file [open "vmd_start.script" "w"]
	    puts $vmdout_file "mol load psf $filename.psf pdb $filename.pdb"
	    puts $vmdout_file "rotate stop"
	    puts $vmdout_file "mol modstyle 0 0 Licorice"
	    # 1.800000 0.300000 8.000000 6.000000"
	    puts $vmdout_file "mol modcolor 0 0 Name"
	    puts $vmdout_file "imd connect $HOSTNAME $port"
	    puts $vmdout_file "imd transfer 1"
	    puts $vmdout_file "imd keep 1"
	    close $vmdout_file
	    if { $start == 0 } {
		::mmsg::send [namespace current]  "Start VMD in the same directory on the machine you with :"
		::mmsg::send [namespace current]  "vmd -e vmd_start.script &"
		imd listen $wait
	    } else {
		exec vmd -e vmd_start.script &
	    }
	}
	"offline" {
	    ::mmsg::send [namespace current]  "offline"
	    variable firstconfignum
	    
	    writepsf_xplor "$outputdir/$ident.vmd.psf" $topology
	    #writepsf "$outputdir/$ident.vmd.psf"
	    
	    set vmd_file [open "$outputdir/vmd_animation.script" "w"]
	    puts $vmd_file "loadseries $ident.vmd 1 $firstconfignum"
	    puts $vmd_file "rotate stop"
	    puts $vmd_file "mol modstyle 0 0 Licorice"
	    puts $vmd_file "mol modcolor 0 0 Name"
	    puts $vmd_file "logfile vmd.log"
	    puts $vmd_file "logfile off"
	    foreach command $params(extracommands) {
		puts $vmd_file $command
	    }
	    close $vmd_file
	    
	    # And for the warmup too
	    writepsf_xplor "$outputdir/warm.vmd.psf"  $topology
	    #writepsf "$outputdir/warm.vmd.psf" 
	    set vmd_file [open "$outputdir/warm_animation.script" "w"]
	    puts $vmd_file "loadseries warm.vmd 1 0"
	    puts $vmd_file "rotate stop"
	    puts $vmd_file "mol modstyle 0 0 Licorice"
	    puts $vmd_file "mol modcolor 0 0 Name"
	    puts $vmd_file "logfile vmd.log"
	    puts $vmd_file "logfile off"
	    foreach command $params(extracommands) {
		puts $vmd_file $command
	    }
	    close $vmd_file
	}


	"default" { 
	    ::mmsg::send [namespace current]  "default"
	    #Do nothing 
	}
    }
}

# ::cgtools::utils::set_std_topology --
#    
# This routine is used when a starting config file is specified but no
# topology information is given.  In that case as a last resort we
# simply assume a homogenous bilayer ie standard topology.
proc ::cgtools::utils::set_std_topology { n_parts beads_per_lipid } {
    ::mmsg::send [namespace current] "assuming flat bilayer with one lipid type"
    set n_lipids_monolayer [expr $n_parts/($beads_per_lipid*2.0)]

    for {set i 0} { $i < $n_lipids_monolayer } {incr i} {
	set mol1 0
	set mol2 1
	for { set b 0 } { $b < $beads_per_lipid } {incr b } {
	    set partnum [expr $i*$beads_per_lipid*2 + $b ]
	    lappend mol1 $partnum
	}
	for { set b $beads_per_lipid } { $b < [expr $beads_per_lipid*2] } { incr b } {
	    set partnum [expr $i*$beads_per_lipid*2 + $b ]
	    lappend mol2 $partnum
	}
	lappend topo $mol1
	lappend topo $mol2
    }
    eval analyze set $topo
    analyze set "topo_part_sync"
}

    
# No longer used ?
proc ::cgtools::utils::prepare_vmd_series { outdir name topology } {

    writepsf_xplor "$outdir/$name.vmd.psf"  $topology
    #writepsf "$outdir/$name.vmd.psf" 
    
    set vmd_file [open "$outdir/vmd_animation.script" "w"]
    puts $vmd_file "loadseries $outdir/$name.vmd 1"
    puts $vmd_file "rotate stop"
    puts $vmd_file "mol modstyle 0 0 Licorice"
    puts $vmd_file "mol modcolor 0 0 Name"
    puts $vmd_file "logfile $outdir/vmd.log"
    puts $vmd_file "logfile off"
    close $vmd_file
}

proc ::cgtools::utils::probe_nparts { f } {
    blockfile $f read variable
    blockfile $f read variable
    
    blockfile $f read interactions
    blockfile $f read particles  
    blockfile $f read bonds
    
    set n_lipids "[expr [ setmd n_part]/3]"
    return $n_lipids
}


proc ::cgtools::utils::cleanup_tabfiles { tablenames tabledir outputdir } {

    
    # If <outputdir/tabdir> doesn't exist then create it
    catch { exec mkdir $outputdir/$tabledir  }

    # Copy forcetables to output/tabdir directory
    for { set i 0 } { $i < [llength $tablenames] } {incr i} {
	set errcode [catch { exec cp [lindex $tablenames $i] $outputdir/$tabledir/ } ]
	if { $errcode } {
	    ::mmsg::err [namespace current]  "Error copying table to $outputdir/$tabledir : $errcode"
	}
    }
    #clean up files
    for { set i 0 } { $i < [llength $tablenames] } {incr i} {
	catch { exec rm [lindex $tablenames $i] }
	catch { exec rm $outputdir/[lindex $tablenames $i] }
    }
}

proc ::cgtools::utils::cleanup_cofffiles { coffnames coffdir outputdir } {

    # If <outputdir/coffdir> doesn't exist then create it
    catch { exec mkdir $outputdir/$coffdir  }

    # Copy coefficient files to output/coffdir directory
    for { set i 0 } { $i < [llength $coffnames] } {incr i} {
        set errcode [catch { exec cp -f [lindex $coffnames $i] $outputdir/$coffdir/ } ]
        if { $errcode } {
            ::mmsg::err [namespace current]  "Error copying overlapped coefficient to $outputdir/$coffdir : $errcode"
        }
    }
    #clean up files
    for { set i 0 } { $i < [llength $coffnames] } {incr i} {
        catch { exec rm [lindex $coffnames $i] }
        catch { exec rm $outputdir/[lindex $coffnames $i] }
    }
}

proc ::cgtools::utils::cleanup_readpdbfiles { readpdbname } {
        catch { exec rm $readpdbname }
}
