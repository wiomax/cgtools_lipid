# cgtools::utils -- 
#
# This package provides a set of routines for basic bookkeeping tasks
# that typically arise in a membrane simulation #
# Author: Zun-Jing Wang 
# Sep. 23 2008
# 

package require ::mmsg
package provide ::cgtools::utils 1.0.0

namespace eval ::cgtools::utils {
    variable moltypeskey

    variable verbosity 0
    variable firstconfignum 0

    # The identity number for written out warmup configurations
    variable warmcfg 0
    variable maxnbinteractiontype 0
}


proc ::cgtools::utils::initmoltypeskey { moltypelists } {
    variable moltypeskey
    set moltypeskey $moltypelists
}
# ::cgtools::utils::matchtype-- 
# Search the molecule type list for a type number and return the type specification
proc ::cgtools::utils::matchtype { mol } {
    variable moltypeskey

    foreach key $moltypeskey {
        if { [lindex $key 0] == [lindex $mol 0] } {
            return $key
        }
    }
    mmsg::err [namespace current] "could not find a matching key to moltype [lindex $mol 0]"
}

# Read histogram file for a given temperature. Used by Parallel tempering
proc ::cgtools::utils::read_histogram {filename} {
        array set histogram ""
        set f [open $filename r]
        set data [read $f]
        set data [split $data "\n"]
        close $f
        foreach line $data {
                if {[lindex $line 0] != "\#" && [llength $line] > 0} {
                        set histogram([lindex $line 0]) [lindex $line 1]
                }
        }
        return [array get histogram]
}


# Write histogram of energies to file at a given temperature. Used by Parallel tempering
proc ::cgtools::utils::write_histogram {filename current_energy} {
	array set histogram [read_histogram $filename]
	set f [open $filename w]
	# Add $current_energy to histogram variable.
	# The bins will have a precision of 10.
        set value [expr round($current_energy/10.)*10]
        if { [catch {incr histogram($value)} ] } {
                set histogram($value) 1
        }

        puts $f "\# Histogram of energies"
        puts $f "\# Potential energy \t Number of hits"
        foreach idx [array names histogram] {
                puts $f "$idx \t[set histogram($idx)]"
        }
        close $f
}

# Append observables to a file
proc ::cgtools::utils::append_obs {filename current_energy {id -1}} {
         set f [open $filename a]
         puts $f "[setmd time] " nonewline
         puts $f " \t$current_energy" nonewline
         if { $id != -1 } {
                 puts $f " \t$id"
         } else { 
                 puts $f ""
         }
         close $f
}

source [file join [file dirname [info script]] warmup.tcl]
source [file join [file dirname [info script]] topo.tcl]
source [file join [file dirname [info script]] math.tcl]
source [file join [file dirname [info script]] misc.tcl]
source [file join [file dirname [info script]] setup.tcl]
source [file join [file dirname [info script]] readpdbcrd.tcl]
source [file join [file dirname [info script]] writepdbcrdpsf.tcl]
source [file join [file dirname [info script]] update_force.tcl]
source [file join [file dirname [info script]] computeparameter.tcl]

