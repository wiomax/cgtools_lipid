# ::cgtools::espresso
#
# Routines for the hybrid MD-MC algorithm
# Author: Zun-Jing Wang  
# Nov 12 2009 Not Done Yet
#
namespace eval cgtools {
    namespace eval espresso {
	#################################################################################
	# proc hybrid_init
	proc hybrid_init { } {
	    #############################################################
	    # local variable: 	
	    #		  	f_topo
	    #                 	starttime
	    #                 	startj
	    #		  	startk
	    #                 	cutoff 
	    #                 	i
	    #                 	timingstart 
	    #                 	timingcurr
	    #                 	elapsedtime
	    #############################################################
		global topology 
	    	global jjjjjj
    		global kkkkkk

    		variable checkpointexists
		variable this
 
	     	set this [namespace current]
                mmsg::send $this "Feeding lipid parameters into ESPResSo running..."

    		# Attempt to read a checkpoint file
    		set checkpointexists [ ::cgtools::utils::readcheckpoint $cgtools::outputdir ]
 
    		# Set the starting time for this run ie override the value in checkpoint file
    		set starttime [clock seconds]

		#----------- Default Parameters set from System Params ----#
		if { $cgtools::warmup_temp == 0 } {
			set cgtools::warmup_temp [lindex $cgtools::systemtemp 0 ]
	   	}

		# ----------- Initialization ------------------ -----------#
		# Start a new computation, if no checkpoint
		if { !$checkpointexists } {
			# No checkpoint exists so we need to setup everything from scratch
    	        	set startj 0
    			set startk 0

    			# Setup the output directory by creating it and copying forcetables and overlapped potcoffs to it
    			::cgtools::utils::setup_outputdir  $cgtools::outputdir -paramsfile $cgtools::paramsfile \
				-tabdir $cgtools::tabledir -tabnames $cgtools::tablenames -coffdir $cgtools::overlapdir \
				-coffnames $cgtools::overlapnames -readpdbdir $cgtools::readpdbdir \
				-readpdbname $cgtools::readpdbname

    			# Set the box dimensions
    			setmd box_l [lindex $cgtools::setbox_l 0] [lindex $cgtools::setbox_l 1] [lindex $cgtools::setbox_l 2]
	
    			# Specify the bonded interactions
    			#puts "bonded_parms $cgtools::bonded_parms"
    			::cgtools::utils::set_bonded_interactions $cgtools::bonded_parms

    			# Specify any other non-bonded interactions
    			if { [ catch { ::cgtools::utils::set_nb_interactions $cgtools::nb_interactions } ] } {
				mmsg::send $this "no non-bonded interactions used"
    			}

    			set cutoff [setmd max_cut] 
    			puts "max_cut is $cutoff"

    			# Initialize variable moltypelists in the namespace ::cgtools::utils
    			::cgtools::utils::initmoltypeskey $cgtools::moltypelists 
    			#puts "$cgtools::moltypelists"

			#puts "$cgtools::ident"
			#puts "$cgtools::system_specs"
			#puts "$cgtools::setbox_l"
    			# Initialize topology
    			set topology [::cgtools::generation::generate_system $cgtools::system_specs $cgtools::setbox_l]
    			# Set the generated topology into the internals of espresso.
    			::cgtools::utils::set_topology $topology
    			#puts "$topology"
	    		
    			# See if there is any fixed molecules 
    			set cgtools::trappedmols [::cgtools::generation::get_trappedmols]
    			# Fix molecules if necessary
    			if { $cgtools::trappedmols != -1 } {
				::cgtools::utils::trap_mols $cgtools::trappedmols
    			}

			# set exclustions for the bonded particles 1: 2-bodyinteraction  2: 3-body interaction
    			part auto_exclusions 1

    			#Initialise Random Number Generator
    			::cgtools::utils::init_random $cgtools::nprocessors

    			# ----------- Integration Parameters before warmup -----------#
                	setmd periodic 1 1 1
    			setmd skin      $cgtools::verlet_skin

    			# Set the topology and molecule information
    			#----------------------------------------------------------#
    			#write topology file
    			set f_topo [open "$cgtools::outputdir/$cgtools::ident.top" w]
    			blockfile_write_topology $f_topo write topology   
    			close $f_topo

    			# Check if there are any extra vmdcommands and if not initialize a default
    			::cgtools::utils::initialize_vmd $cgtools::use_vmd $cgtools::outputdir $cgtools::ident \
				$topology -extracommands $cgtools::vmdcommands

    			#Perform the warm up integration
    
    			#----------------------------------------------------------#
    			# Warm up containing fixed particles 
    			setmd time_step $cgtools::warm_time_step
    			thermostat langevin $cgtools::warmup_temp $cgtools::langevin_gamma
    			mmsg::send $this "warming up at [setmd temp]"
    			::cgtools::utils::warmup  $cgtools::warmsteps $cgtools::warmtimes $topology -cfgs 10 \
				-outputdir $cgtools::outputdir

    			# If the membrane have any fixed particles, unfix them after warmup
    			set cgtools::userfixedparts [::cgtools::generation::get_userfixedparts ]
    			for {set i 0} { $i <  [setmd n_part] } {incr i} {
				if { [lsearch $cgtools::userfixedparts $i ] == -1 } {
	    				part [expr $i] fix 0 0 0
				}
    			}

    			# Set MD step, themostat before warm up without any fixed article
    			setmd time_step $cgtools::main_time_step
    			thermostat langevin  [lindex $cgtools::systemtemp 0] $cgtools::langevin_gamma

    			# Warm up without any fixed particle 
    			::mmsg::send $this "warming up without fixed particles at  [setmd temp]"
    			::cgtools::utils::warmup $cgtools::free_warmsteps $cgtools::free_warmtimes $topology \
				-startcap 1000 -outputdir $cgtools::outputdir
    
    			# ----------- Integration Parameters after warmup -----------#
    			# Set MD step, themostat after warm up
    			setmd time_step $cgtools::main_time_step
    			thermostat langevin  [lindex $cgtools::systemtemp 0] $cgtools::langevin_gamma

    			# Setup analysis
    			::cgtools::analysis::setup_analysis $cgtools::analysis_flags -outputdir  $cgtools::outputdir \
				 -g $cgtools::mgrid -str $cgtools::stray_cut_off
    
    			mmsg::send $this "starting integration: run $cgtools::int_n_times times $cgtools::int_steps steps"

    			# Reset the time to a starttime (usually zero) after warmup
    			setmd time $cgtools::startmdtime   
		}

		# Resume a computation, if exists checkpoint
		if { $checkpointexists } {
    			# A checkpoint exists so all we need to do is reset the moltypelists, topology and setup analysis again
    
    			::cgtools::utils::initmoltypeskey $cgtools::moltypelists 
    			::cgtools::utils::read_topology "$cgtools::outputdir/$cgtools::topofile"
			
			# set exclustions for the bonded particles 1: 2-bodyinteraction  2: 3-body interaction
    			part auto_exclusions 1


   			::cgtools::analysis::setup_analysis $cgtools::analysis_flags -outputdir  $cgtools::outputdir \
			 	-g $cgtools::mgrid -str $cgtools::stray_cut_off
    			
			set topology [analyze set] 
    			::cgtools::utils::initialize_vmd $cgtools::use_vmd $cgtools::outputdir $cgtools::ident \
				$topology

    			# Make sure that we start exactly from where the checkpoint was written
    			set startj [set jjjjjj]
    			set startk [expr [set kkkkkk] + 1]
		}

	    	set timingstart [clock clicks -milliseconds]
	    	set jjjjjj $startj

		#Main Integration                                          #
		#----------------------------------------------------------#
		thermostat langevin $cgtools::systemtemp $cgtools::langevin_gamma

		if { $cgtools::thermo == "DPD" } {
			thermostat off
    			set dpd_r_cut [setmd max_cut]
		        thermostat set dpd $cgtools::systemtemp $cgtools::dpd_gamma $dpd_r_cut
		        mmsg::send $this "DPD thermostat has been set"
		        mmsg::send $this "Thermostat is: [thermostat]"
		}
		if { $cgtools::npt == "on" } {
		        integrate set npt_isotropic $cgtools::p_ext $cgtools::piston_mass 1 1 0
		        mmsg::send $this "npt integrator has been set"
		        flush stdout
		        #-cubic_box
		        thermostat set npt_isotropic $cgtools::systemtemp  $cgtools::gamma_0  $cgtools::gamma_v
		}

		set timingstart [clock clicks -milliseconds]
		for {set kkkkkk $startk } { $kkkkkk <  $cgtools::int_n_times } { incr kkkkkk} {
		    	mmsg::send $this "run $kkkkkk at time=[setmd time]"

		    	# Call all of the analyze routines that we specified when setting up our analysis
		    	::cgtools::analysis::do_analysis

		    	# If kkkkkk is a multiple of analysis_write_frequency then write the analysis results to file
		    	if { [expr $kkkkkk + 1] % $cgtools::analysis_write_frequency ==0 } {
				::cgtools::analysis::print_averages
				#::cgtools::utils::update_force $rdfcglist $rdfaalist $tabledir $tablenames
		    	}

			# If kkkkkk is a multiple of write_frequency then write out a full particle configuration
			if { [expr $kkkkkk + 1] % $cgtools::write_frequency ==0 } {
				polyBlockWrite "$cgtools::outputdir/$cgtools::ident.[format %04d $jjjjjj].out" \
					{time box_l npt_p_diff } \
					{id pos type mass v f molecule} 
				mmsg::send $this "wrote file $cgtools::outputdir/$cgtools::ident.[format %04d $jjjjjj].out " 
				flush stdout

				if { $cgtools::use_vmd == "offline" } {
			    		::cgtools::utils::writecrd_charmm \
						"$cgtools::outputdir/$cgtools::ident.vmd[format %04d $jjjjjj].crd" $topology
    					::cgtools::utils::writepdb_charmm \
						"$cgtools::outputdir/$cgtools::ident.vmd[format %04d $jjjjjj].pdb" $topology
			    		#writepdbfoldtopo "$cgtools::outputdir/$cgtools::ident.vmd[format %04d $jjjjjj].pdb"  
				}

				incr jjjjjj

			    	# Write a checkpoint to allow restarting.  Overwrites previous checkpoint
			   	mmsg::send $this "setting checkpoint $kkkkkk [setmd time] $jjjjjj"    
				catch { exec rm -f $cgtools::outputdir/checkpoint.latest.chk }
			        checkpoint_set "$cgtools::outputdir/checkpoint.latest.out"
			        # Try to copy a checkpoint to the backup checkpoint folder
			        # Usefull if the program crashes while writing a checkpoint
			        if { [ catch { exec cp -f $cgtools::outputdir/checkpoint.latest.out \
					$cgtools::outputdir/checkpoint_bak/checkpoint.latest.out } ] } {
					mmsg::warn $this "warning: couldn't copy backup checkpoint"
			    	}

			}
			#end of if { [expr $kkkkkk + 1] % $cgtools::write_frequency ==0 }
			# MC: mode moves
			if { [expr $kkkkkk + 1] % $cgtools::mc_frequency ==0 } {
				set number_accept 0
				for {set i_mc 0 } { $i_mc <  $cgtools::mc_moves } { incr i_mc} {
					set acceptflag [mc_run]
					if { $acceptflag == 1 } {
						incr number_accept 
					}
				}
				::mmsg::send $this "mc moves accepted ratio is  [expr $number_accept/$cgtools::mc_moves]"



				#----------------------------------------------------------#
                        	# Warm up containing fixed particles 
                        	setmd time_step $cgtools::warm_time_step
                    		thermostat langevin $cgtools::warmup_temp $cgtools::langevin_gamma
                        	mmsg::send $this "warming up at [setmd temp]"
                        	::cgtools::utils::warmup  $cgtools::warmsteps $cgtools::warmtimes $topology -cfgs 10 \
                                	-outputdir $cgtools::outputdir

                       		# Warm up without any fixed particle 
                       		setmd time_step $cgtools::main_time_step
                       		thermostat langevin  [lindex $cgtools::systemtemp 0] $cgtools::langevin_gamma
                       		::mmsg::send $this "warming up without fixed particles at  [setmd temp]"
                       		::cgtools::utils::warmup $cgtools::free_warmsteps $cgtools::free_warmtimes $topology \
                               		 -startcap 1000 -outputdir $cgtools::outputdir

			}

			# MD: integrat equations of motion
                        integrate $cgtools::int_steps

			# Set the elapsed CPU time in computation, do not count that used for warm up
			set timingcurr [clock clicks -milliseconds]
			set elapsedtime [expr  $timingcurr - $timingstart]
			    ::mmsg::send $this "elapsed time: $elapsedtime"
		}
		#end of MD integration for {set kkkkkk $startk } { $kkkkkk <  $cgtools::int_n_times } { incr k}

		# Clean up the table, coefficent, readpdb files if they are in [pwd] or $outputdir
		#::cgtools::utils::cleanup_tabfiles $cgtools::tablenames $cgtools::tabledir $cgtools::outputdir
		#::cgtools::utils::cleanup_cofffiles $cgtools::overlapnames $cgtools::overlapdir $cgtools::outputdir
		#::cgtools::utils::cleanup_readpdbfiles $cgtools::readpdbname
		return
	}

	#################################################################################
	# proc mc_run
	# Perform Monte Carlo moves
	proc mc_run { } {
	    	variable this
	    	::mmsg::send $this "Starting MC moves"

	    	#periodic boundary conditions for (x, y) of all the lipids to have them in [0: box_l], not necessary??

	    	#random generate a mode vector, (qx, qy) where qx/qy=2*PI*n/box_l , n={0, 1, 2, ....} 
	    	set nmax $cgtools::mgrid
    	    	set boxx [lindex $cgtools::setbox_l 0]
    	    	set boxy [lindex $cgtools::setbox_l 1]
	     	if {$boxx != $boxy} {
			::mmsg::err $this "the side length of the box is set wrong: box_x != box_y"
	     	}
    	     	set nnow 0
	     	while {($nnow >= 8) || ($nnow == 0)} {
	    		set nnow_x [expr int($nmax*[t_random])]
    			set nnow_y [expr int($nmax*[t_random])]
	    		#choose mode vector with probobility of exp(aaa*1/q4), where aaa is proportional to kc from experiment
	    		#or choose n^2<=8
            		set nnow [expr $nnow_x*$nnow_x + $nnow_y*$nnow_y]
	     	}
	     	set nvector_now 0
	     	unset nvector_now
	     	lappend nvector_now $nnow_x
		lappend nvector_now $nnow_y
			
		set pos_list_old 0
		unset pos_list_old	

		#set initial_amplitude [expr $boxx/200.]
		set initial_amplitude 0.2 

		#change configurations of all the particles
		set pos_list_old [reset_z $nvector_now $initial_amplitude $boxx]

		#compute potential energy before a trial move
	        set epotold [expr [analyze energy total] - [analyze energy kinetic]]

		::mmsg::send $this "epotold: $epotold"
		#**************** Not sure if I will use  this part ************************
		#Biased MC
		#generate new configurations by moving center of the lipids compute Utrans, irelevant to tilts
		#generate k trial tilts from current center of lipids
		#for each of trial tilts, j={1, ...., k}, calculate the energy Utilt,j
		#Rosenbluth factor: W(n) = sumj exp(-Utilt,j/kBT) j=1,2...k
		#select one with a probability p(bn)=exp(-Utilt,n/kBT)/W(n)
		#for old configuration, generate k-1 trial tilts from old center of lipids
                #compute W(o) = exp(-Utilt,o/kBT) + sumj exp(-Utilt,j/kBT) j=1,2,...k-1
                #acceptance probability acc(o->n) = min(1, W(n)/W(o)exp{-(Utrans(n)-U(o)trans)/kBT}
		#***************************************************************************
		#compute potential energy after trial move
	        set epotnew [expr [analyze energy total] - [analyze energy kinetic]]
		::mmsg::send $this "epotnew: $epotnew"
		#Metropolis rool for accepting the trial move acc(o->n) = min(1, W(n)/W(o)exp()
	        set boltzmannweight [expr exp( -1. * ([set epotnew] - [set epotold]) ) ]
		if {$boltzmannweight > 1.0} {
			return_z $pos_list_old	
                        mmsg::send $this "The MC move is rejected"
			return 0

		} else {
                        mmsg::send $this "The MC move is accepted"
	    		return  1
		}
		
	}

	#################################################################################
	# Resetting z of particles according to x, y 
	proc reset_z { nvector initial_amplitude boxlength} {

		global topology 

		set nx [lindex $nvector 0]
		set ny [lindex $nvector 1]

		set pos_list_old 0
		unset pos_list_old

   		foreach mol $topology {
   			set moltype [lindex $mol 0]
		    	set typeinfo [::cgtools::utils::matchtype $moltype ]
  
			set partbondlists [lindex $typeinfo 2]
		        set partbondtypelists [lindex $typeinfo 3]

			set beadlists [lindex $partbondlists 0]
			set nbeads_mol [llength $beadlists]

		        set beadtypelists [lindex $partbondtypelists 0]

			set sum_mass 0.
			set sum_massx 0.
			set sum_massy 0.
			#set sum_massz 0.

			for { set b 0 } { $b < $nbeads_mol } {incr b } {
			        set partnum [lindex $mol [ expr $b + 1] ]
				set parttype [lindex $beadlists $b]
		        	set parttypeinfo [lindex $beadtypelists $parttype]
			        set partmass [lindex $parttypeinfo 2]

        			#current positions of particles
				set curpos [part $partnum print id pos]
				#save old configurations 
			        lappend pos_list_old $curpos

			        set posx [lindex $curpos 1]
			        set posy [lindex $curpos 2]
			        #set posz [lindex $curpos 3]	

			        #center of mass of the molecule
				set sum_massx [expr $sum_massx + $partmass * $posx]
				set sum_massy [expr $sum_massy + $partmass * $posy]
				#set sum_massz [expr $sum_massz + $partmass * $posz]
				set sum_mass [expr $sum_mass + $partmass]
				
			}
			
			set sum_massx [expr $sum_massx/$sum_mass]
			set sum_massy [expr $sum_massy/$sum_mass]
			#set sum_massz [expr $sum_massz/$sum_mass]

			#compute dz according to x,y
			set amplitude_max [expr $initial_amplitude/($nx*$nx+$ny*$ny)]
	    		set amplitude [expr $amplitude_max*([t_random] - 0.5)*2.0 ]
			set phase [expr 2*3.1415926*[t_random]]
			set qnow_x [expr 2*3.1415926/$boxlength*$nx]
			set qnow_y  [expr 2*3.1415926/$boxlength*$ny] 
			set dz [expr $amplitude*cos($qnow_x*$sum_massx+$qnow_y*$sum_massy+$phase)]

			for { set b 0 } { $b < $nbeads_mol } {incr b } {

			        set partnum [lindex $mol [ expr $b + 1] ]

        			#current positions of particles
				set curpos [part $partnum print id pos]
				#save old configurations 
			        lappend pos_list_old $curpos

			        #set posx [lindex $curpos 1]
			        #set posy [lindex $curpos 2]
			        set posz [lindex $curpos 3]	

			        #reset z of each particle of the molecule
				set posz [expr $posz + $dz]

				#update configuraiton from dz
				part $partnum pos $posx $posy $posz 

			}

	   	}
		return $pos_list_old

	}

	#################################################################################
	# Returning z of particles 
	proc return_z { pos_list_old } {

		global topology 

   		foreach mol $topology {
   			set moltype [lindex $mol 0]
		    	set typeinfo [::cgtools::utils::matchtype $moltype ]

			set partbondlists [lindex $typeinfo 2]
		    	set beadlists [lindex $partbondlists 0]
			set nbeads_mol [llength $beadlists]

			for { set b 0 } { $b < $nbeads_mol } {incr b } {
        			#current positions of particles
			        set partnum [lindex $mol [ expr $b + 1] ]

				set curpos [lindex $pos_list_old $b]
				set partnum_old [lindex $curpos 0]	
				if{ $partnum_old != $partnum }{
					mmsg::err [namespace $this] " wrong when reset_z"
				}
				
			        set posx [lindex $curpos 1]
			        set posy [lindex $curpos 2]
			        set posz [lindex $curpos 3]	

				#return to the old configuraiton
				part $partnum pos $posx $posy $posz 
			}
	   	}
	}
#end "namespace eval espresso"
    } 
#end "namespace eval cgtools"
}

