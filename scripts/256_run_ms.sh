#!/bin/bash

FLORIDA=true
TEXAS=false
ALABAMA=false
# generate 1000 samples of size number of inds * 2 with:
# pop of size 48 (ALFW) = 96
# pop of size 47 (ALST) = 94
# pop of size 41 (FLCC) = 82
# pop of size 47 (FLFW) = 94
# pop of size 48 (LAFW) = 96
# pop of size 41 (TXCC) = 82
# pop of size 31 (TXFW) = 62


if [ "$FLORIDA" = true ]; then
	# FLORIDA
	# nsam = 176 (82 + 94)
	# best model = asym_mig
	#   Split into two populations, with different migration rates.

    # 	nu1: Size of population 1 after split.
    # 	nu2: Size of population 2 after split.
    # 	T: Time in the past of split (in units of 2*Na generations) 
    # 	m12: Migration from pop 2 to pop 1 (2*Na*m12)
    # 	m21: Migration from pop 1 to pop 2
    # params: nu1=0.1317,nu2=8.4225,m12=0.7777,m21=0.0561,T=0.1441


	# 1000 reps of 176 samples, with 82 of those from pop 1 and 94 of those from pop 2
	# -t 1 means that the mutation parameter is set to 1 (which is dadi's default)
	# I'm specifying the migration rates m12 and m21 which are 2x the estimates from dadi
	# I also specify the population sizes after the split (there's no growth)
	# and the time when the two populations split (in coalescence this is when they merge)
	import os
	import dadi
	
	# python version
core="-m 1 2 1.5554 -m 2 1 0.1122 -n 1 0.1317 -n 2 8.4225 -ej 0.07205 2 1"
command=dadi.Misc.ms_command(theta=0.1026112, ns=(60,70), core=core, iter=12103)
ms_fs=dadi.Spectrum.from_ms_file(os.popen(command))
folded_sfs=dadi.Spectrum.fold(ms_fs)
dadi.Plotting.plot_single_2d_sfs(folded_sfs,vmin=0.005)


# dadi model without 

	#command line version
	# ms 176 12103 -t 0.1026112 -I 2 82 94 -m 1 2 1.5554 -m 2 1 0.1122 -n 1 0.1317 -n 2 8.4225 -ej 0.07205 2 1

fi

if [ "$TEXAS" = true ]; then
	# TEXAS
	# nsam = 144 (82 + 62)
	# best model = sym_mig_size
	# params: nu1a=5.65,nu2a=4.46,nu1b=1.01,nu2b=19.94,m=0.46,T1=1.34,T2=0.40
	# different growth rates with -g i GROWTHRATE
	# set past events with e (like -eg t i GROWTHRATE)


	ms 144 1000 -t THETA -I -n 


fi

if [ "$ALABAMA" = true ]; then
	# ALABAMA/LOUISIANA
	# nsam = 286 (96 + 94 + 96)
	ms 286 1000 -t THETA -I -n 

fi

