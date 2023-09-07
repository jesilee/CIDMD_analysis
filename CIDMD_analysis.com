#!/bin/bash
#
# Purpose: this script is to post process the CIDMD simulations #
# Dependencies : python 3.6 
#		 postproc_Feb23.py
#		 get_gputime.py
#		 plot_ar_vel.py
#		 Molecule/mol_info.in


##= Removing old results.

rm ../results/*
rmdir ../results
mkdir ../results
rm plot_ar_vel.in



##= Extracting ar_velocities.

for i in $(ls -d ../calcs/*/chunk_0000/scr); do 
	echo ${i}/ar_vel.dat | tee -a plot_ar_vel.in
	(cd $i
	grep Ar vel.log  | awk '{ printf("%16.10f\n", sqrt($2*$2 + $3*$3 + $4*$4)) }' | nl > ar_vel.dat
	)
done

python3 plot_ar_vel.py 1




##= Extracting processing time.

find ../calcs/*/ -name 'run.out' | xargs egrep -i 'Total processing time' | awk '{print $4}' >> ../results/gputime.log





##= Postprocessing CIDMD simulations.

python3 postproc.py


rm gputime.log





