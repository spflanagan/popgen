#!/bin/bash


cd "${0%/*}" # move to location of script
dir="dadi_scripts_abacus"

if [ ! -d "$dir" ] ; then mkdir $dir ; fi 

############ SET THESE PARAMS ############
pops=('FLLG' 'FLCC' 'ALFW' 'ALST' 'LAFW' 'TXFW' 'TXCC')
projs=(70 60 72 70 72 46 61)
models=('SI' 'IM' 'AM' 'SC' 'SI2N' 'SIG' 'SI2NG' 'IMG' 'IM2N' 'IM2m' 'IM2mG' 'IM2NG' 'AM2N' 'AMG' 'AM2m' 'AM2NG' 'AM2N2m' 'AM2mG' 'AM2N2mG' 'SCG' 'SC2N' 'SC2m' 'SC2NG' 'SC2N2m' 'SC2mG' 'SC2N2mG')


############ RUN THE ANALYSIS ############
# the python script takes these arguments in this order:
# snps_file, popi, prji, popj, prjj, x, y

for ((i=0; i<(${#pops[@]}-1); ++i)); do
	for ((j=(i+1); j<${#pops[@]}; ++j)); do
		for ((mod=0; mod<(${#models[@]}); ++mod)); do
			echo "#!/bin/bash" > ${dir}/${pops[$i]}_${pops[$j]}_${models[$mod]}.sh
			echo "cd \"\${0%/*}\"" >> ${dir}/${pops[$i]}_${pops[$j]}_${models[$mod]}.sh
			echo "cd /share/data/people/spf50/fwsw_results/dadi_results/" >> ${dir}/${pops[$i]}_${pops[$j]}_${models[$mod]}.sh
			echo "rangeX=\$1" >> ${dir}/${pops[$i]}_${pops[$j]}_${models[$mod]}.sh
			echo "rangeY=\$2" >> ${dir}/${pops[$i]}_${pops[$j]}_${models[$mod]}.sh
			echo "/share/apps/Miniconda2/bin/python /home/spf50/popgen/scripts/252_pairwise_dadi.py fwsw75.dadi.snps '${pops[$i]}' ${projs[$i]} '${pops[$j]}' ${projs[$j]} \$rangeX \$rangeY ${models[$mod]} 2>&1 > ${pops[$i]}_${pops[$j]}/${pops[$i]}-${pops[$j]}_${models[$mod]}_\${rangeX}_\${rangeY}.log" >> ${dir}/${pops[$i]}_${pops[$j]}_${models[$mod]}.sh
		done
	done
done

for file in ./${dir}/*sh; do chmod u+x ${file}; done



