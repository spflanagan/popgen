#!/bin/bash


cd "${0%/*}" # move to location of script
dir="dadi_scripts_abacus"

if [ ! -d "$dir" ] ; then mkdir $dir ; fi 

############ SET THESE PARAMS ############
fw_pops=('TXFW' 'LAFW' 'ALFW' 'FLLG')
sw_pops=('TXCC' 'ALST' 'ALST' 'FLCC')
fw_projs=(30 30 30 30) #(54 52 60 76)
sw_projs=(30 30 30 30) #(86 62 62 72)
models=('SI' 'IM' 'AM' 'SC' 'SI2N' 'SIG' 'SI2NG' 'IMG' 'IM2N' 'IM2m' 'IM2mG' 'IM2NG' 'AM2N' 'AMG' 'AM2m' 'AM2NG' 'AM2N2m' 'AM2mG' 'AM2N2mG' 'SCG' 'SC2N' 'SC2m' 'SC2NG' 'SC2N2m' 'SC2mG' 'SC2N2mG')


############ RUN THE ANALYSIS ############
# the python script takes these arguments in this order:
# snps_file, popi, prji, popj, prjj, x, y

for ((i=0; i<(${#fw_pops[@]}); ++i)); do
	for ((mod=0; mod<(${#models[@]}); ++mod)); do
		echo "#!/bin/bash" > ${dir}/${fw_pops[$i]}_${sw_pops[$i]}_${models[$mod]}.sh
		echo "cd \"\${0%/*}\"" >> ${dir}/${fw_pops[$i]}_${sw_pops[$i]}_${models[$mod]}.sh
		echo "cd /share/data/people/spf50/fwsw_results/dadi_results/" >> ${dir}/${fw_pops[$i]}_${sw_pops[$i]}_${models[$mod]}.sh
		echo "rangeX=\$1" >> ${dir}/${fw_pops[$i]}_${sw_pops[$i]}_${models[$mod]}.sh
		echo "rangeY=\$2" >> ${dir}/${fw_pops[$i]}_${sw_pops[$i]}_${models[$mod]}.sh
		echo "/share/apps/Miniconda2/bin/python /home/spf50/popgen/scripts/252_pairwise_dadi.py ${fw_pops[$i]}_${sw_pops[$i]}.dadi.snps '${fw_pops[$i]}' ${fw_projs[$i]} '${sw_pops[$i]}' ${fw_projs[$i]} \$rangeX \$rangeY ${models[$mod]} 2>&1 > ${fw_pops[$i]}_${sw_pops[$i]}/${fw_pops[$i]}_${sw_pops[$i]}_${models[$mod]}_\${rangeX}_\${rangeY}.log" >> ${dir}/${fw_pops[$i]}_${sw_pops[$i]}_${models[$mod]}.sh
	done
done

for file in ./${dir}/*sh; do chmod u+x ${file}; done



