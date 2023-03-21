#!/bin/bash
# script to combine and average DIN, DIP (winter) and Chl (growing season)
# written by Daniel Thewes, 10th Jan 2022
# example calls:
# ./DT_combavg.sh ## specify infile and outfile below!
# ./DT_combavg.sh <simid>
# ./DT_combavg.sh <simid> remove_hist(0/1) # removes global attribute "history" from outfile

# configure input
if (( "$#" > 0));then
  simid=$1
else
  simid="sns144-4g-CS-NEC"
fi

if (( "$#" > 1));then
  biophys_switch=$2
else
  biophys_switch=1 # 1: bio, 2: phys
fi
if (( "$#" > 2));then
  remove_hist=$3
else
  remove_hist=1
fi

# get filenames
if [ $biophys_switch == 1 ];then
  filetype="skillCS"
else
  filetype="MphysCS"
fi

infile='/work/ku0646/g260105/IR/sns-oe/simout-gpmeh/'${simid}'/extract_'${filetype}'_'${simid}'.2017-mm.nc'
outfile='/work/ku0646/g260105/IR/sns-oe/simout-gpmeh/'${simid}'/extract_'${filetype}'_'${simid}'.2017-avgout.nc'

echo $simid
if [ $biophys_switch == 1 ];then
  echo 'biology'
else
  echo 'physics'
fi

# for biology: extract time slabs (not necessary for physics)
if [ $biophys_switch == 1 ];then
  echo "extracting time slabs"
  ncks -O -v total_chlorophyll_calculator_result -d time,2,8 ${infile} chlaux.nc
  ncks -O -v EH_abioP_DIP -d time,11,1 ${infile} dipaux.nc
  ncks -O -v EH_abioP_DINO3 -d time,11,1 ${infile} dino3aux.nc
  ncks -O -v EH_abioP_DINH4 -d time,11,1 ${infile} dinh4aux.nc
fi

# rename variables
echo "renaming variables"
if [ $biophys_switch == 1 ];then
  ncrename -O -v EH_abioP_DINH4,DIN dinh4aux.nc
  ncrename -O -v EH_abioP_DINO3,DIN dino3aux.nc
  ncrename -O -v EH_abioP_DIP,DIP dipaux.nc
  ncrename -O -v total_chlorophyll_calculator_result,Chl chlaux.nc
else
  ncrename -O -v tempmean,temp ${infile} tempaux.nc
  ncrename -O -v saltmean,salt ${infile} saltaux.nc
fi

# add NH4 and NO3 to get DIN
if [ $biophys_switch == 1 ];then
echo "adding DIN together"
ncbo -O --op_typ=add dino3aux.nc dinh4aux.nc dinaux.nc
fi

# average over time
echo "averaging over time"
if [ $biophys_switch == 1 ];then
  ncwa -O -a time -v Chl chlaux.nc ${outfile}
  ncwa -A -a time -v DIN dinaux.nc ${outfile}
  ncwa -A -a time -v DIP dipaux.nc ${outfile}
else
  ncwa -O -a time -v temp tempaux.nc ${outfile}
  ncwa -O -a time -v salt saltaux.nc ${outfile}
fi

# add coordinates and bathymetry
echo "adding coordinates and bathymetry"
ncks -A -v lat ${infile} ${outfile}
ncks -A -v lon ${infile} ${outfile}
ncks -A -v bathymetry ${infile} ${outfile}

#remove nco history, if wanted (default yes)
if [ $remove_hist == 1 ];then
  echo "removing history"
  ncatted -O -h -a history_of_appended_files,global,d,, ${outfile}
  ncatted -O -h -a history,global,d,, ${outfile}
fi

# clean up
echo "cleaning up"
rm *aux.nc
echo ${outfile}
echo "script finished"
