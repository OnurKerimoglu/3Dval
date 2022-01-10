#!/bin/bash
# script to combine and average DIN, DIP (winter) and Chl (growing season)
# written by Daniel Thewes, 10th Jan 2022
# example calls:
# ./DT_combavg.sh ## specify infile and outfile below!
# ./DT_combavg.sh <infile> <outfile>
# ./DT_combavg.sh <infile> <outfile> remove_hist(0/1) # removes global attribute "history" from outfile

if (( "$#" > 0));then
  infile=$1
else
  infile='/work/ku0646/g260105/IR/sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-ICGEMO-CS-BCdcsmP-rivWS/extract_skillCS_sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-ICGEMO-CS-BCdcsmP-rivWS.2014-mm.nc'
fi

if (( "$#" > 1));then
  outfile=$2
else
  outfile='/work/ku0646/g260105/IR/sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-ICGEMO-CS-BCdcsmP-rivWS/extract_skillCS_sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-ICGEMO-CS-BCdcsmP-rivWS.2014-avgout.nc'
fi

if (( "$#" > 2));then
  remove_hist=$3
else
  remove_hist=1
fi

echo "extracting time slabs"
ncks -O -v total_chlorophyll_calculator_result -d time,2,8 ${infile} chlaux.nc
ncks -O -v EH_abioP_DIP -d time,11,1 ${infile} dipaux.nc
ncks -O -v EH_abioP_DINO3 -d time,11,1 ${infile} dino3aux.nc
ncks -O -v EH_abioP_DINH4 -d time,11,1 ${infile} dinh4aux.nc

echo "renaming variables"
ncrename -O -v EH_abioP_DINH4,DIN dinh4aux.nc
ncrename -O -v EH_abioP_DINO3,DIN dino3aux.nc
ncrename -O -v EH_abioP_DIP,DIP dipaux.nc
ncrename -O -v total_chlorophyll_calculator_result,Chl chlaux.nc

echo "adding DIN together"
ncbo -O --op_typ=add dino3aux.nc dinh4aux.nc dinaux.nc

echo "averaging over time"
ncwa -O -a time -v Chl chlaux.nc ${outfile}
ncwa -A -a time -v DIN dinaux.nc ${outfile}
ncwa -A -a time -v DIP dipaux.nc ${outfile}

echo "adding coordinates and bathymetry"
ncks -A -v lat ${infile} ${outfile}
ncks -A -v lon ${infile} ${outfile}
ncks -A -v bathymetry ${infile} ${outfile}

if [ $remove_hist == 1 ];then
  echo "removing history"
  ncatted -O -h -a history_of_appended_files,global,d,, ${outfile}
  ncatted -O -h -a history,global,d,, ${outfile}
fi
echo "removing auxilaries"
rm *aux.nc

echo "script finished"
