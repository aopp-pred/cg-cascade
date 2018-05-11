#!/bin/bash
#
# mkinput produces coarse-grained fields over whole region. postproc_geostrophic_winds
# sets the geostrophic winds to zero or to the wind field to test sensitivity
# 
# Script loops over all files in full/ directory and processes winds.
# Saves to file of same name in new directory.
#
#============================================

# define path to data
datadir=/group_workspaces/jasmin2/aopp/cg-cascade/cg07/data/scm_in
#datadir=/group_workspaces/jasmin2/aopp/cg-cascade/cg07/system/mkinput/scm_in

# define region

#============================
# Set geostrophic winds to zu or v respectively
region=gswuv
#============================

# move to data
cd $datadir

pwd

# make directory for new data
mkdir $region

FILES=*.nc
for f in $FILES;
do
 echo "Processing $f"
 # do something on $f

 ncks -x -v ug "$f" "$region/tmp1"
 ncks -x -v vg "$region/tmp1" "$region/tmp2"
 ncap2 -s "ug=u" "$region/tmp2" "$region/tmp3"
 ncap2 -s "vg=v" "$region/tmp3" "$region/$f"

 rm "$region/tmp1" "$region/tmp2" "$region/tmp3"

done


