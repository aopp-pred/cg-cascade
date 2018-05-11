#!/bin/bash
#
# mkinput produces coarse-grained fields over whole region. extract_block selects and saves a subregion for processing.
#
# User defines region of interest.
# Script loops over all files in full/ directory and extracts region of interest.
# Saves to file of same name in new directory.
#
#============================================

# define path to data
datadir=/group_workspaces/jasmin2/aopp/cg-cascade/cg07/data/scm_in/gswuv

# define region

#============================
### MARITIME CONTINENT ##
#region=mc
##lat bounds - must be float
#minimum_lat=-8.0
#maximum_lat=2.0
##lon bounds - must be float
#minimum_lon=100.0
#maximum_lon=123.0
#============================
### DODGY POINT
#region=dodgy
##lat bounds - must be float
#minimum_lat=-7.0
#maximum_lat=-6.7
##lon bounds - must be float
#minimum_lon=106.5
#maximum_lon=106.7
#============================
### WEST PACIFIC
#region=wp
##lat bounds - must be float
#minimum_lat=-8.0
#maximum_lat=2.0
##lon bounds - must be float
#minimum_lon=150.0
#maximum_lon=173.0
#============================
## WEST PACIFIC 2 - no land
region=wp2
#lat bounds - must be float
minimum_lat=-5.0
maximum_lat=5.0
#lon bounds - must be float
minimum_lon=155.0
maximum_lon=178.0
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
 ncea -d lat,$minimum_lat,$maximum_lat -d lon,$minimum_lon,$maximum_lon "$f" "$region/$f"
done


