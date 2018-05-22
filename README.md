# Coarse Graining Input Generator

NCL scripts that take CASCADE data as input, together with a reference ECMWF IFS file.
The scripts produce a coarse-grained version of CASCADE, on the grid of the ECMWF IFS reference file, suitable for driving the OpenIFS Single Column Model.

The chosen CASCADE dataset is the Tropical Pacific/Warm Pool dataset at 4km resolution, with 3D Smagorinsky mixing. The ID for this dataset is xfhfc


SCRIPTS:

coarsen_cascade_manyt.ncl

The file coarsen_cascade_manyt.ncl is the driver script. It takes two consecutive CASCADE data dumps (separated by one hour), coarse-grains these data, and interpolates from the 1-hr CASCADE separation to the 15-minute timestep required by the SCM.

The script is called with variables to indicate which two CASCADE time steps are to be coarse-grained and interpolated between:
ncl coarsen_cascade_manyt.ncl pindex1=I pindex2=J tindex1=i tindex2=j lcloud=True/False

The indices refer to the structure of CASCADE files. A typical CASCADE filename looks like this:
CASCADE_WarmPool-4km_xfhfc_p1_10.nc
where the name of the dataset and dataset ID are followed by two identifiers: px_y. The number x indicates the timestamp for the file. Each file has two consectutive fields, separated by one hour. The number y refers to the variable contained in the file, following the UK Met Office data conventions.

In the usage of the driver script, pindex (which can take interger values between 1 and 120) refers to the timestamp number in the CASCADE filename, while tindex (which can take the value of 0 or 1) refers to which of the two timesteps in the CASCADE file is required.

The final variable input is the logical flag lcloud. This is necessary because of an idiosyncracy in the CASCADE data, whereby cloud variables are missing every sixth hour. If lcloud is 'true', the script treats the special case of missing cloud data, and instead selects the next nearest file and interpolates over a 2-hour window.

Within the ncl script, the user can specify the bounds of the region to be coarse-grained. The default is to coarse-grain the entire CASCADE region, though subregions can be selected for testing if required.

coarsen_diag_cascade_manyt.ncl

This behaves as for coarsen_cascade_manyt.ncl, except it coarse-grains CASCADE fields for diagnostic purposes that are not needed to drive the SCM. For example, producing coarse-grained precipitation, OLR, ..., fields, for comparison with the output of the SCM.


AUTOMATION SCRIPTS

It is convenient to automate the calling of the driver script. The following scripts automate the procedure of calling coarsen_cascade_manyt.ncl. The user supplies a list of time IDs which are converted by the shell script into p and t indices, and the correct lcloud flag set (according to the time).

The time IDs simply correspond to the number of hours after the start of the CASCADE run. e.g.

./u_run myexp 1

Will produce the coarse-grained fields for times between 01:00 and 02:00 after the start of the CASCADE simulation. The argument 'myexp' is a user-defined name to keep track of the coarse-graining. It will be included in log files. On JASMIN it is reasonable to call 12 timesteps at a time, to be run in parallel:
./u_run myexp 1-12

The similar file u_run_diag.sh can be used to call the coarsen_diag_cascade_manyt.ncl file.


POSTPROCESSING SCRIPTS

Having coarse-grained the entire CASCADE region, it is convenient to select sub-regions over which the SCM can be run. This can be achieved using:

extract_block.sh

Within the script, define a path to the input data, the min/max lat/lon bounds of the region of interest, and give the region a name. The region will be extracted and written to a folder with this name.

Furthermore, it was interesting to consider the impact of using 'incorrect' geostrophic winds on the evolution of the SCM, to emphasise the benefits of deriving the geostrophic winds from the CASCADE data. This was achieved by creating a copy of the data, and simply overwriting the correct geostrophic winds with approximations.

postproc_geostrophic_winds_0.sh
postproc_geostrophic_winds_uv.sh


