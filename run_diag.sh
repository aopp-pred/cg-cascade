#!/bin/bash
#BSUB -q short-serial
#BSUB -n 1
#BSUB -R "select[maxmem >= 32000] rusage[mem=32000]"
#BSUB -M 32000000
#BSUB -W 04:00

set -e
set -u


#-----------------------------------------------------------------------
# Calculate a model start time from a job array index
#-----------------------------------------------------------------------
calc_processing_time () {
    local array_index="$1"
    local base_time=$(date -u -d "2009-04-06 01:00:00")
    local offset_hours=$(( ${array_index} - 1 ))
    date -u -d "${base_time}+${offset_hours} hours" -u +'%F %T'
}

#-----------------------------------------------------------------------
# Convert a date/time to a cascade index pair.
#
# Args:
#   datetime
#     The required date/time in the format "%Y-%m-%d %H:%M:%S".
# Globals:
#   None
# Returns:
#   path_index
#     Path component of the index.
#   time_index
#     Time component of the index.
# Exits:
#   0 on success.
#-----------------------------------------------------------------------
cg_datetime_to_indices () {
    local target_time=$(date -d "$1" -u +'%Y-%m-%d %H:%M:%S')
    local base_time='2009-04-06 01:00:00'
    local diff_seconds=$(( $(date -d "${target_time}" -u +'%s') - \
    $(date -d "${base_time}" -u +'%s') ))
    if [[ $(( ${diff_seconds} % 3600 )) -ne 0 ]]; then
        echo "target time must be a whole number of hours after ${base_time}" 1>&2
        return 1
    fi
    local pidx=$(( ${diff_seconds} / 3600 / 2 + 1))
    local tidx=0
    if [[ $(( (${diff_seconds} / 3600) % 2 )) -eq 1 ]]; then
        tidx=1
    fi
    echo ${pidx} ${tidx}
}

#-----------------------------------------------------------------------
# Convert a cascade index pair to a date/time.
#
# Args:
#   path_index
#     Path component of the index (1 <= pindex <= 120).
#   time_index
#     Time component of the index (0 or 1).
# Globals:
#   None
# Returns:
#   datetime
#     The date/time corresponding to the index pair in the format
#     "%Y-%m-%d %H:%M:%S".
# Exits:
#   0 on success.
#-----------------------------------------------------------------------
cg_indices_to_datetime () {
    local pindex="$1"
    local tindex="$2"
    local base_time=$(date -u -d '2009-04-06 01:00:00')
    date -u -d "${base_time} +$(( ((${pindex} - 1) * 2) + ${tindex} ))hours" -u +'%F %T'
}

#-----------------------------------------------------------------------
# Determine if a workaround for missing cloud fields is necessary for
# a particular start time.
#
# Args:
#   start_time
#     A start time in the format "YYYY-MM-DD HH:MM:SS".
# Returns:
#   None
# Exits:
#   0 if a workaround is required, 1 otherwise.
#-----------------------------------------------------------------------
cloud_workaround () {
    local cloud_missing=("2009-04-06 06:00:00"
                         "2009-04-06 12:00:00"
                         "2009-04-06 18:00:00"
                         "2009-04-07 00:00:00"
                         "2009-04-08 00:00:00"
                         "2009-04-09 00:00:00"
                         "2009-04-10 00:00:00"
                         "2009-04-11 00:00:00"
                         "2009-04-12 00:00:00"
                         "2009-04-13 00:00:00"
                         "2009-04-14 00:00:00"
                         "2009-04-15 00:00:00"
                         "2009-04-16 00:00:00")
    local start_time="$1"
    for workaround_time in "${cloud_missing[@]}"; do
        if [[ $workaround_time == $start_time ]]; then
            return 0
        fi
    done
    return 1
}

# Define the base directory for running the NCL programs:
readonly RUN_DIR="/group_workspaces/jasmin2/aopp/cg-cascade/cg07/system/mkinput"

# Change into the run directory and run the main program using MPI:
cd "$RUN_DIR"

# Run the NCL program with the required file indices:
start_time=$(calc_processing_time $LSB_JOBINDEX)
end_time=$(calc_processing_time $(( $LSB_JOBINDEX + 1 )))
start_index=($(cg_datetime_to_indices "$start_time"))
end_index=($(cg_datetime_to_indices "$end_time"))
if cloud_workaround "$start_time"; then
    lcloud="True"
else
    lcloud="False"
fi
echo "Running processor for cascade diagnostics for time range: ${start_time} - ${end_time} (lcloud=${lcloud})"
ncl coarsen_diag_cascade_manyt.ncl pindex1=${start_index[0]} tindex1=${start_index[1]} pindex2=${end_index[0]} tindex2=${end_index[1]} lcloud=${lcloud}

# Split the output file into one file per time-step:
output_dir="/group_workspaces/jasmin2/aopp/cg-cascade/cg07/data/diag_cas"
output_file_ncl="$output_dir/raw/CASCADE_WarmPool-4km_xfhfc_p${start_index[0]}.${start_index[1]}-p${end_index[0]}.${end_index[1]}_SCM_T639.nc"
if ! [[ -f "$output_file_ncl" ]]; then
    echo "error: expected output file does not exist"
    exit 1
fi
base_time=$(date -ud "$start_time")
if [[ ${end_index[0]} -eq 120 ]] && [[ ${end_index[1]} -eq 1 ]]; then
    # For the last time-step we should extract all slices:
    end_slice=4
else
    # For intermediate steps we only take the first 4 slices:
    end_slice=3
fi
for i in $(seq 0 $end_slice); do
    valid_time=$(date -u -d "${base_time}+$(( $i * 15 )) minutes" -u +'%Y%m%d%H%M%S')
    echo "Extracting time-step: $valid_time"
    ncks -d time,${i},${i} "$output_file_ncl" "$output_dir/diag_cas.${valid_time}.nc"
done
echo "Deleting original file: $output_file_ncl"
rm -f "$output_file_ncl"
