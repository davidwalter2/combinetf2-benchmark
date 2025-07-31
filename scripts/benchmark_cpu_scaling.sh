#!/usr/bin/env bash

input=$1    # input directory for hdf5 files
output=$2   # output directory for csv files

# following are optional arguments passed into the fitter, e.g. '--noHessian'
shift 2
options=("$@")

###
# e.g. source scripts/benchmark_cpu_scaling.sh /ceph/submit/data/user/d/david_w/WMassAnalysis/combineResults/250502_mw/WMass_eta_pt_charge/ /ceph/submit/data/user/d/david_w/combineTF2_benchmark/250516_scaling_cpu/
###

scripts="$BENCHMARK_BASE/scripts"

if [[ ! -d "$output" ]]; then
    echo "Create output directory $output"
    mkdir -p "$output"
fi

export APPTAINER_BIND="/tmp,/home/submit,/work/submit,/ceph/submit,/scratch/submit,/cvmfs,/etc/grid-security,/run" 

# Create an array of CPUs to test, skip numbers larger than available CPUs
available_cpus=$(nproc)
echo "$available_cpus CPUs available"
numbers_scans=(768 512 384 256 128 64 32 16 8 4 2 1)
# numbers_scans=(32 16 8 4 2 1)
# numbers_scans=(256 128 64 32)
numbers=()
for n in "${numbers_scans[@]}"; do
    if (( n <= available_cpus )); then
        numbers+=("$n")
    fi
done


echo "==========================================="
echo "Benchmark Rabbit"
results=$output/timing_cpu_scaling_rabbit_submit80_jax.csv
echo "nCPU,time" > $results

# source env_rabbit/bin/activate
# source setup.sh

for n in "${numbers[@]}"; do
    echo "Now at $n CPUs"

    source $scripts/cgroup_restrict_cpus.sh "$n"

    start=$(date +%s.%N)
    # rabbit_fit.py inputs/rabbit.hdf5 -t 0 -o test/
    rabbit_fit.py $input/rabbit.hdf5 -t 0 --unblind --postfix submit80_jax -o $output/ $options
    end=$(date +%s.%N)
    duration=$(echo "$end - $start" | bc)
    echo "Operation took $duration seconds"

    # write out results
    echo "$n,$duration" >> $results
done

source $scripts/cgroup_restrict_cpus.sh numbers[0]

# deactivate


# echo "==========================================="
# echo "Benchmark Rabbit(singularity)"
# results=$output/timing_cpu_scaling_rabbit_singularity.csv
# echo "nCPU,time" > $results

# CONTAINER=/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling\:v44

# # time needed to set up environment will be subtracted
# start=$(date +%s.%N)
# singularity run $CONTAINER $scripts/setup_rabbit.sh
# end=$(date +%s.%N)
# offset=$(echo "$end - $start" | bc)
# echo "Time for setting up environment took $offset seconds" 

# for n in "${numbers[@]}"; do
#     echo "Now at $n CPUs"

#     source $scripts/cgroup_restrict_cpus.sh "$n"

#     start=$(date +%s.%N)
#     singularity run $CONTAINER $scripts/setup_and_run_rabbit.sh $input $output $options
#     end=$(date +%s.%N)

#     duration=$(echo "$end - $start - $offset" | bc)
#     echo "Operation took $duration seconds" 
#     echo "$n,$duration" >> $results
# done


# echo "==========================================="
# echo "Benchmark CombineTF"
# results=$output/timing_cpu_scaling_combinet_submit80.csv
# echo "nCPU,time" > $results

# source /cvmfs/cms.cern.ch/cmsset_default.sh

# # time needed to set up environment will be subtracted
# start=$(date +%s.%N)
# cmssw-cc7 --command-to-run "cd CMSSW_10_6_19_patch2/src/ ; cmsenv ; cd -"
# end=$(date +%s.%N)
# offset=$(echo "$end - $start" | bc)
# echo "Time for setting up environment took $offset seconds" 

# for n in "${numbers[@]}"; do
#     echo "Now at $n CPUs"

#     source $scripts/cgroup_restrict_cpus.sh "$n"

#     start=$(date +%s.%N)
    
#     cmssw-cc7 --command-to-run "cd CMSSW_10_6_19_patch2/src/ ; cmsenv ; python HiggsAnalysis/CombinedLimit/scripts/combinetf.py $input/combinetf.hdf5 --postfix submit80 --outputDir $output --binByBinStat -t 0 --unblind-value --unblind-fit-result --yes-i-really-really-mean-it $options; cd -"
    
#     end=$(date +%s.%N)
#     duration=$(echo "$end - $start - $offset" | bc)
#     echo "Operation took $duration seconds" 
#     echo "$n,$duration" >> $results

# done

# source $scripts/cgroup_restrict_cpus.sh 768
