#!/usr/bin/env bash

input=$1    # input directory for csv files
input="/ceph/submit/data/user/d/david_w/WMassAnalysis/combineResults/250502_mw/WMass_eta_pt_charge/"
output=$2   # output directory for csv files

###
# e.g. source scripts/benchmark_cpu_scaling.sh /ceph/submit/data/user/d/david_w/WMassAnalysis/combineResults/250502_mw/WMass_eta_pt_charge/ /ceph/submit/data/user/d/david_w/combineTF2_benchmark/250516_scaling_cpu/
###

scripts=$BENCHMARK_BASE/scripts

if [[ ! -d "$output" ]]; then
    echo "Create output directory $output"
    mkdir -p "$output"
fi

export APPTAINER_BIND="/tmp,/home/submit,/work/submit,/ceph/submit,/scratch/submit,/cvmfs,/etc/grid-security,/run" 
numbers=(768 512 256 128 64 32 16 8 4 2 1)

echo "==========================================="
echo "Benchmark CombineTF2"
results=$output/timing_cpu_scaling_combinetf2.csv
echo "nCPU,time" > $results

source env_combinetf2/bin/activate
source setup.sh

for n in "${numbers[@]}"; do
    echo "Now at $n CPUs"

    if [[ "$n" -ne 768 ]]; then
        source $scripts/cgroup_restrict_cpus.sh "$n"
    fi

    start=$(date +%s.%N)
    # combinetf2_fit.py inputs/combinetf2.hdf5 -t 0 -o test/
    combinetf2_fit.py $input/combinetf2.hdf5 -t 0 --unblind 'r:.*' -o test/
    end=$(date +%s.%N)
    duration=$(echo "$end - $start" | bc)
    echo "Operation took $duration seconds"

    # write out results
    echo "$n,$duration" >> $results
done

deactivate

echo "==========================================="
echo "Benchmark CombineTF2(singularity)"
results=$output/timing_cpu_scaling_combinetf2_singularity.csv
echo "nCPU,time" > $results

CONTAINER=/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling\:v44
for n in "${numbers[@]}"; do
    echo "Now at $n CPUs"

    if [[ "$n" -ne 768 ]]; then
        source $scripts/cgroup_restrict_cpus.sh "$n"
    fi

    singularity run $CONTAINER $scripts/setup_and_run_combinetf2.sh $n $input $results
done


echo "==========================================="
echo "Benchmark CombineTF1"
results=$output/timing_cpu_scaling_combinetf1.csv
echo "nCPU,time" > $results

source /cvmfs/cms.cern.ch/cmsset_default.sh

# time needed to set up environment will be subtracted
start=$(date +%s.%N)
cmssw-cc7 --command-to-run "cd CMSSW_10_6_19_patch2/src/ ; cmsenv ; cd -"
end=$(date +%s.%N)
offset=$(echo "$end - $start" | bc)
echo "Time for setting up environment took $offset seconds" 

for n in "${numbers[@]}"; do
    echo "Now at $n CPUs"

    if [[ "$n" -ne 768 ]]; then
        source $scripts/cgroup_restrict_cpus.sh "$n"
    fi

    start=$(date +%s.%N)
    
    cmssw-cc7 --command-to-run "cd CMSSW_10_6_19_patch2/src/ ; cmsenv ; combinetf.py $input/combinetf1.hdf5 --binByBinStat -t 0 --unblind-value --unblind-fit-result --yes-i-really-really-mean-it --noHessian; cd -"
    
    end=$(date +%s.%N)
    duration=$(echo "$end - $start - $offset" | bc)
    echo "Operation took $duration seconds" 
    echo "$n,$duration" >> $results

done
