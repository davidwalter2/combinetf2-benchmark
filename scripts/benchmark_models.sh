
project=$1

# following are optional arguments passed into the fitter, e.g. '--noHessian'
shift 1
options=("$@")

if [[ ! -d "$project" ]]; then
    echo "Create project directory $project"
    mkdir -p "$project"
fi

source setup.sh
scripts="$BENCHMARK_BASE/scripts"

export APPTAINER_BIND="/tmp,/home/submit,/work/submit,/ceph/submit,/scratch/submit,/cvmfs,/etc/grid-security,/run" 

scan_bins=(10 100 1000) # 10000 100000)
scan_systs=(1 2 4 8 16 32 64 128 256 512 1024) # 2048 4096 8192)

source env_combinetf2/bin/activate


benchmark_command() {
    local command="$1"
    local memlog="$2"
    local nbins="$3"
    local nsysts="$4"
    local results="$5"
    local offset="${6:-0}"      # Optional offset, default 0
    local extra_csv="${7:-}"    # Optional extra CSV string

    start=$(date +%s.%N)
    psrecord "$command" --interval 0.1 --include-children --log "$memlog"
    end=$(date +%s.%N)

    duration=$(echo "$end - $start - $offset" | bc)
    max_real=$(awk 'NR > 1 { if ($3 > max) max = $3 } END { print max }' "$memlog")

    echo "Operation took $duration seconds, peak memory ${max_real} MB"
    echo "$nbins,$nsysts,$duration,$max_real${extra_csv:+,$extra_csv}" >> "$results"

    # If the run exceeded 600 s, return 1 so the caller can react.
    if (( $(echo "$duration > 600" | bc -l) )); then
        return 1              # non‑zero ⇒ “failure” for the caller
    fi
    return 0
}


### make models
source env_combinetf2/bin/activate
source setup.sh
for nbins in "${scan_bins[@]}"; do
    for nsysts in "${scan_systs[@]}"; do
        echo "Make model with $nbins" bins and $nsysts systematics

        model=$project/model_nBins${nbins}_nSysts${nsysts}

        python generate.py -o $model --nBins $nbins --nSystematics $nsysts

    done
done


### combineTF 2 (virt. env.)
echo "==========================================="
results=$project/timing_model_scaling_combinetf2.csv
echo "nBins,nSyst,fit,mem_real" > $results
source setup.sh
for nbins in "${scan_bins[@]}"; do
    for nsysts in "${scan_systs[@]}"; do
        echo "Benchmark CombineTF 2 (virt. env.) with $nbins" bins and $nsysts systematics

        model=$project/model_nBins${nbins}_nSysts${nsysts}
        memlog=$model/memlog_combinetf2.txt
        command="combinetf2_fit.py $model/combinetf2.hdf5 -t 0 --noBinByBinStat --unblind '.*' --allowNegativePOI -o $model/ --postfix combinetf2 $options"

        # run; if duration >600 s the function returns 1 → break out
        if ! benchmark_command "$command" "$memlog" "$nbins" "$nsysts" "$results"; then
            echo "Last run took more than 600s, stopping this scan."
            break          
        fi
    done
done


### combineTF 2 (singularity)
echo "==========================================="
results=$project/timing_model_scaling_combinetf2_singularity.csv
echo "nBins,nSyst,fit,mem_real" > $results

CONTAINER=/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling\:v44

# time needed to set up environment will be subtracted
start=$(date +%s.%N)
singularity run $CONTAINER $scripts/setup_combinetf2.sh
end=$(date +%s.%N)
offset=$(echo "$end - $start" | bc)
echo "Time for setting up environment took $offset seconds" 

for nbins in "${scan_bins[@]}"; do
    for nsysts in "${scan_systs[@]}"; do
        echo "Benchmark CombineTF 2 (singularity) with $nbins" bins and $nsysts systematics

        model=$project/model_nBins${nbins}_nSysts${nsysts}
        memlog=$model/memlog_combinetf2_singularity.txt
        command="singularity run \"$CONTAINER\" $scripts/setup_and_run_combinetf2.sh $model $model $options"

        # run; if duration >600 s the function returns 1 → break out
        if ! benchmark_command "$command" "$memlog" "$nbins" "$nsysts" "$results" "$offset"; then
            echo "Last run took more than 600s, stopping this scan."
            break          
        fi
    done
done


### combineTF 1
echo "==========================================="
results=$project/timing_model_scaling_combinetf1.csv
echo "nBins,nSyst,fit,mem_real" > $results
source /cvmfs/cms.cern.ch/cmsset_default.sh
# time needed to set up environment will be subtracted
start=$(date +%s.%N)
cmssw-cc7 --command-to-run "cd CMSSW_10_6_19_patch2/src/ ; cmsenv ; cd -"
end=$(date +%s.%N)
offset=$(echo "$end - $start" | bc)
echo "Time for setting up environment took $offset seconds" 
for nbins in "${scan_bins[@]}"; do
    for nsysts in "${scan_systs[@]}"; do
        echo "Benchmark CombineTF 1 with $nbins" bins and $nsysts systematics

        model=$project/model_nBins${nbins}_nSysts${nsysts}
        memlog=$model/memlog_combinetf1.txt
        command="cmssw-cc7 --command-to-run \"cd CMSSW_10_6_19_patch2/src/ ; cmsenv ; combinetf.py $model/combinetf1.hdf5 --outputDir $model -t 0 --unblind-value --unblind-fit-result --yes-i-really-really-mean-it $options; cd -\""

        # run; if duration >600 s the function returns 1 → break out
        if ! benchmark_command "$command" "$memlog" "$nbins" "$nsysts" "$results" "$offset"; then
            echo "Last run took more than 600s, stopping this scan."
            break          
        fi
    done
done


### combine 9.2.1 (via podman)
echo "==========================================="
results=$project/timing_model_scaling_combine_9p2p1_podman.csv
echo "nBins,nSyst,fit,mem_real,preparation" > $results

for nbins in "${scan_bins[@]}"; do
    for nsysts in "${scan_systs[@]}"; do
        echo "Benchmark Combine with $nbins" bins and $nsysts systematics

        model=$project/model_nBins${nbins}_nSysts${nsysts}

        cd $model/combine
        start=$(date +%s.%N)
        $scripts/setup_and_run_combine.sh text2workspace.py datacard.txt -m 125
        end=$(date +%s.%N)
        preparation=$(echo "$end - $start" | bc)

        memlog=$model/memlog_combine.txt
        command="$scripts/setup_and_run_combine.sh combine -M MultiDimFit -t 0 -d datacard.root --algo singles"

        # run; if duration >600 s the function returns 1 → break out
        if ! benchmark_command "$command" "$memlog" "$nbins" "$nsysts" "$results" "$0" "$preparation"; then
            echo "Last run took more than 600s, stopping this scan."
            break          
        fi

    done
done

cd -


### combine 10.2.1 (via setting environment variables)
echo "==========================================="
results=$project/timing_model_scaling_combine_10p2p0.csv

for nbins in "${scan_bins[@]}"; do
    for nsysts in "${scan_systs[@]}"; do
        echo "Benchmark Combine with $nbins" bins and $nsysts systematics

        model=$project/model_nBins${nbins}_nSysts${nsysts}
        memlog=$model/memlog_combine.txt
        
        cd $mode  l/combine
        start=$(date +%s.%N)
        $scripts/setup_and_run_combine.sh text2workspace.py datacard.txt -m 125
        end=$(date +%s.%N)
        preparation=$(echo "$end - $start" | bc)

        command="env -i bash -c 'source $scripts/setup_combine.sh; cd $model/combine; combine -M MultiDimFit -t 0 -d datacard.root --algo singles'"

        # run; if duration >600 s the function returns 1 → break out
        if ! benchmark_command "$command" "$memlog" "$nbins" "$nsysts" "$results" "0" "$preparation"; then
            echo "Last run took more than 600s, stopping this scan."
            break     
        fi
    done
done