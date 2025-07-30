# run in singularity


project=$1

# scan_bins=(1000 10000 100000) 
scan_bins=(100000) 
scan_systs=(1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192) # 16384)
# scan_systs=(256 512 1024 2048 4096 8192 16384)
benchmark_command() {
    local command="$1"
    local memlog="$2"

    start=$(date +%s.%N)
    timeout 600s psrecord "$command" --interval 0.1 --include-children --log "$memlog"
    exit_code=$?

    # If it was killed by timeout (exit code 124)
    if [[ $exit_code -eq 124 ]]; then
        echo "Run killed after 600s, stopping this scan."
        return 1
    fi

    return 0
}

### rabbit
echo "==========================================="
suffix="rabbit_timing_$2"
# echo "nBins,nSyst,fit,mem_real" > $results
# source setup.sh
for nbins in "${scan_bins[@]}"; do
    for nsysts in "${scan_systs[@]}"; do
        echo "Benchmark Rabbit with $nbins" bins and $nsysts systematics

        model=$project/model_nBins${nbins}_nSysts${nsysts}
        memlog=$model/memlog_$suffix.txt
        command="rabbit_fit.py $model/rabbit.hdf5 -t 0 --unblind --allowNegativePOI -o $model/ --postfix $suffix"

        # rabbit_fit.py $model/rabbit.hdf5 -t 0 --unblind --allowNegativePOI -o $model/ --postfix $suffix # --useTFMinimizer  --minimizerMethod trust-exact

        # run; if duration >600 s the function returns 1 → break out
        if ! benchmark_command "$command" "$memlog"; then
            break          
        fi
    done
done
