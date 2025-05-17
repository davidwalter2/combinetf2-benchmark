SCRIPT_FILE_REL_PATH="${BASH_SOURCE[0]}"
if [[ "$SCRIPT_FILE_REL_PATH" == "" ]]; then
  SCRIPT_FILE_REL_PATH="${(%):-%N}"
fi

source combinetf2/setup.sh
source wums/setup.sh

export OMP_NUM_THREADS="1"
export OPENBLAS_NUM_THREADS=$(($1>64 ? 64 : $1))

start=$(date +%s%3N) 
combinetf2_fit.py $2/combinetf2.hdf5 -t 0 --unblind '*' -o $3 $5
end=$(date +%s%3N) 
elapsed_ms=$(($end - $start))
elapsed_sec="${elapsed_ms:0:-3}.${elapsed_ms: -3}"
echo "Operation took $elapsed_sec seconds" 
echo "$1,$elapsed_sec" >> $4
