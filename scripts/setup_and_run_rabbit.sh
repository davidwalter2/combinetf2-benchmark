SCRIPT_FILE_REL_PATH="${BASH_SOURCE[0]}"
if [[ "$SCRIPT_FILE_REL_PATH" == "" ]]; then
  SCRIPT_FILE_REL_PATH="${(%):-%N}"
fi

source rabbit/setup.sh
source wums/setup.sh

export OMP_NUM_THREADS="1"
export OPENBLAS_NUM_THREADS=$((`nproc`>64 ? 64 : `nproc`))

echo "rabbit_fit.py $1/rabbit.hdf5 -t 0 --unblind --allowNegativePOI -o $2 $3"

rabbit_fit.py $1/rabbit.hdf5 -t 0 --unblind --allowNegativePOI -o $2 $3
