SCRIPT_FILE_REL_PATH="${BASH_SOURCE[0]}"
if [[ "$SCRIPT_FILE_REL_PATH" == "" ]]; then
  SCRIPT_FILE_REL_PATH="${(%):-%N}"
fi
BENCHMARK_BASE=$( cd "$( dirname "${SCRIPT_FILE_REL_PATH}" )" && pwd )

source jax_rabbit/setup.sh
source wums/setup.sh

export PYTHONPATH="${BENCHMARK_BASE}:$PYTHONPATH"
export OPENBLAS_NUM_THREADS=$((`nproc --all`>64 ? 64 : `nproc --all`))

# source HiggsAnalysis/CombinedLimit/env_standalone.sh
