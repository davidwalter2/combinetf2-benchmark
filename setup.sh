SCRIPT_FILE_REL_PATH="${BASH_SOURCE[0]}"
if [[ "$SCRIPT_FILE_REL_PATH" == "" ]]; then
  SCRIPT_FILE_REL_PATH="${(%):-%N}"
fi
BENCHMARK_BASE=$( cd "$( dirname "${SCRIPT_FILE_REL_PATH}" )" && pwd )

source combinetf2/setup.sh
source wums/setup.sh

# source HiggsAnalysis/CombinedLimit/env_standalone.sh