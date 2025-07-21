# Benchmark tools for binned profile maximum likelihood fitting

The following was done on submit82.mit.edu

## Make synthetic inputs
Run 
```bash
python generate.py --nBins 100 --nSystematics 100 -o inputs/
```

## Rabbit
Running the fit
``` bash 
rabbit_fit.py inputs/rabbit.hdf5 -o inputs/ -t 0 --unblind '*' --allowNegativePOI --noBinByBinStat  --saveHists --saveHistsPerProcess --computeHistErrors -m Basemodel
```
The `--allowNegativePOI` is needed to not square the signal strength parameter (to be fixed in rabbit)

Parameter pulls and constraints
```bash
rabbit_print_pulls_and_constraints.py inputs/fitresults.hdf5
```

Plots
```bash
rabbit_plot_hists.py inputs/fitresults.hdf5 -m Basemodel -o ~/public_html/combinetf2-benchmark/250517_test --titlePos 0 --extraTextLoc 0.03 0.97 --subtitle Preliminary --config style_config.py --yscale 1.4 --prefit --rrange 0.8 1.2
rabbit_plot_hists.py inputs/fitresults.hdf5 -m Basemodel -o ~/public_html/combinetf2-benchmark/250517_test --titlePos 0 --extraTextLoc 0.03 0.97 --subtitle Preliminary --config style_config.py --yscale 1.4
```

## Combine

Setting up environment (outside singularity)
```bash
source scripts/setup_combine.sh
```

```bash
cd test/combine
text2workspace.py datacard.txt -m 125
combine -M MultiDimFit -t 0 -d datacard.root --algo singles
```

Investigate the fit parameters
```bash
combine -M FitDiagnostics -t 0 -d datacard.root
root -l fitDiagnosticsTest.root
root [1] fit_s->Print()
```

## CombineTF
Setting up environment (outside singularity)
```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
APPTAINER_BIND="/tmp,/home/submit,/work/submit,/scratch/submit,/ceph/submit/,/cvmfs,/etc/grid-security,/run" cmssw-cc7
export SCRAM_ARCH="slc7_amd64_gcc700"
cmsrel CMSSW_10_6_19_patch2
cd CMSSW_10_6_19_patch2/src/
cmsenv
git clone -b tensorflowfit git@github.com:bendavid/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
scram b -j 8
```

Converting the combine datacards to combineTF 
```bash
cd text/
# text2hdf5 datacard.txt 
combinetf.py combinetf.hdf5 --binByBinStat -t 0 --unblind-value --unblind-fit-result --yes-i-really-really-mean-it
```
