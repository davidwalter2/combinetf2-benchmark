import argparse

import numpy as np
from scipy.stats import truncnorm
from scipy.stats import norm
import hist
import uproot
import csv

from combinetf2 import tensorwriter

from wums import logging 

# from wrappers.to_combine import generate_datacard

# range is always [0,1] -> length of 1
# 

def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v",
        "--verbose",
        type=int,
        default=3,
        choices=[0, 1, 2, 3, 4],
        help="Set verbosity level with logging, the larger the more verbose",
    )
    parser.add_argument(
        "--noColorLogger", action="store_true", help="Do not use logging with colors"
    )
    parser.add_argument("-o", "--output", default="./", help="output directory")
    parser.add_argument(
        "--nEventsPerBin", 
        type=int, 
        help="Number of events for the toy model",
        default=1000
        )
    parser.add_argument(
        "--nBins", 
        type=int, 
        help="Number of bins for the toy model",
        default=10
        )
    parser.add_argument(
        "--nSystematics", 
        type=int, 
        help="Number of systematics for the toy model",
        default=1
        )
    parser.add_argument(
        "--seed", 
        type=int, 
        help="Seed for random number generation",
        default=42
        )
    parser.add_argument(
        "--binByBinStat", 
        action="store_true",
        help="Vary limited MC size",
        )

    return parser.parse_args()


def random_normal(sigma, mu=0, lo=-0.1, hi=0.1, size=1):
    a_, b_ = (lo-mu) / sigma, (hi-mu) / sigma
    rv = truncnorm(a_, b_, loc=mu, scale=sigma)
    return rv.rvs(size=size)


def local_gaussian_basis(mu=None, sigma=None, amplitude=1.0):
    """
    Returns a Gaussian bump centered at mu with width sigma.
    
    Parameters:
    - x : np.ndarray
        Array of bin centers or feature values.
    - mu : float
        Center of the Gaussian.
    - sigma : float
        Width of the Gaussian.
    - amplitude : float
        Amplitude of the Gaussian bump (default 1.0).

    Returns:
    - np.ndarray
        Gaussian function values evaluated at x.
    """
    if mu is None:
        mu = np.random.uniform()
    if sigma is None:
        sigma = np.random.uniform()
    return lambda x, a=amplitude: a * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def fourier_basis(frequency, phase):
    return lambda x, f=frequency, p=phase: np.sin(2*np.pi*f*(x + p))


class Model:

    # keep track of all systematics
    systematics = set()
    norm_systematics = set()

    def __init__(self, name, events_per_bin=1000, nbins=10, cdf=None, inv_cdf=None, mu=None, sigma=None, asimov=False):
        self.name=name

        self.nbins = int(nbins)
        self.nevents = int(events_per_bin) * nbins
        self.events_per_bin = int(events_per_bin)

        self.edges = np.linspace(0,1,self.nbins+1)
        self.centers = (self.edges[1:] - self.edges[:-1])/2. + self.edges[:-1]

        # models can be either specified by giving the inverse cdf transform to generate them from, in this case also the pdf has to be provided
        #   xor by giving a gaussian mu and sigma

        if cdf is None:
            # normalize to integral between [0,1]
            p = norm.cdf(1, loc=mu, scale=sigma) - norm.cdf(0, loc=mu, scale=sigma)
            self.cdf = lambda l,h: 1/p * (norm.cdf(h, loc=mu, scale=sigma) - norm.cdf(l, loc=mu, scale=sigma))
        else:
            self.cdf = lambda l,h: cdf(h) - cdf(l)

        self.inv_cdf = inv_cdf

        self.mu = mu
        self.sigma = sigma

        self.asimov=asimov # if asimov the MC events are just the expectations

        self.norm_uncertainties = {}
        self.shape_systematics = {}
        
    def get_asimov(self):
        return self.nevents * np.array([self.cdf(l,h) for l,h in zip(self.edges[:-1],self.edges[1:])])

    def get_samples(self):
        if self.inv_cdf:
            if self.asimov:
                u = np.linspace(0, 1, self.nevents, endpoint=False)
            else:
                u = np.random.uniform(0,1,self.nevents)
            return self.inv_cdf(u)
        else:
            u = random_normal(self.sigma, self.mu, lo=0, hi=1,size=self.nevents)
            return u

    def generate_samples(self):
        self.samples = self.get_samples()

    def get_prediction(self, samples=None, weights=None, syst_name=None):
        
        if self.asimov:
            # don't throw samples, just take expected event count
            arr = self.get_asimov()
            if syst_name is not None:
                var = 1+self.shape_systematics[syst_name](self.centers)
                return arr*var, arr/var
            return arr
        
        if samples is None:
            if self.samples is not None:
                samples = self.samples
            else:
                samples = self.get_samples()
        if syst_name is not None:
            var = 1+self.shape_systematics[syst_name](samples)

            arr_up = np.histogram(samples, weights=var, bins=self.nbins, range=[0,1])[0]
            arr_dn = np.histogram(samples, weights=1./var, bins=self.nbins, range=[0,1])[0]
            return arr_up, arr_dn

        return np.histogram(samples, weights=weights, bins=self.nbins, range=[0,1])[0]

    def add_norm_uncertainty(self, name, size):
        self.norm_uncertainties[name] = size
        Model.norm_systematics.add(name)

    def add_systematic(self, name, basis, delta=None):
        # weight variations always have a coefficient 'norm' and a basis function f
        if delta is None:
            # pick norm random normal distributed
            delta = random_normal()[0]
        self.shape_systematics[name] = lambda x, d=delta, f=basis: d * f(x)
        Model.systematics.add(name)
        

def get_data(models, mu, theta, data_stat=True, vary_systematics=True):
    # data follows poisson distributed 
    central = 0
    for m in models:
        pred = m.get_asimov()
        logger.debug(f"sum({m.name})={np.sum(pred)}")

        if m.name=="sig":
            pred *= mu

        if vary_systematics:
            # norm uncertainties
            morm_vars = [(1+size)**theta[syst] for syst, size in m.norm_uncertainties.items()]
            shape_vars = [(1+syst(m.centers))**theta[name] for name, syst in m.shape_systematics.items()]
            pred = pred * np.prod(morm_vars) * np.prod(shape_vars, axis=0)
        central += pred
    if data_stat:
        counts = np.random.poisson(central)
    else:
        counts = central

    return counts


def main():
    args = make_parser()

    global logger
    logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

    directory = args.output
    if directory == "":
        directory = "./"

    # Generate random data for filling
    np.random.seed(args.seed)  # For reproducibility

    ### define signal and background models
    kwargs = dict(
        events_per_bin=0.2*args.nEventsPerBin, 
        nbins= args.nBins, 
        cdf=lambda x: x, 
        inv_cdf=lambda x: x,
        asimov=not args.binByBinStat
    )
    models = [
        Model("sig", **kwargs),
        Model("bkg0", **kwargs),
        Model("bkg1", **kwargs),
        Model("bkg2", **kwargs),
        Model("bkg3", **kwargs),
    ]

    # ### define systematic uncertainties
    # # keep track of systematics to assign systematic shifts
    # # define normalization uncertainties
    # #   1 fully correlated,
    # #   1 correlated across all backgrounds
    # #   8, 1 for each background
    # for m in models:
    #     m.add_norm_uncertainty("norm", 0.01)
        
    #     if m.name=="sig":
    #         continue

    #     m.add_norm_uncertainty("norm_bkg", 0.02)
    #     m.add_norm_uncertainty(f"norm_{m.name}", 0.05)

    # define shape uncertainties
    # make systematic uncertainties with size of std=1%, but maximum 50%
    for i in range(args.nSystematics):

        # basis = local_gaussian_basis()
        # frequency
        frequency =  np.random.uniform(0.25, args.nBins/2.)
        phase = np.random.uniform(0, 1)
        basis = fourier_basis(frequency, phase)

        for m in models:
            # fully correlated across all processes, 
            #   the detla is random meaning that each process is affected differently
            delta = random_normal(sigma=0.02, mu=0, lo=-0.1, hi=0.1)[0]
            m.add_systematic(f"syst_{i}", basis, delta)


    # shift nuisances in data generation
    theta = {syst: np.random.normal(0, 1) for syst in Model.systematics}
    theta.update({syst: np.random.normal(0, 1) for syst in Model.norm_systematics})
    # random signal strength
    mu = random_normal(0.05, 1, 0.5, 1.5)[0]

    ### compute data, predictions, and systematic variations, and write them out
    if args.binByBinStat:
        for m in models:
            m.generate_samples()

    logger.info("=== add data ===")
    data = get_data(models, mu, theta)

    ## combineTF1/2
    writer = tensorwriter.TensorWriter(
        sparse=False,
        systematic_type="log_normal",
    )
    writer.add_channel([hist.axis.Regular(args.nBins, 0,1, overflow=False, underflow=False, name="x")], "ch0")
    writer.add_data(data, "ch0")

    logger.info("=== add processes ===")
    for m in models:
        pred = m.get_prediction()

        logger.debug(f"sum({m.name})={np.sum(pred)}")

        writer.add_process(pred, m.name, "ch0", signal=m.name=="sig")

        for n,s in m.norm_uncertainties.items():
            writer.add_lnN_systematic(n, m.name, "ch0", 1+s)

        for syst_name, syst in m.shape_systematics.items():
            syst_pred = m.get_prediction(syst_name=syst_name)
            writer.add_systematic(
                [syst_pred[0], syst_pred[1]],
                syst_name,
                m.name,
                "ch0",
            )

    writer.write(outfolder=directory, outfilename="combinetf2")

    writer.symmetric_tensor = False # for combinetf1
    writer.write(outfolder=directory, outfilename="combinetf1")

    ## pyHF


    ## Combine
    # https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/

    # generate root file
    axis = hist.axis.Regular(args.nBins, 0,1, overflow=False, underflow=False)
    with uproot.recreate(f"{directory}/combine/shapes.root") as f:
        h_data = hist.Hist(axis, storage=hist.storage.Double(), data=data)

        f[f"data_obs"] = h_data
        
        for m in models:

            h_pred = hist.Hist(
                axis, 
                storage=hist.storage.Weight() if args.binByBinStat else hist.storage.Double(), 
                data=m.get_prediction()
            )

            f[f"{m.name}/{m.name}"] = h_pred

            for syst_name in m.shape_systematics.keys():
                syst_pred = m.get_prediction(syst_name=syst_name)
                
                h_up = hist.Hist(
                    axis, 
                    storage=hist.storage.Double(), 
                    data=syst_pred[0]
                )
                
                h_down = hist.Hist(
                    axis, 
                    storage=hist.storage.Double(), 
                    data=syst_pred[1]
                )

                f[f"{m.name}/{syst_name}Up"] = h_up
                f[f"{m.name}/{syst_name}Down"] = h_down

    # generate data card
    with open(f"{directory}/combine/datacard.txt", "w") as f:

        # Header
        f.write(f"imax 1  number of bins\n")
        f.write(f"jmax {len(models)-1}  number of processes minus 1\n")
        f.write(f"kmax {len(theta)}  number of nuisance parameters (explicitly defined below)\n")
        f.write(f"shapes data_obs * shapes.root $PROCESS\n\n")
        f.write(f"shapes * * shapes.root $PROCESS/$PROCESS $PROCESS/$SYSTEMATIC\n\n")

        # Observations
        f.write("bin         bin1 \n")
        f.write(f"observation {sum(data)}\n\n")

        # Processes and rates
        rows = {
            "bin": [],
            "process": [],
            "process_index": [],
            "rate": []
        }

        # Processes and rates
        rows = {
            "bin": [],
            "process": [],
            "process_index": [],
            "rate": []
        }

        for i, m in enumerate(models):
            rows["bin"].append("bin1")
            rows["process"].append(m.name)
            rows["process_index"].append(str(i))
            rows["rate"].append("-1")  # -1 tells Combine to get shape rate from ROOT

        for key in ["bin", "process", "process_index", "rate"]:
            f.write(f"{key.replace('process_index','process'):<12} {' '.join(rows[key])}\n")
        f.write("\n")

        # Systematics
        for name in Model.norm_systematics:
            systype = "lnN"
            f.write(f"{name:<20} {systype:<10} ")
            for m in models:
                if name in m.norm_uncertainties.keys():
                    f.write(f"{1 + m.norm_uncertainties[name]} ")
                else:
                    f.write("- ")
            f.write("\n")

        for name in Model.systematics:
            systype = "shapeN2" # shapeN2 corresponds to bin by bin lnN variations (what is done in combineTF1/2)
            f.write(f"{name:<20} {systype:<10} ")
            for m in models:
                if name in m.shape_systematics.keys():
                    f.write("1 ")
                else:
                    f.write("- ")
            f.write("\n")

        if args.binByBinStat:
            f.write(f"* autoMCStats -1 1\n\n")

    # subprocess.run(f"text2workspace.py {directory}/combine/datacard.txt -o {directory}combine/datacard.hdf5 -m 125", shell=True, check=True)

    ## Write out parameters
    params = {"sig": mu, **theta}

    with open(f"{directory}/params.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["param", "value"])  # Header
        for key, value in params.items():
            writer.writerow([key, value])


if __name__ == "__main__":
    main()