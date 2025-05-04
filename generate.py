import argparse

import numpy as np
from scipy.stats import truncnorm
from scipy.stats import norm
import hist
import uproot

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
        "--nEvents", 
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
    return parser.parse_args()


def random_normal(sigma, mu=0, lo=0.1, hi=0.9, size=1):
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
    return lambda x: amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)



class Model:

    # keep track of all systematics
    systematics = set()
    norm_systematics = set()

    def __init__(self, name, nevents, nbins, cdf=None, inv_cdf=None, mu=None, sigma=None):
        self.name=name

        self.nbins = int(nbins)
        self.nevents = int(nevents)

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

        self.norm_uncertainties = {}
        self.shape_systematics = {}
        
    def get_asimov(self):
        edges = np.linspace(0,1,self.nbins+1)
        return self.nevents * np.array([self.cdf(l,h) for l,h in zip(edges[:-1],edges[1:])])

    def get_samples(self):
        if self.inv_cdf:
            u = np.random.uniform(0,1,self.nevents)
            return self.inv_cdf(u)
        else:
            u = random_normal(self.sigma, self.mu, lo=0, hi=1,size=self.nevents)
            return u

    def generate_samples(self):
        self.samples = self.get_samples()

    def get_prediction(self, samples=None, weights=None, syst_name=None):
        if samples is None:
            if self.samples is not None:
                samples = self.samples
            else:
                samples = self.get_samples()
        if syst_name is not None:
            var = self.shape_systematics[syst_name](samples)

            arr_up = np.histogram(samples, weights=1+var, bins=self.nbins, range=[0,1])
            arr_dn = np.histogram(samples, weights=1-var, bins=self.nbins, range=[0,1])
            return arr_up, arr_dn

        return np.histogram(samples, weights=weights, bins=self.nbins, range=[0,1])

    def add_norm_uncertainty(self, name, size):
        self.norm_uncertainties[name] = size
        Model.norm_systematics.add(name)

    def add_systematic(self, name, basis, delta=None):
        # weight variations always have a coefficient 'norm' and a basis function f
        if delta is None:
            # pick norm random normal distributed
            delta = random_normal()
        self.shape_systematics[name] = lambda x, d=delta, f=basis: d * f(x)
        Model.systematics.add(name)
        

def get_data(models, theta, data_stat=True, vary_systematics=True):
    # data follows poisson distributed 
    central = 0
    for m in models:
        pred = m.get_asimov()
        logger.debug(f"sum({m.name})={np.sum(pred)}")


        if vary_systematics:
            # norm uncertainties
            morm_vars = [(1+size)**theta[syst] for syst, size in m.norm_uncertainties.items()]
            edges = np.linspace(0,1,m.nbins+1)
            centers = (edges[1:] - edges[:-1])/2. + edges[:-1]
            pred = pred * np.prod(morm_vars) * np.prod([(1+syst(centers))**theta[name] for name, syst in m.shape_systematics.items()], axis=0)
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
    np.random.seed(42)  # For reproducibility

    ### define signal and background models
    models = [
        Model("signal", 0.2*args.nEvents, args.nBins, mu=0.5,sigma=0.2),
        Model("bkg0", 0.1*args.nEvents, args.nBins, lambda x: x, lambda x: x),
        Model("bkg1", 0.1*args.nEvents, args.nBins, lambda x: x**2, lambda x: np.sqrt(x)),
        Model("bkg2", 0.1*args.nEvents, args.nBins, lambda x: 2*x-x**2, lambda x: 1 - np.sqrt(1-x)),
        Model("bkg3", 0.1*args.nEvents, args.nBins, mu=0.5, sigma=0.5),
        Model("bkg4", 0.1*args.nEvents, args.nBins, mu=0.2, sigma=0.5),
        Model("bkg5", 0.1*args.nEvents, args.nBins, mu=0.8, sigma=0.5),
        Model("bkg6", 0.1*args.nEvents, args.nBins, mu=0.4, sigma=0.25),
        Model("bkg7", 0.1*args.nEvents, args.nBins, mu=0.6, sigma=0.25),

    ]

    ### define systematic uncertainties
    # keep track of systematics to assign systematic shifts
    # systematics = ["norm", "norm_bkg"]
    # define normalization uncertainties
    #   1 fully correlated,
    #   1 correlated across all backgrounds
    #   8, 1 for each background
    for m in models:
        m.add_norm_uncertainty("norm", 0.01)
        
        if m.name=="signal":
            continue

        m.add_norm_uncertainty("norm_bkg", 0.02)
        m.add_norm_uncertainty(f"norm_{m.name}", 0.05)

        # systematics.append(f"norm_{m.name}")


    # define shape uncertainties
    # make systematic uncertainties with size of std=1%, but maximum 50%
    sigma, mu, lo, hi = 0.01, 0, 0, 0.5
    for i in range(args.nSystematics):
        # systematics.append(f"syst_{i}_local")
        # systematics.append(f"syst_{i}_local_bkg")

        basis = local_gaussian_basis()
        delta = random_normal(sigma, mu, lo, hi)

        basis_bkg = local_gaussian_basis()
        delta_bkg = random_normal(sigma, mu, lo, hi)
        for m in models:
            # systematics.append(f"syst_{i}_local_{m.name}")

            # fully correlated across all backgrounds -> 1
            m.add_systematic(f"syst_{i}_local", basis, delta)

            # uncorrelated for each process -> 8
            m.add_systematic(f"syst_{i}_local_{m.name}", local_gaussian_basis(), random_normal(sigma, mu, lo, hi),)

            if m.name=="signal":
                continue

            # correlated across backgrounds -> 1
            m.add_systematic(f"syst_{i}_local_bkg", basis_bkg, delta_bkg)


    # shift nuisances in data generation
    theta = {syst: np.random.normal(0, 1) for syst in Model.systematics}
    theta.update({syst: np.random.normal(0, 1) for syst in Model.norm_systematics})
    # theta = {syst: 0 for syst in Model.systematics}

    ### compute data, predictions, and systematic variations, and write them out
    for m in models:
        m.generate_samples()

    logger.info("=== add data ===")
    data = get_data(models, theta)

    ## combineTF2
    writer = tensorwriter.TensorWriter(
        sparse=False,
        systematic_type="log_normal",
    )
    writer.add_channel([hist.axis.Regular(args.nBins, 0,1, overflow=False, underflow=False)], "ch0")
    writer.add_data(data, "ch0")

    logger.info("=== add processes ===")
    for m in models:
        pred = m.get_prediction()[0]

        logger.debug(f"sum({m.name})={np.sum(pred)}")

        writer.add_process(pred, m.name, "ch0", signal=m.name=="signal")

        for n,s in m.norm_uncertainties.items():
            writer.add_lnN_systematic(n, m.name, "ch0", 1+s)

        for syst_name, syst in m.shape_systematics.items():
            syst_pred = m.get_prediction(syst_name=syst_name)
            writer.add_systematic(
                [syst_pred[0][0], syst_pred[1][0]],
                syst_name,
                m.name,
                "ch0",
            )

    writer.write(outfolder=directory, outfilename="combinetf2")


    ## CombineTF1

    ## Combine

    # generate root file
    with uproot.recreate(f"{directory}/combine/shapes.root") as f:
        h_data = hist.Hist(hist.axis.Regular(args.nBins, 0,1, overflow=False, underflow=False), storage=hist.storage.Double(), data=data)

        f[f"data_obs"] = h_data
        
        for m in models:

            f[f"{m.name}/{m.name}"] = m.get_prediction()

            for syst_name, syst in m.shape_systematics.items():
                syst_pred = m.get_prediction(syst_name=syst_name)
                
                f[f"{m.name}/{m.name}_{syst_name}up"] = syst_pred[0]
                f[f"{m.name}/{m.name}_{syst_name}down"] = syst_pred[1]

    # generate data card
    with open(f"{directory}/combine/datacard.txt", "w") as f:

        # Header
        f.write(f"imax 1  number of bins\n")
        f.write(f"jmax {len(models)-1}  number of processes minus 1\n")
        f.write(f"kmax {len(theta)}  number of nuisance parameters (explicitly defined below)\n")
        f.write(f"shapes data_obs * shapes.root $PROCESS\n\n")
        f.write(f"shapes * * shapes.root $PROCESS/$PROCESS $PROCESS/$PROCESS__$SYSTEMATIC\n\n")

        # Observations
        f.write("bin         ch0 \n")
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
            f.write(f"{key:<12} {' '.join(rows[key])}\n")
        f.write("\n")

        # Systematics
        for name in Model.norm_systematics:
            systype = "lnN"
            f.write(f"{name:<20} {systype:<10} ")
            for m in models:
                if name in m.norm_uncertainties.keys():
                    f.write(f"{m.norm_uncertainties[name]} ")
                else:
                    f.write("- ")
            f.write("\n")

        for name in Model.systematics:
            systype = "shape"
            f.write(f"{name:<20} {systype:<10} ")
            for m in models:
                if name in m.shape_systematics.keys():
                    f.write("1 ")
                else:
                    f.write("- ")
            f.write("\n")




if __name__ == "__main__":
    main()