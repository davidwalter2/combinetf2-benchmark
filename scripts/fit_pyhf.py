import argparse

import json

import pyhf
import numpy as np
import time

from wums import logging 

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
    parser.add_argument("-i", "--input", default="./", help="Input .json file with pyhf specs")
    # parser.add_argument("-o", "--output", default="./", help="output directory")

    parser.add_argument(
        "--backend", 
        type=str,
        default="numpy",
        choices=["numpy", "tensorflow", "pytorch", "jax"],
    )
    parser.add_argument(
        "--optimizer", 
        type=str,
        default=None,
        choices=["scipy", "minuit"],
    )
    parser.add_argument(
        "--method", 
        type=str,
        default=None,
        choices=["L-BFGS-B", "trust-krylov", "trust-exact"],
    )
    return parser.parse_args()

def main():
    args = make_parser()

    global logger
    logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)
    
    pyhf.set_backend(args.backend, custom_optimizer=args.optimizer)

    start_time = time.time()

    if args.input.endswith(".hdf5"):
        import h5py
        from wums import ioutils
        h5file = h5py.File(args.input, "r")
        spec = ioutils.pickle_load_h5py(h5file["spec"])
        for i, o in enumerate(spec["observations"]):
            spec["observations"][i]["data"] = o["data"].get()

        for i, o in enumerate(spec["channels"]):
            for j, p in enumerate(o["samples"]):
                spec["channels"][i]["samples"][j]["data"] = p["data"].get()
                for k, q in enumerate(p["modifiers"]):
                    if q["data"] is None:
                        continue
                    if isinstance(q["data"], dict):
                        spec["channels"][i]["samples"][j]["modifiers"][k]["data"]["hi_data"] = q["data"]["hi_data"].get()
                        spec["channels"][i]["samples"][j]["modifiers"][k]["data"]["lo_data"] = q["data"]["lo_data"].get()
                    else:
                        spec["channels"][i]["samples"][j]["modifiers"][k]["data"] = q["data"].get()

    elif args.input.endswith(".json.gz"):
        import gzip
        with gzip.open(args.input, "rt") as f:
            spec = json.load(f)
    else:
        with open(args.input, "rt") as f:
            spec = json.load(f)

    end_time = time.time()
    dtime = end_time - start_time
    logger.info(f"Time to load input: {dtime}s")

    start_time = time.time()

    workspace = pyhf.Workspace(spec, validate=False)
    model = workspace.model()

    end_time = time.time()
    dtime = end_time - start_time
    logger.info(f"Time to create model: {dtime}s")

    # does not work for np arrays
    # data = workspace.data(model)

    data = np.sum((workspace.observations[c] for c in model.config.channels))
    data = np.append(data, np.array(model.config.auxdata,dtype=int))

    optimizer_args = dict()
    if args.optimizer == "scipy":
        optimizer_args["method"] = args.method

    bestfit, results = pyhf.infer.mle.fit(
        data, 
        model, 
        return_result_obj=True,
        # return_uncertainties=True,
        **optimizer_args
    )

    if args.optimizer == "scipy":
        inv_hess = results.hess_inv
        n_params = len(bestfit)

        param_variances = np.array([inv_hess.matvec(np.eye(n_params)[i])[i] 
                                for i in range(n_params)])

    logger.info(f"Result: {results}")

    tensorlib, _ = pyhf.get_backend()
    pars, data = tensorlib.astensor(bestfit), tensorlib.astensor(data)
    # exp = model.expected_data(pars)

    # logger.info(f"Backend: {backend}")
    # logger.info(f"Optimizer: {optimizer}")
    nll = results['fun'] #-2*model.logpdf(bestfit, data)[0]
    logger.info(f"-2*log(L) = {nll}")


    # if len(pars) > 1:
    #     auxdata = model.expected_auxdata(pars)
    #     nll_constraint = -2*model.constraint_logpdf(auxdata, pars)
    # else:
    #     nll_constraint = 0

    # nll_data = nll - nll_constraint
    # # nll_main = model.mainlogpdf(data, bestfit)

    # # logger.info(f"-2*log(L) = {results['fun']}")
    # logger.info(f"-2*log(L) = {nll}")
    # logger.info(f"-2*log(L_constraint) = {nll_constraint} ")
    # logger.info(f"-2*log(L_data) = {nll_data} ")

    # logger.info(f"-2*log(L) = {-2*model.constraint_logpdf(bestfit, data)} from model.constraint_logpdf")

    param_names = [p for p in model.config._par_order if p != "MCstat"]

    for name, val, err in zip(param_names, bestfit, results["unc"]):
        logger.info(f"{name:20s} = {val:.4f} ± {err:.4f}")

if __name__ == "__main__":
    main()