import argparse

import json

import pyhf
import numpy
import gzip

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

    with gzip.open(args.input, "rt") as f:
        spec = json.load(f)
    
    pyhf.set_backend(args.backend, custom_optimizer=args.optimizer)

    workspace = pyhf.Workspace(spec)
    model = workspace.model()
    data = workspace.data(model)

    optimizer_args = dict()
    if args.optimizer == "scipy":
        optimizer_args["method"] = args.method

    bestfit, results = pyhf.infer.mle.fit(
        data, 
        model, 
        return_result_obj=True,
        **optimizer_args
    )
    logger.info(f"Result: {results}")

    # logger.info(f"Backend: {backend}")
    # logger.info(f"Optimizer: {optimizer}")

    logger.info(f"-2*log(L) = {results['fun']}")

    param_names = model.config.parameters
    for name, val, err in zip(param_names, bestfit, results["unc"]):
        logger.info(f"{name:20s} = {val:.4f} ± {err:.4f}")

if __name__ == "__main__":
    main()