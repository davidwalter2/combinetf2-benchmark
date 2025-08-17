import pandas as pd
import matplotlib.pyplot as plt
import sys
import argparse
import os
from wums import logging, output_tools, plot_tools  # isort: skip
from matplotlib.lines import Line2D
import glob

ylabels = {
    "fit":"Time in second",
    "mem_real": "Peak memory in MB",
}
def create_plot(
        csv_file_dirs, 
        nBins,
        key,
        args,
        xlabel="Number of systematics",
        ):
    
    # Create the plot
    fig, ax1 = plot_tools.figure([0,1], xlabel, ylabels[key], xlim=args.xlim, ylim=args.ylim, logy=True, logx=True, logx_base=2)

    # Load the CSV files
    linestyles = ["-", "--"]
    for n, l, c, m in (
        # ("pyhf_numpy_minuit", "PyHF+Numpy+Minuit", "darkblue", "P"),
        # ("pyhf_jax_minuit", "PyHF+JAX+Minuit", "blue", "P"),
        # ("pyhf_pytorch_minuit", "PyHFt+pyTorch+Minuit", "deepskyblue", "P"),
        # ("pyhf_tensorflow_minuit", "PyHF+TF+Minuit", "cyan", "P"),
        # ("pyhf_numpy_scipy", "PyHF+Numpy+Scipy", "deepskyblue", "P"),
        # ("pyhf_jax_scipy", "PyHF+JAX+Scipy", "cyan", "P"),
        ("pyhf_numpy_minuit", "PyHF (np, minuit)", "navy", "^"),
        ("pyhf_jax_scipy", "PyHF (jax, scipy)", "dodgerblue", "v"),
        ("text2workspace", "text2workspace", "palegreen", "P"),
        # ("combine_10p2p0", "Combine", "green", "."),
        ("combine_10p2p0_v2", "combine", "green", "."),
        # ("pyhf_numpy_minuit", "PyHF", "blue", "P"),
        ("combinetf", "CombineTF", "brown", "o"),
        ("rabbit_singularity", "Rabbit", "orange", "*"),
        # ("rabbit", "Rabbit (virt. env.)", "orange", "x"),
        # ("rabbit_singularity", "Rabbit (sing.)", "orange", "*"),
        # ("rabbit_singularity_eager", "Rabbit (sing., eager)", "purple", "*"),
    ):
        for i, csv_file_dir in enumerate(csv_file_dirs):
            csv_file = f"{csv_file_dir}/timing_model_scaling_{n}.csv"

            files = glob.glob(f"{csv_file_dir}/*/pyhf/timing_model_scaling_{n}.csv")

            if len(files) > 0:
                dfs = [pd.read_csv(f) for f in files]
                dfs = [d for d in dfs if not d.empty]
                df = pd.concat(dfs, ignore_index=True)
            elif os.path.isfile(csv_file):
                df = pd.read_csv(csv_file)
            else:
                continue            
            
            # print(f"Successfully loaded {csv_file}")
            # print(f"Columns found: {', '.join(df.columns)}")

            df = df.loc[df["nBins"]==nBins]

            if len(df) == 0:
                continue

            df = df.loc[df["fit"] > 1.8]

            df = df.sort_values(by="nSyst")

            ### check if fits have actually succeeded
            if n.startswith("text2workspace"):
                # for text2workspace check if datacard.root is there
                def path_exists(row):
                    path = f"{csv_file_dir}/model_nBins{int(row['nBins'])}_nSysts{int(row['nSyst'])}/combine/datacard.root"
                    return os.path.isfile(path)

                df = df[df.apply(path_exists, axis=1)]
            elif n.startswith("combinetf"):
                def check_success(row):
                    path = f"{csv_file_dir}/model_nBins{int(row['nBins'])}_nSysts{int(row['nSyst'])}/fitresults_123456789.root"
                    if not os.path.isfile(path):
                        return False
                    import ROOT
                    try:
                        tfile = ROOT.TFile(path, "READ")
                        fitresults = tfile.Get("fitresults")
                        edmval = fitresults.edmval
                        if edmval > 10e-5:
                            print(f"Critical EDM value = {edmval}, fit may not have succeeded!")
                            return False
                        return True
                    except Exception as e:
                        return False
                df = df[df.apply(check_success, axis=1)]
            elif n.startswith("combine_"):
                def check_success(row):
                    path = f"{csv_file_dir}/model_nBins{int(row['nBins'])}_nSysts{int(row['nSyst'])}/combine/combine_logger.out"
                    if not os.path.isfile(path):
                        return False
                    try:
                        with open(path, 'rb') as f:
                            # Go to the end and read backwards to find the last line
                            f.seek(-2, os.SEEK_END)
                            while f.read(1) != b'\n':
                                f.seek(-2, os.SEEK_CUR)
                            last_line = f.readline().decode()
                    except Exception as e:
                        return False

                    return "Minimization success! status=0" in last_line
                df = df[df.apply(check_success, axis=1)]
            elif n.startswith("pyhf_") and n == "pyhf_jax_scipy":
                def check_success(row):
                    path = f"{csv_file_dir}/model_nBins{int(row['nBins'])}_nSysts{int(row['nSyst'])}/pyhf/out.log"
                    if not os.path.isfile(path):
                        return False
                    try:
                        with open(path, 'r') as f:
                            status_ok = False
                            success_ok = False
                            for line in f:
                                line = line.strip()
                                if line == "status: 0":
                                    status_ok = True
                                elif line == "success: True":
                                    success_ok = True
                                # early exit if both found
                                if status_ok and success_ok:
                                    return True
                    except Exception as e:
                        return False

                    return False
                df = df[df.apply(check_success, axis=1)]
            if len(df) == 0:
                continue

            # Ensure required columns exist
            required_cols = ['nSyst', key]
            for col in required_cols:
                if col not in df.columns:
                    print(f"Error: Required column '{col}' not found in CSV. Available columns: {df.columns.tolist()}")
                    sys.exit(1)
            
            # Convert columns to appropriate types
            df['nSyst'] = pd.to_numeric(df['nSyst'], errors='coerce')
            df[key] = pd.to_numeric(df[key], errors='coerce')
            
            # Remove any rows with missing values after conversion
            orig_len = len(df)
            df = df.dropna(subset=['nSyst', key])
            if len(df) < orig_len:
                print(f"Warning: Removed {orig_len - len(df)} rows with non-numeric values")
            
            # Group by nSyst (although seaborn will handle this for us in the plot)
            cpu_groups = df.groupby('nSyst')
            # print(f"Number of unique CPU values: {len(cpu_groups)}")
            # for cpu, group in cpu_groups:
            #     print(f"nSyst={cpu}: {len(group)} measurements, Avg time: {group[key].mean():.4f}")
            
            # if key == "fit" and "preparation" in df.keys():
            #     # ax1.plot(df["nSyst"].values, df["preparation"].values, label=l+"preparation" if i==0 else None, marker=m, linestyle="dotted", color=c)
            #     df["fit"] = df["fit"] + df["preparation"]

            ax1.plot(df["nSyst"].values, df[key].values, label=l if i==0 else None, marker=m, linestyle=linestyles[i], color=c)

    plot_tools.add_decor(
        ax1,
        args.title,
        args.subtitle,
        data=True,
        lumi=None,
        loc=args.titlePos,
        text_size=args.legSize,
        no_energy=True,
    )

    legend1 = ax1.legend(loc="upper left", ncols=args.legCols, columnspacing=1.4)

    if len(csv_file_dirs) > 1:
        line_solid = Line2D([0], [0], color='black', linestyle='-', label="no BB")#label='no bin by bin stat.')
        line_dashed = Line2D([0], [0], color='black', linestyle='--', label="BB-lite")#label='bin by bin stat.')

        # Second legend
        ax1.legend(handles=[line_solid, line_dashed], loc='lower right')
        # ax1.legend(handles=[line_solid, line_dashed], loc='upper center')

    ax1.add_artist(legend1)

    ax1.text(0.67, 0.02, r"$N^\mathrm{bins} = "+f"{nBins}$",  transform=ax1.transAxes, va="bottom", ha="right")

    outdir = args.output
    outfile = f"model_{key}_nBins{nBins}"
    if args.postfix:
        outfile += f"_{args.postfix}"
    if args.subtitle == "Preliminary":
        outfile += "_preliminary"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    plot_tools.save_pdf_and_png(outdir, outfile)

    output_tools.write_index_and_log(
        outdir,
        outfile,
        args=args,
    )


def main():
    # Default CSV file path
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--inputDirs",
        type=str,
        nargs="+",
        default="results/",
        help="Input csv file",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="results/",
        help="Output directory",
    )
    parser.add_argument(
        "--postfix",
        default=None,
        type=str,
        help="Postfix to append on output file name",
    )
    parser.add_argument(
        "--keys",
        type=str,
        nargs="+",
        default=["fit", "mem_real"],
        help="Key to take from the csv file to plot on as y coordinates",
    )
    parser.add_argument(
        "--xlim",
        type=float,
        nargs=2,
        help="Min and max values for x axis (if not specified, range set automatically)",
    )
    parser.add_argument(
        "--ylim",
        type=float,
        nargs=2,
        help="Min and max values for y axis (if not specified, range set automatically)",
    )
    parser.add_argument(
        "--title",
        default="Rabbit",
        type=str,
        help="Title to be printed in upper left",
    )
    parser.add_argument(
        "--subtitle",
        default="",
        type=str,
        help="Subtitle to be printed after title",
    )
    parser.add_argument("--titlePos", type=int, default=2, help="title position")
    parser.add_argument(
        "--legSize",
        type=str,
        default="small",
        help="Legend text size (small: axis ticks size, large: axis label size, number)",
    )
    parser.add_argument(
        "--legCols",
        type=int,
        default=1,
        help="Number of columns in the legend",
    )
    parser.add_argument(
        "--nBins",
        type=int,
        nargs="+",
        default=[1_000, 10_000, 100_000],
        help="Number of columns in the legend",
    )
    args = parser.parse_args()

    for key in args.keys:
        for nBins in args.nBins:
            create_plot(args.inputDirs, nBins, key, args)

if __name__ == "__main__":
    main()
