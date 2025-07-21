import pandas as pd
import matplotlib.pyplot as plt
import sys
import argparse
import os
from wums import logging, output_tools, plot_tools  # isort: skip
from matplotlib.lines import Line2D

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
        ("pyhf_numpy_minuit", "PyHF", "blue", "P"),
        # ("pyhf_numpy_minuit", "PyHF+Numpy+Minuit", "darkblue", "P"),
        # ("pyhf_jax_minuit", "PyHF+JAX+Minuit", "blue", "P"),
        # ("pyhf_pytorch_minuit", "PyHFt+pyTorch+Minuit", "deepskyblue", "P"),
        # ("pyhf_tensorflow_minuit", "PyHF+TF+Minuit", "cyan", "P"),
        # ("pyhf_numpy_scipy", "PyHF+Numpy+Scipy", "deepskyblue", "P"),
        # ("pyhf_jax_scipy", "PyHF+JAX+Scipy", "cyan", "P"),
        ("text2workspace", "text2workspace", "olive", "."),
        # ("combine_10p2p0", "Combine", "green", "."),
        ("combine_10p2p0_v2", "combine", "green", "."),
        ("combinetf", "CombineTF", "brown", "o"),
        ("rabbit_singularity", "Rabbit", "orange", "*"),
        # ("rabbit", "Rabbit (virt. env.)", "orange", "x"),
        # ("rabbit_singularity", "Rabbit (sing.)", "orange", "*"),
        # ("rabbit_singularity_eager", "Rabbit (sing., eager)", "purple", "*"),
    ):
        for i, csv_file_dir in enumerate(csv_file_dirs):
            csv_file = f"{csv_file_dir}/timing_model_scaling_{n}.csv"

            if not os.path.isfile(csv_file):
                continue
            
            df = pd.read_csv(csv_file)
            # print(f"Successfully loaded {csv_file}")
            # print(f"Columns found: {', '.join(df.columns)}")

            df = df.loc[df["nBins"]==nBins]

            if len(df) == 0:
                continue

            df = df.loc[df["fit"] > 1.8]

            ### check if fits have actually succeeded
            if n.startswith("text2workspace"):
                # for text2workspace check if datacard.root is there
                def path_exists(row):
                    path = f"{csv_file_dir}/model_nBins{int(row['nBins'])}_nSysts{int(row['nSyst'])}/combine/datacard.root"
                    return os.path.isfile(path)

                df = df[df.apply(path_exists, axis=1)]
            elif n.startswith("combine"):
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

    legend1 = ax1.legend(loc="upper right")

    if len(csv_file_dirs) > 1:
        line_solid = Line2D([0], [0], color='black', linestyle='-', label="no BB-lite")#label='no bin by bin stat.')
        line_dashed = Line2D([0], [0], color='black', linestyle='--', label="BB-lite")#label='bin by bin stat.')

        # Second legend
        # ax1.legend(handles=[line_solid, line_dashed], loc='lower right')
        ax1.legend(handles=[line_solid, line_dashed], loc='upper left')

    ax1.add_artist(legend1)

    ax1.text(0.97, 0.02, "$N^\mathrm{bins} = "+f"{nBins}$",  transform=ax1.transAxes, va="bottom", ha="right")

    outdir = args.output
    outfile = f"model_{key}_nBins{nBins}"
    if args.postfix:
        outfile += f"_{args.postfix}"

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
        type=float,
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
    args = parser.parse_args()

    for key in args.keys:
        for nBins in (10, 100, 1000, 10000, 100000):
            create_plot(args.inputDirs, nBins, key, args)

if __name__ == "__main__":
    main()
