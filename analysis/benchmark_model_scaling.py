import pandas as pd
import matplotlib.pyplot as plt
import sys
import argparse
import os
from wums import logging, output_tools, plot_tools  # isort: skip
from matplotlib.lines import Line2D

def create_plot(
        csv_file_dirs, 
        nBins,
        args,
        xlabel="Number of systematics",
        ylabel="Time in second",
        ):
    
    # Create the plot
    fig, ax1 = plot_tools.figure([0,1], xlabel, ylabel, xlim=args.xlim, ylim=args.ylim, logy=True, logx=True, logx_base=2)

    # Load the CSV files
    linestyles = ["-", "--"]
    for n, l, c, m in (
        ("combine", "Combine", "green", "."),
        ("combinetf1", "CombineTF", "blue", "o"),
        # ("combinetf2_singularity", "CombineTF 2 (sing.)", "orange", "*"),
        ("combinetf2", "CombineTF 2 (virt. env.)", "red", "x"),
    ):
        for i, csv_file_dir in enumerate(csv_file_dirs):
            csv_file = f"{csv_file_dir}/timing_model_scaling_{n}.csv"

            if not os.path.isfile(csv_file):
                continue
            
            df = pd.read_csv(csv_file)
            print(f"Successfully loaded {csv_file}")
            print(f"Columns found: {', '.join(df.columns)}")

            df = df.loc[df["nBins"]==nBins]

            # Ensure required columns exist
            required_cols = ['nSyst', 'fit']
            for col in required_cols:
                if col not in df.columns:
                    print(f"Error: Required column '{col}' not found in CSV. Available columns: {df.columns.tolist()}")
                    sys.exit(1)
            
            # Convert columns to appropriate types
            df['nSyst'] = pd.to_numeric(df['nSyst'], errors='coerce')
            df['fit'] = pd.to_numeric(df['fit'], errors='coerce')
            
            # Remove any rows with missing values after conversion
            orig_len = len(df)
            df = df.dropna(subset=['nSyst', 'fit'])
            if len(df) < orig_len:
                print(f"Warning: Removed {orig_len - len(df)} rows with non-numeric values")
            
            # Group by nSyst (although seaborn will handle this for us in the plot)
            cpu_groups = df.groupby('nSyst')
            print(f"Number of unique CPU values: {len(cpu_groups)}")
            for cpu, group in cpu_groups:
                print(f"nSyst={cpu}: {len(group)} measurements, Avg time: {group['fit'].mean():.4f}")

            ax1.plot(df["nSyst"].values, df["fit"].values, label=l if i==0 else None, marker=m, linestyle=linestyles[i], color=c)

    legend1 = ax1.legend(loc="upper left")

    if len(csv_file_dirs) > 1:
        line_solid = Line2D([0], [0], color='black', linestyle='-', label='2 EPYC 9965')
        line_dashed = Line2D([0], [0], color='black', linestyle='--', label='2 EPYC 9654')

        # Second legend
        ax1.legend(handles=[line_solid, line_dashed], loc='upper right')

    ax1.add_artist(legend1)

    ax1.text(0.95, 0.9, f"N(bins) = {nBins}",  transform=ax1.transAxes, ha="right")

    outdir = args.output
    outfile = f"model_time_nBins{nBins}"
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

    for nBins in (10, 100, 1000):
        create_plot(args.inputDirs, nBins, args)

if __name__ == "__main__":
    main()