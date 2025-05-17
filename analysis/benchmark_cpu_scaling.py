import pandas as pd
import matplotlib.pyplot as plt
import sys
import argparse
import os
from wums import logging, output_tools, plot_tools  # isort: skip
from matplotlib.lines import Line2D

def create_plot(
        csv_file_dirs, 
        args,
        xlabel="Number of CPU threads",
        ylabel="Time in second",
        ):
    
    # Create the plot
    fig, ax1 = plot_tools.figure([0,1], xlabel, ylabel, xlim=(0.8, 1000), ylim=args.ylim, logy=True, logx=True, logx_base=2)

    # Load the CSV files
    linestyles = ["-", "--"]
    for n, l, c, m in (
        ("combinetf1", "CombineTF", "blue", "o"),
        ("combinetf2_singularity", "CombineTF 2 (sing.)", "orange", "*"),
        ("combinetf2", "CombineTF 2 (virt. env.)", "red", "x"),
    ):
        for i, csv_file_dir in enumerate(csv_file_dirs):
            csv_file = f"{csv_file_dir}/timing_cpu_scaling_{n}.csv"

            if not os.path.isfile(csv_file):
                continue
            
            df = pd.read_csv(csv_file)
            print(f"Successfully loaded {csv_file}")
            print(f"Columns found: {', '.join(df.columns)}")

            # Ensure required columns exist
            required_cols = ['nCPU', 'time']
            for col in required_cols:
                if col not in df.columns:
                    print(f"Error: Required column '{col}' not found in CSV. Available columns: {df.columns.tolist()}")
                    sys.exit(1)
            
            # Convert columns to appropriate types
            df['nCPU'] = pd.to_numeric(df['nCPU'], errors='coerce')
            df['time'] = pd.to_numeric(df['time'], errors='coerce')
            
            # Remove any rows with missing values after conversion
            orig_len = len(df)
            df = df.dropna(subset=['nCPU', 'time'])
            if len(df) < orig_len:
                print(f"Warning: Removed {orig_len - len(df)} rows with non-numeric values")
            
            # Group by nCPU (although seaborn will handle this for us in the plot)
            cpu_groups = df.groupby('nCPU')
            print(f"Number of unique CPU values: {len(cpu_groups)}")
            for cpu, group in cpu_groups:
                print(f"nCPU={cpu}: {len(group)} measurements, Avg time: {group['time'].mean():.4f}")

            ax1.plot(df["nCPU"].values, df["time"].values, label=l if i==0 else None, marker=m, linestyle=linestyles[i], color=c)

    legend1 = ax1.legend(loc="lower left")

    if len(csv_file_dirs) > 1:
        line_solid = Line2D([0], [0], color='black', linestyle='-', label='2 EPYC 9965')
        line_dashed = Line2D([0], [0], color='black', linestyle='--', label='2 EPYC 9654')

        # Second legend
        ax1.legend(handles=[line_solid, line_dashed], loc='upper right')

    ax1.add_artist(legend1)

    outdir = args.output
    outfile = "cpu_time"
    if args.postfix:
        outfile += f"_{args.postfix}"

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
        "--ylim",
        type=float,
        nargs=2,
        help="Min and max values for y axis (if not specified, range set automatically)",
    )
    args = parser.parse_args()
    
    create_plot(args.inputDirs, args)

if __name__ == "__main__":
    main()