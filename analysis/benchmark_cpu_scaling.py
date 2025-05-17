import pandas as pd
import matplotlib.pyplot as plt
import sys
import argparse
import os
from wums import logging, output_tools, plot_tools  # isort: skip

def create_plot(
        csv_file_dir, 
        args,
        xlabel="Number of CPU threads",
        ylabel="Time in second",
        ):
    
    # Create the plot
    fig, ax1 = plot_tools.figure([], xlabel, ylabel, xlim=(0.8, 1000), logx=True)

    # Load the CSV files
    for l, c, m, s in (
        ("combinetf1", "blue", "o", "-"),
        ("combinetf2", "red", "x", "-"),
    ):
        csv_file = f"{csv_file_dir}/timing_cpu_scaling_{l}.csv"

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

        ax1.plot(df["nCPU"].values, df["time"].values, label=l, marker=m, linestyle=s, color=c)

    ax1.legend()

    outdir = args.output
    outfile = "cpu_time"

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
        "--inputDir",
        type=str,
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
    args = parser.parse_args()
    
    create_plot(args.inputDir, args)

if __name__ == "__main__":
    main()