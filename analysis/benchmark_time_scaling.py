import pandas as pd
import matplotlib.pyplot as plt
import sys
import argparse
import os
from wums import output_tools, plot_tools  # isort: skip
from matplotlib.lines import Line2D
import glob
import rabbit.io_tools
import numpy as np
import pandas as pd

colors = {
    "Loss": "tab:blue",
    "HVP": "tab:red",
    "EDM": "tab:orange",
    "Minimizer": "tab:green",
    "Overhead": "tab:purple",
    "Copy": "tab:brown",
    "Init": "tab:pink",
    "Total": "tab:grey"
}

def create_plot(
        files, 
        args,
        xlabel="Number of systematics",
        ylabel="Time in second",
        xlim=(0.8, 16_384),
        ):
    
    outdir = output_tools.make_plot_dir(args.outpath, eoscp=False)

    mem_r = []
    mem_v = []
    is_jax = []

    n_systs = []
    n_bins = []

    n_grad = []
    n_hvp = []
    t_grad = []
    t_hvp = []
    t_edm = []
    t_min = []
    t_init = []
    t_tot = []

    t_c_l1 = []
    t_c_l2 = []
    t_c_h1 = []
    t_c_h2 = []
    for i, fitresult_file in enumerate(files):
        try:
            fitresult = rabbit.io_tools.get_fitresult(fitresult_file)
        except:
            continue
        if fitresult is None:
            continue
        if "time_init" not in fitresult.keys():
            continue

        filename = fitresult_file.split("/")[-2]
        parts = filename.split("_")

        is_jax.append("jax" in fitresult_file.split("/")[-1])

        n_systs.append(int(parts[-1].split("nSysts")[-1]))
        n_bins.append(int(parts[-2].split("nBins")[-1]))

        n_grad.append(fitresult["n_grad"])
        n_hvp.append(fitresult["n_hvp"])
        t_grad.append(fitresult["time_grad"])
        t_hvp.append(fitresult["time_hvp"])
        t_edm.append(fitresult["time_edm"])
        t_min.append(fitresult["time_minimizer"])
        t_init.append(fitresult["time_init"])
        t_tot.append(fitresult["time_total"])

        t_c_l1.append(fitresult["time_grad_copy_on"])
        t_c_l2.append(fitresult["time_grad_copy_off"])
        t_c_h1.append(fitresult["time_hvp_copy_on"])
        t_c_h2.append(fitresult["time_hvp_copy_off"])

        memlog_file = fitresult_file.replace("fitresults", "memlog").replace(".hdf5",".txt")

        if os.path.isfile(memlog_file):
        
            df = pd.read_csv(memlog_file, delimiter="\s+", skiprows=1, names=["time","CPU", "Real", "Virtual"])
            mem_r.append(df["Real"].max()/1000.)
            mem_v.append(df["Virtual"].max()/1000.)
        else:
            mem_r.append(0)
            mem_v.append(0)

    mem_r = np.array(mem_r)
    mem_v = np.array(mem_v)

    is_jax = np.array(is_jax)
    n_systs = np.array(n_systs)
    n_bins = np.array(n_bins)

    n_grad = np.array(n_grad)
    n_hvp = np.array(n_hvp)
    t_grad = np.array(t_grad)
    t_hvp = np.array(t_hvp)
    t_edm = np.array(t_edm)
    t_min = np.array(t_min)
    t_init = np.array(t_init)
    t_tot = np.array(t_tot)

    t_c_l1 = np.array(t_c_l1)
    t_c_l2 = np.array(t_c_l2)
    t_c_h1 = np.array(t_c_h1)
    t_c_h2 = np.array(t_c_h2)

    for nb in list(set(n_bins)):
        if args.extraText is not None:
            if args.extraTextCols==2:
                text = args.extraText+",\n"+r" $N^\mathrm{bins} = "+f"{nb}$"
            else:
                text = args.extraText+r", $N^\mathrm{bins} = "+f"{nb}$"
        else:
            text = r"$N^\mathrm{bins} = "+f"{nb}$"

        idxs = np.where(n_bins == nb)
        sort_idx = np.argsort(n_systs[idxs])

        m_r = mem_r[idxs][sort_idx]
        m_v = mem_v[idxs][sort_idx]

        is_j = is_jax[idxs][sort_idx]

        n_s = n_systs[idxs][sort_idx]
        n_g = n_grad[idxs][sort_idx]
        n_h = n_hvp[idxs][sort_idx]
        t_g = t_grad[idxs][sort_idx]
        t_h = t_hvp[idxs][sort_idx]
        t_e = t_edm[idxs][sort_idx]
        t_m = t_min[idxs][sort_idx]
        t_i = t_init[idxs][sort_idx]
        t_t = t_tot[idxs][sort_idx]

        t_l1 = t_c_l1[idxs][sort_idx]
        t_l2 = t_c_l2[idxs][sort_idx]
        t_h1 = t_c_h1[idxs][sort_idx]
        t_h2 = t_c_h2[idxs][sort_idx]

        # copy
        t_c = t_l1 + t_l2 + t_h1 + t_h2
        # overhead
        t_o = t_t - t_m - t_e - t_i
        # minimizer
        t_m = t_m - t_h - t_g - t_c

        # if any(t_t - (t_m + t_o + t_c + t_e + t_h + t_g) != 0):
        #     print("Breakdown not correct")

        xlim = (0.78, max(n_s)*1.3)

        # Create the plot
        fig, ax1 = plot_tools.figure([0,1], xlabel, ylabel, xlim=xlim, ylim=args.ylim, logx=True, logx_base=2)

        ax1.plot(n_s, n_g, label="Loss & Grad", color=colors["Loss"])
        ax1.plot(n_s, n_h, label="EDM", color=colors["EDM"])

        ax1.legend(loc="upper left")
        ax1.text(0.97, 0.02, text,  transform=ax1.transAxes, va="bottom", ha="right")


        outfile = f"nBins{nb}_calls"
        if args.postfix:
            outfile += f"_{args.postfix}"
        if args.subtitle == "Preliminary":
            outfile += "_preliminary"
        plot_tools.save_pdf_and_png(outdir, outfile)

        output_tools.write_index_and_log(
            outdir,
            outfile,
            args=args,
        )

        fig, ax1 = plot_tools.figure([0,1], xlabel, ylabel, xlim=xlim, ylim=(0.06, 500), logy=True, logx=True, logx_base=2)

        custom_lines = [
            Line2D([0], [0], color='black', linestyle='-', label='TF'),
            Line2D([0], [0], color='black', linestyle='--', label='JAX'),
        ]
        legend2 = ax1.legend(handles=custom_lines, loc=args.secondColPos)
        ax1.add_artist(legend2)

        for j in (0, 1):

            for y, l, c in (
                (t_t, "Total", "Total"),
                (t_i, "Init", "Init"),
                (t_h, "HVP" if not args.hess else "Hess", "HVP"),
                (t_g, "Loss & Grad", "Loss"),
                (t_e, "EDM", "EDM"),
                (t_m, "Minimizer", "Minimizer"),
                (t_o, "Overhead", "Overhead"),
                (t_c, "Copy", "Copy"),
            ):
                if max(y) < 1:
                    continue
                ax1.plot(n_s[is_j==j], y[is_j==j], label=l if j==0 else None, color=colors[c], linestyle="solid" if j==0 else "dashed")

            # ax1.plot(n_s, t_g, label="Loss & Grad", color=colors["Loss"])
            # ax1.plot(n_s, t_h, label="HVP" if not args.hess else "Hess", color=colors["HVP"])
            # ax1.plot(n_s, t_e, label="EDM", color=colors["EDM"])
            # ax1.plot(n_s, t_m, label="Minimizer", color=colors["Minimizer"])
            # ax1.plot(n_s, t_o, label="Overhead", color=colors["Overhead"])
            # ax1.plot(n_s, t_c, label="Copy", color=colors["Copy"])
            # ax1.plot(n_s, t_i, label="Init", color=colors["Init"])
            # ax1.plot(n_s, t_t, label="Total", color=colors["Total"])


        # ax1.plot(n_s, t_l1, label="Copy loss & grad on GPU")
        # ax1.plot(n_s, t_l2, label="Copy loss & grad on CPU")
        # ax1.plot(n_s, t_h1, label="Copy HVP on GPU")
        # ax1.plot(n_s, t_h2, label="Copy HVP on CPU")


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

        ax1.legend(loc="upper left", ncols=args.legCols)
        ax1.text(0.97, 0.02, text,  transform=ax1.transAxes, va="bottom", ha="right")

        outfile = f"nBins{nb}_time"
        if args.postfix:
            outfile += f"_{args.postfix}"
        if args.subtitle == "Preliminary":
            outfile += "_preliminary"
        plot_tools.save_pdf_and_png(outdir, outfile)

        output_tools.write_index_and_log(
            outdir,
            outfile,
            args=args,
        )

        fig, ax1 = plot_tools.figure([0,1], xlabel, ylabel, xlim=xlim, ylim=(0.01, 20), logy=True, logx=True, logx_base=2)

        ax1.plot(n_s, t_g/n_g, label="Loss & Grad", color=colors["Loss"])
        ax1.plot(n_s, t_h/n_h, label="HVP" if not args.hess else "Hess", color=colors["HVP"])
        ax1.plot(n_s, t_e, label="EDM", color=colors["EDM"])

        ax1.legend(loc="upper left")
        ax1.text(0.97, 0.02, text,  transform=ax1.transAxes, va="bottom", ha="right")

        outfile = f"nBins{nb}_time_per_call"
        if args.postfix:
            outfile += f"_{args.postfix}"
        if args.subtitle == "Preliminary":
            outfile += "_preliminary"
        plot_tools.save_pdf_and_png(outdir, outfile)

        output_tools.write_index_and_log(
            outdir,
            outfile,
            args=args,
        )

        fig, ax1 = plot_tools.figure([0,1], xlabel, "Max memory [GB]", xlim=xlim, ylim=(1, 2*max(m_v)), logy=True, logx=True, logx_base=2)

        custom_lines = [
            Line2D([0], [0], color='black', linestyle='-', label='TF'),
            Line2D([0], [0], color='black', linestyle='--', label='JAX'),
        ]
        legend2 = ax1.legend(handles=custom_lines, loc='upper center')
        ax1.add_artist(legend2)

        for j in (0, 1):
            ax1.plot(n_s[is_j==j], m_r[is_j==j], label="Real" if j==0 else None, color="red", linestyle="solid" if j==0 else "dashed")
            ax1.plot(n_s[is_j==j], m_v[is_j==j], label="Virtual" if j==0 else None, color="orange", linestyle="solid" if j==0 else "dashed")

        ax1.legend(loc="upper left")
        ax1.text(0.97, 0.02, text,  transform=ax1.transAxes, va="bottom", ha="right")

        outfile = f"nBins{nb}_mem"
        if args.postfix:
            outfile += f"_{args.postfix}"
        if args.subtitle == "Preliminary":
            outfile += "_preliminary"
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
        "--input",
        type=str,
        nargs="*",
        default="results",
        help="Input directory to fitresult files",
    )
    parser.add_argument(
        "-o",
        "--outpath",
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
    parser.add_argument(
        "--extraText",
        type=str,
        default="",
        help="Extra text to be shown in plot",
    )
    parser.add_argument(
        "--hess",
        action="store_true",
        help="If hessian was used in loss",
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
        "--legCols", type=int, default=1, help="Number of columns in legend"
    )
    parser.add_argument(
        "--secondColPos", type=str, default="upper center", help="Location of seonc legend"
    )
    parser.add_argument(
        "--extraTextCols", type=int, default=1, help="Extra text columns"
    )
    args = parser.parse_args()
    
    create_plot(args.input, args)

if __name__ == "__main__":
    main()