import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from statsmodels.stats.proportion import proportions_ztest

# make the default font size point 7
plt.rcParams.update({"font.size": 7})


# define function for calculating sensitivity, specificity, PPV and NPV
# calculate sensitivity: number of phenotypically resistant samples that are found to have a resistant mutation / total number of phenotypically resistant samples
# calculate specificity: number of phenotypically sensitive samples that do not have a resistance mutation / total number of phenotypically sensitive samples
# calculate PPV: number of samples with a resistance mutation that are phenotypically resistant / total number of samples with a resistance mutation
# calculate NPV: number of samples without a resistance mutation that are phenotypically sensitive / total number of samples without a resistance mutation
def calculate_statistics(df):
    sensitivity = (
        df[(df["IS_RESISTANT"]) & (df["HAS_RESISTANT_MUTATION"])].shape[0]
        / df[df["IS_RESISTANT"]].shape[0]
    )
    specificity = (
        df[(~df["IS_RESISTANT"]) & (~df["HAS_RESISTANT_MUTATION"])].shape[0]
        / df[~df["IS_RESISTANT"]].shape[0]
    )
    PPV = (
        df[(df["HAS_RESISTANT_MUTATION"]) & (df["IS_RESISTANT"])].shape[0]
        / df[df["HAS_RESISTANT_MUTATION"] == True].shape[0]
    )
    NPV = (
        df[(~df["HAS_RESISTANT_MUTATION"]) & (~df["IS_RESISTANT"])].shape[0]
        / df[df["HAS_RESISTANT_MUTATION"] == False].shape[0]
    )
    return sensitivity, specificity, PPV, NPV


# plot sensitivity and specificity for no FRS cutoff and for FRS cutoff >=0.9 (major allele definition) as a histogram
def create_stats_df(sen_1, sen_2, spe_1, spe_2):
    stats = pd.DataFrame(
        {
            "Metric": ["Sensitivity", "Specificity"],
            "FRS cutoff >=0.9 (major allele definition)": [sen_1, spe_1],
            "No FRS cutoff": [sen_2, spe_2],
        }
    )
    return stats


def plot_stats_histogram(
    stats,
    significant=True,
    save_figure=False,
    save_name="sensitivity_specificity_histogram",
    colours=[
        "#ff7f0e",
        "#1f77b4",
    ],
    metric_labels=["senstivity", "specificity"],
    legend_label_1="FRS cutoff >=0.9",
    legen_label_2="No FRS cutoff",
):
    # Adjusting the plot to group bars by criteria and display percentage inside the bar
    fig, ax = plt.subplots(figsize=(2.5, 2.75))

    # Plotting bars next to each other
    x = [0, 0.7]  # positions for FRS cutoff >=0.9 (major allele definition)
    bar_width = 0.25

    bars1 = ax.bar(
        x,
        stats["FRS cutoff >=0.9 (major allele definition)"],
        width=bar_width,
        color=colours,
        edgecolor=None,
        linewidth=1.2,
        label=legend_label_1,
    )
    bars2 = ax.bar(
        [p + bar_width for p in x],
        stats["No FRS cutoff"],
        width=bar_width,
        color=["white", "white"],
        edgecolor=colours,
        linewidth=1.2,
        label=legen_label_2,
    )

    # bars1[0].set_color('#ff7f0e')

    # Adding percentage values inside the bars
    for bars in [bars1, bars2]:
        for bar in bars:
            yval = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                yval + 0.02,
                f"{yval:.1%}",
                ha="center",
                va="bottom",
                fontsize=7,
                color="k",
            )

    # Customize the plot
    # Set y-axis to display values as percentages without the percent sign
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: "{:.0f}".format(y * 100)))

    # ax.set_title('Sensitivity and Specificity for Rifampicin Resistance', fontsize=18)
    # ax.set_ylabel('%', fontsize=7)
    ax.set_ylim(0, 1.1)

    # Setting custom x-axis labels
    ax.set_xticks([r + bar_width / 2 for r in x])
    ax.set_xticklabels(metric_labels, fontsize=7)

    # Increase axis fontsizes
    ax.tick_params(axis="both", which="major", labelsize=7)

    # Adding the legend
    # ax.legend(loc='lower right', fontsize=12)

    # remove unnecessary spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # add significance indicator with horizontal bar between bars of sensitivity
    if significant:
        ax.plot([0, 0], [1.05, 1.07], color="black", lw=2)
        ax.plot([0.2, 0.2], [1.05, 1.07], color="black", lw=2)
        ax.plot([0, 0.2], [1.07, 1.07], color="black", lw=2)
        ax.text(0.1, 1.08, "*", fontsize=20, color="black", ha="center", va="center")

    plt.tight_layout()
    plt.show()

    # save figure if print_stats is True
    if save_figure:
        fig.savefig(f"pdf/{save_name}.pdf", bbox_inches="tight", transparent=True)


def plot_performance_frs(
    df,
    save_figure=False,
    save_name="fig-sensitivity_specificity_FRS_cutoff",
    metrics=["sensitivity", "specificity"],
    colours=["#ff7f0e", "#1f77b4"],
):

    fig, ax1 = plt.subplots(figsize=(5, 3))

    # plot sensitivity and specificity using one y axis
    ax1.plot(df.min_FRS, df[metrics[0]], label=metrics[0], color=colours[0])
    ax1.plot(df.min_FRS, df[metrics[1]], label=metrics[1], color=colours[1])
    # plt.legend(fontsize=15)
    ax1.set_xlabel("fraction of reads (FRS) supporting RAV (%)")
    # plt.ylabel('%', fontsize=15)
    # plt.title('Sensitivity and Specificity for different FRS cutoffs')

    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.set_ylim(0.94, 0.99)
    ax1.set_yticks(np.arange(0.94, 0.99, 0.01))
    ax1.set_xlim(0, 1)
    # increase tick label sizes
    # plt.xticks(fontsize=15)
    # plt.yticks(fontsize=15)

    # Set axes to display values as percentages without the percent sign
    plt.gca().yaxis.set_major_formatter(
        FuncFormatter(lambda y, _: "{:.0f}".format(y * 100))
    )
    plt.gca().xaxis.set_major_formatter(
        FuncFormatter(lambda x, _: "{:.0f}".format(x * 100))
    )

    if save_figure:
        fig.savefig(
            f"pdf/{save_name}.pdf",
            bbox_inches="tight",
            transparent=True,
        )

    plt.show()


def plot_histogram(
    df,
    save_figure=False,
    save_name="sensitivity_specificity_histogram",
    colours=[
        "#ff7f0e",
        "#1f77b4",
    ],
    metrics=None,
):
    # Adjusting the plot to group bars by criteria and display percentage inside the bar
    fig, ax = plt.subplots(figsize=(2.5, 2.75))

    # Plotting bars next to each other
    x = [0, 0.7]  # positions for FRS cutoff >=0.9 (major allele definition)
    bar_width = 0.25
    b1 = ax.bar(0, df.iloc[0][metrics[0]], width=bar_width, color=colours[0])
    b2 = ax.bar(0.7, df.iloc[0][metrics[1]], width=bar_width, color=colours[1])
    b3 = ax.bar(
        0.25,
        df.iloc[1][metrics[0]],
        width=bar_width,
        color="white",
        edgecolor=colours[0],
    )
    b4 = ax.bar(
        0.95,
        df.iloc[1][metrics[1]],
        width=bar_width,
        color="white",
        edgecolor=colours[1],
    )

    # Adding percentage values inside the bars
    for bar in [b1, b2, b3, b4]:
        yval = bar[0].get_height()
        ax.text(
            bar[0].get_x() + bar[0].get_width() / 2,
            yval + 0.02,
            f"{yval:.1%}",
            ha="center",
            va="bottom",
            fontsize=7,
            color="k",
        )

    # Customize the plot
    # Set y-axis to display values as percentages without the percent sign
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: "{:.0f}".format(y * 100)))

    # ax.set_title('Sensitivity and Specificity for Rifampicin Resistance', fontsize=18)
    # ax.set_ylabel('%', fontsize=7)
    ax.set_ylim(0, 1.1)

    # Setting custom x-axis labels
    ax.set_xticks([r + bar_width / 2 for r in x])
    ax.set_xticklabels(metrics, fontsize=7)

    # Increase axis fontsizes
    ax.tick_params(axis="both", which="major", labelsize=7)

    # Adding the legend
    # ax.legend(loc='lower right', fontsize=12)

    # remove unnecessary spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    a = df.iloc[0]
    b = df.iloc[1]
    for i in range(len(metrics)):
        metric = metrics[i]
        if metric == "sensitivity":
            counts = [int(a.TP), int(b.TP)]
            nobs = [int(a.P), int(b.P)]
        elif metric == "specificity":
            counts = [int(a.TN), int(b.TN)]
            nobs = [int(a.N), int(b.N)]
        elif metric == "PPV":
            counts = [int(a.TP), int(b.TP)]
            nobs = [int(a.TP) + int(a.FP), int(b.TP) + int(b.FP)]
        elif metric == "NPV":
            counts = [int(a.TN), int(b.TN)]
            nobs = [int(a.TN) + int(a.FN), int(b.TN) + int(b.FN)]

        z, p = proportions_ztest(counts, nobs)
        if z < 0.05:
            print(metric, z, p)
            x_value = x[i]
            ax.plot([x_value, x_value], [1.07, 1.09], color="black", lw=2)
            ax.plot([x_value + 0.2, x_value + 0.2], [1.07, 1.09], color="black", lw=2)
            ax.plot([x_value, x_value + 0.2], [1.09, 1.09], color="black", lw=2)
            ax.text(
                x_value + 0.1,
                1.1,
                "*",
                fontsize=20,
                color="black",
                ha="center",
                va="center",
            )

    plt.tight_layout()
    plt.show()

    # save figure if print_stats is True
    if save_figure:
        fig.savefig(f"pdf/{save_name}.pdf", bbox_inches="tight", transparent=True)


def plot_frs_histogram(df, savefig=False, save_name="hist-FRS-RAV"):

    axes = df.FRS.hist(bins=20, color="grey")

    axes.set_xlabel("Fraction of reads supporting RAV (%)")
    axes.spines["top"].set_visible(False)
    axes.spines["right"].set_visible(False)
    axes.set_xlim(0, 1)
    axes.grid(False)
    fig = axes.get_figure()
    fig.set_size_inches(3, 2)
    if savefig:
        fig.savefig(f"pdf/{save_name}.pdf", bbox_inches="tight", transparent=True)
