import numpy as np

from Evaluation.ReadResults import ReadResults
import os.path
import os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

class IndividualAnalysis:
    def __init__(self, read_results, out_folder):
        self.ReadResults: ReadResults = read_results
        self.tab = self.ReadResults.processed_table
        self.method = self.ReadResults.method_type
        self.score = self.ReadResults.score_type
        self.score_threshold = self.ReadResults.score_threshold
        self.out_folder = os.path.join(out_folder, "Individual_Analysis_" + self.method)
        os.makedirs(self.out_folder, exist_ok=True)
        self.get_lengths()
        if self.method != "Mouse_True" and self.method != "Mouse_Known":
            if self.method != "Simulated":
                self.score_histograms()
            self.get_Nr_CpGs()
        self.Nr_DMRs = self.get_Nr_DMRs()


    def get_Nr_CpGs(self):
        plt.figure(figsize=(10, 6))
        sns.set(style="whitegrid")

        # Create the histogram
        sns.histplot(self.tab['num_CpGs'], bins=10, kde=False, color='skyblue')

        # Add title and subtitle
        plt.title('Distribution of CpGs in ' + self.method, fontsize=16)
        plt.suptitle('Histogram of number of CpGs across found DMRs', fontsize=12)

        # Label the axes
        plt.xlabel('Number of CpGs', fontsize=14)
        plt.ylabel('Frequency', fontsize=14)

        # Improve layout
        plt.tight_layout()

        # Show & save the plot
        #plt.show()
        savepath = os.path.join(self.out_folder, self.method +'_CpGs_histogram.png')
        plt.savefig(savepath, dpi=300)

    def get_lengths(self):
        plt.figure(figsize=(8, 5))
        sns.set(style="ticks", palette="muted", font_scale=1.2)

        # Create the histogram
        #print("Bug: ")
        #print(np.array(self.tab['width'],dtype=float))
        sns.histplot(np.array(self.tab['width'],dtype=float), bins=10, kde=True, color='gray', edgecolor='black')

        # Add title and subtitle
        plt.title('Distribution of DMR length in ' + self.method, fontsize=16, weight='bold')
        plt.suptitle('Histogram of length across DMRs', fontsize=12, y=0.93)

        # Label the axes
        plt.xlabel('Genomic length', fontsize=14)
        plt.ylabel('Frequency', fontsize=14)

        # Remove top and right spines for a cleaner look
        sns.despine()

        # Improve layout
        plt.tight_layout()

        # Save the plot
        savepath = os.path.join(self.out_folder, self.method +'_length_histogram.png')
        plt.savefig(savepath, dpi=300, bbox_inches='tight')  # Save as a PNG file with 300 DPI


    def score_histograms(self):
        if self.score == 'p-value' or self.score == 'fdr':
            # Set up the figure size and style
            plt.figure(figsize=(8, 5))
            sns.set(style="whitegrid", font_scale=1.2)

            # Create the histogram
            sns.histplot(self.tab['score'], bins=10, kde=False, color='dodgerblue', edgecolor='black')

            # Add title and subtitle
            plt.title('P-value Distribution of DMRs in ' + self.method, fontsize=16, weight='bold')
            plt.suptitle('Histogram of p-values (score) across DMRs', fontsize=12, y=0.93)

            # Label the axes
            plt.xlabel('P-value', fontsize=14)
            plt.ylabel('Frequency', fontsize=14)

            # Improve layout
            plt.tight_layout()

            # Save the plot
            savepath = os.path.join(self.out_folder, self.method + '_pvalue_histogram.png')
            plt.savefig(savepath, dpi=300, bbox_inches='tight')  # Save as a PNG file with 300 DPI
            self.p_value_ECDF_plot()
        elif self.score == 'meanDiff':
            plt.figure(figsize=(8, 5))
            sns.set(style="whitegrid", font_scale=1.2)

            # Create the histogram
            sns.histplot(self.tab['score'], bins=10, kde=False, color='mediumseagreen', edgecolor='black')

            # Add title and subtitle
            plt.title('Mean Difference Distribution of DMRs in ' + self.method, fontsize=16, weight='bold')
            plt.suptitle('Histogram of meanDiff across DMRs', fontsize=12, y=0.93)

            # Label the axes
            plt.xlabel('Mean Difference', fontsize=14)
            plt.ylabel('Frequency', fontsize=14)

            # Improve layout
            plt.tight_layout()

            # Save the plot
            savepath = os.path.join(self.out_folder, self.method + '_meanDiff_histogram.png')
            plt.savefig(savepath, dpi=300, bbox_inches='tight')  # Save as a PNG file with 300 DPI
        elif self.score == 'Fisher':
            # Set up the figure size and style
            plt.figure(figsize=(8, 5))
            sns.set(style="whitegrid", font_scale=1.2)

            # Create the histogram
            sns.histplot(self.tab['score'], bins=10, kde=False, color='mediumseagreen', edgecolor='black')

            # Add title and subtitle
            plt.title('Distribution of Fisher\'s Multiple Comparison Statistic of DMRs in ' + self.method, fontsize=16, weight='bold')
            plt.suptitle('Histogram of Fisher\'s Statistic across comparisons in DMRs', fontsize=12, y=0.93)

            # Label the axes
            plt.xlabel('Fisher\'s Statistic', fontsize=14)
            plt.ylabel('Frequency', fontsize=14)

            # Improve layout
            plt.tight_layout()

            # Save the plot
            savepath = os.path.join(self.out_folder, self.method + '_Fisher_histogram.png')
            plt.savefig(savepath, dpi=300,
                        bbox_inches='tight')  # Save as a PNG file with 300 DPI
            self.Fisher_CDF_plot()


    def get_Nr_DMRs(self):#
        #TODO: significant number of DMRs #FRAGE: output dimmer already signifikant??unique nicht gecheckt
        nr_DMRs = len(self.tab)
        return nr_DMRs

    def p_value_ECDF_plot(self):
        # Set up the figure size and style
        plt.figure(figsize=(8, 5))
        sns.set(style="ticks", font_scale=1.2)

        # Create the ECDF plot
        sns.ecdfplot(self.tab['score'], color='darkred')

        # Add title and subtitle
        plt.title('Empirical Cumulative Distribution of P-values (' + self.method + ")", fontsize=16, weight='bold')
        plt.suptitle('ECDF of p-values across DMRs', fontsize=12, y=0.93)

        # Label the axes
        plt.xlabel('P-value', fontsize=14)
        plt.ylabel('Proportion', fontsize=14)

        # Remove top and right spines for a cleaner look
        sns.despine()

        # Improve layout
        plt.tight_layout()

        # Save the plot
        savepath = os.path.join(self.out_folder, self.method + '_Score_ECDF.png')
        plt.savefig(savepath, dpi=300, bbox_inches='tight')  # Save as a PNG file with 300 DPI


    def Fisher_CDF_plot(self): #bringt das was
        # Create the ECDF plot
        sns.ecdfplot(self.tab['score'], color='darkblue')

        # Add title and subtitle
        plt.title('Cumulative Distribution of Fisher\'s Statistic in DMRs in ' + self.method, fontsize=16, weight='bold')
        plt.suptitle('ECDF of Fisher\'s Statistic across comparisons ', fontsize=12, y=0.93)

        # Label the axes
        plt.xlabel('Fisher\'s Statistic', fontsize=14)
        plt.ylabel('Proportion', fontsize=14)

        # Remove top and right spines for a cleaner look
        sns.despine()

        # Improve layout
        plt.tight_layout()

        # Save the plot
        savepath = os.path.join(self.out_folder, self.method + '_Score_Fisher_CDF.png')
        plt.savefig(savepath, dpi=300, bbox_inches='tight')  # Save as a PNG file with 300 DPI

##runtimes

##plots