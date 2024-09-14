import pandas as pd
import os
import sys
from FindData import FindData
import seaborn as sns
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


class Plot_Findings:
    def __init__(self,path, output_folder):
        self.path = path
        self.ConfMatrixObj = FindData(path,"Compare_conf_matrix.csv", "deep")
        self.GSB_CM = self.ConfMatrixObj.GSB_df
        self.DET_CM = self.ConfMatrixObj.DET_df
        self.NR_DMRs = FindData(path,"Compare_Nr_DMRs.csv", "general")
        self.NR_DMRs_df = self.NR_DMRs.GEN_df
        self.GSB_CM = self.GSB_CM.reset_index()
        self.DET_CM = self.DET_CM.reset_index()
        self.NR_DMRs_df = self.NR_DMRs_df.reset_index()
        self.PerformanceObj = FindData(path,"Compare_performance.csv", "deep")
        self.GSB_P = self.PerformanceObj.GSB_df
        self.DET_P = self.PerformanceObj.DET_df
        self.GSB_P = self.GSB_P.reset_index()
        self.DET_P = self.DET_P.reset_index()
        self.output_folder = output_folder



    def plot_NrDMRs(self):
        print(self.NR_DMRs_df)
        df = self.NR_DMRs_df[["DiMmer","dmrseq", "DMRcate", "BSmooth", "Gold-Standard", "NR_DMR"]]
        # Ensure 'NR_DMR' is treated as a categorical variable for x-axis
        df['NR_DMR'] = pd.Categorical(df['NR_DMR'], categories=[1500, 3000, 4500, 6000], ordered=True)
        print(df)

        # Pivot the data for plotting
        pivot_df = df.melt(id_vars='NR_DMR', value_vars=['DiMmer', 'dmrseq', 'DMRcate', 'BSmooth', 'Gold-Standard'],
                             var_name='Method', value_name='Detected_DMRs')
        print(pivot_df)

        plt.figure(figsize=(10, 6))
        sns.lineplot(data=pivot_df, x='NR_DMR', y='Detected_DMRs', hue='Method', marker='o')

        plt.xlabel('Number of Simulated DMRs')
        plt.ylabel('Number of Detected DMRs')
        plt.title('Number of Detected DMRs by Method')
        plt.grid(True)
        plt.xticks([1500, 3000, 4500, 6000])
        plt.legend(title='Method')
        #plt.show()
        savepath = os.path.join(self.output_folder,'NR_DMRs.png')
        plt.savefig(savepath, dpi=300)

    def plot_CM(self, CM_df):
        CM_df['NR_DMR'] = pd.Categorical(CM_df['NR_DMR'], categories=[1500, 3000, 4500, 6000], ordered=True)
        CM_df_long = pd.melt(CM_df, id_vars=["type", "PERC", "NR_DMR"],
                          value_vars=["DiMmer", "dmrseq", "DMRcate", "BSmooth"],
                          var_name="Method", value_name="Value")
        # Set the plot style
        sns.set(style="whitegrid")

        # Initialize a grid of plots with seaborn's FacetGrid
        g = sns.FacetGrid(CM_df_long,row="PERC", col="type",   hue="Method", margin_titles=True, height=4)

        # Map a lineplot onto the grid for each tool
        g.map(sns.lineplot, "NR_DMR", "Value",marker='o')

        # Add legends
        g.add_legend()

        # Adjust titles
        g.set_axis_labels("Number of Simulated DMRs", "Values")
        g.set_titles(row_template="Higher than {row_name} Overlap Fraction", col_template="{col_name} Count")

        # Show the plot
        #plt.show()
        savepath = os.path.join(self.output_folder,'ConfMatrix.png')
        plt.savefig(savepath, dpi=300)

    def plot_Perf(self, Perf_df):
        Perf_df['NR_DMR'] = pd.Categorical(Perf_df['NR_DMR'], categories=[1500, 3000, 4500, 6000], ordered=True)
        Perf_df_long = pd.melt(Perf_df, id_vars=["type", "PERC", "NR_DMR"],
                             value_vars=["DiMmer", "dmrseq", "DMRcate", "BSmooth"],
                             var_name="Method", value_name="Value")
        # Set the plot style
        sns.set(style="whitegrid")

        # Initialize a grid of plots with seaborn's FacetGrid
        g = sns.FacetGrid(Perf_df_long, row="PERC", col="type", hue="Method", margin_titles=True, height=4)

        # Map a lineplot onto the grid for each tool
        g.map(sns.lineplot, "NR_DMR", "Value",marker='o')

        # Add legends
        g.add_legend()

        # Adjust titles
        g.set_axis_labels("Number of Simulated DMRs", "Values")
        g.set_titles(row_template="Higher than {row_name} Overlap Fraction", col_template="{col_name}")

        # Show the plot
        #plt.show()
        savepath = os.path.join(self.output_folder,'Performance.png')
        plt.savefig(savepath, dpi=300)




if __name__ == '__main__':
    setup = Plot_Findings('/Users/leonoraka/Desktop/MA-project/Snakemake/Analysis/05_Evaluation/', '/Users/leonoraka/Desktop/MA-project/Snakemake/Analysis/07_Overall/')
    setup.plot_NrDMRs()
    setup.plot_CM(setup.GSB_CM)
    setup.plot_Perf(setup.GSB_P)