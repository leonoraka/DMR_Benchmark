#from Evaluation import *
from Evaluation.ReadResults import ReadResults
from Evaluation.Individual_Analysis import IndividualAnalysis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, roc_curve, auc
import os.path
import os
import seaborn as sns
from math import pi
from tqdm import tqdm
from scipy.stats import spearmanr, kendalltau

pd.options.mode.chained_assignment = None
class CompleteAnalysis:
    def __init__(self, input_dict, output_folder, DS):
        self.results_DiMmer = ReadResults(input_dict['DiMmer'],"DiMmer")
        self.results_dmrseq = ReadResults(input_dict['dmrseq'],"dmrseq")
        if DS == "sim":
            self.results_DMRcate = ReadResults(input_dict['DMRcate'],"DMRcate")
            self.results_BSmooth = ReadResults(input_dict['BSmooth'], "BSmooth")
            self.results_all = {"DiMmer": self.results_DiMmer, "dmrseq": self.results_dmrseq,
                            "DMRcate": self.results_DMRcate, "BSmooth": self.results_BSmooth}
            self.simulated = ReadResults(input_dict['Simulated'],"Simulated") #in ReadResults
            self.input_simulated = ReadResults(input_dict['input_simulated'],"input_simulated")
            self.mouse = False
        elif DS == "Mouse_True" or DS == "Mouse_Known":
            self.results_all = {"DiMmer": self.results_DiMmer, "dmrseq": self.results_dmrseq}
            self.simulated = ReadResults(input_dict['Mouse'], DS)  # in ReadResults
            self.mouse = True
        self.output_folder = output_folder
        self.analysis_folder = os.path.join(output_folder, "CompleteAnalysis")

        self.base_ref = None #["det";"gsb"]
        self.filter_perc = 0

        self.complete_analysis_folder = None

        if self.base_ref == "det":
            self.complete_analysis_folder = os.path.join(self.analysis_folder, "Detection-based" + str(self.filter_perc))
        elif self.base_ref == "gsb":
            self.complete_analysis_folder = os.path.join(self.analysis_folder, "Gold-Standard-Based" + str(self.filter_perc))

        if not os.path.exists(self.analysis_folder):
            os.makedirs(self.analysis_folder)
        self.individual_analysis_results = None


        self.TP = {}
        self.FP = {}
        self.FN = {}

        self.precision = {}
        self.recall = {}
        self.F1 = {}

        self.TP_CpG = {}
        self.FP_CpG = {}
        self.FN_CpG= {}
        self.TN_CpG = {}

        self.precision_CpG = {}
        self.recall_CpG = {}
        self.F1_CpG = {}
        self.accuracy_CpG = {}



    def compute_individual_analysis(self):
        self.individual_analysis_results = {method: IndividualAnalysis(self.results_all[method],self.output_folder) for method in self.results_all.keys()}
        if not self.mouse:
            self.individual_analysis_results['Gold-Standard'] = IndividualAnalysis(self.simulated,self.output_folder)


    def get_performance_all(self, mode, method):
        if mode == "DMR based":
            NEW = self.find_overlaps_optimized(self.results_all[method].processed_table,
                                                            self.simulated.processed_table)
            NEW_TP, NEW_FP = NEW["chosen_TP"], NEW["FP"]
            NEW_per_sim = NEW["TP_per_sim"]
            NEW_all_sim = NEW["all_sim"]
            NEW_FN = self.getFN(NEW_TP,self.simulated.processed_table)

            #self.plot_roc_curve(NEW_TP, NEW_FP, method)
            if self.base_ref == "det":
                self.TP[method] = len(NEW_TP)
            elif self.base_ref == "gsb":
                self.TP[method] = len(NEW_per_sim)
            self.FP[method] = len(NEW_FP)
            self.FN[method] = len(NEW_FN)
            self.precision[method] = self.TP[method] / (self.TP[method] + self.FP[method])
            self.recall[method] = self.TP[method] / (self.TP[method] + self.FN[method])
            self.F1[method] = (2 * self.precision[method] * self.recall[method]) / (self.precision[method] + self.recall[method])
            self.plot_roc_curve(NEW_TP, NEW_FP, method, mode)
            self.plot_precision_recall_curve(NEW_TP,NEW_FP, NEW_FN, method, mode)
            #self.compare_performance(mode)
            #self.compare_conf_matrix(mode)
        return {"all_overlaps" : NEW_TP, "overlap_per_sim" : NEW_per_sim, "all_sim" : NEW_all_sim}


      #return {"TP": NEW_TP, "FP": NEW_FP, "FN": NEW_FN}
    def analysis(self): #TODO write out tables
        print("STEP 1: Compute Individual Analysis")
        self.compute_individual_analysis()
        print("STEP 2: Plot Comparisons - Nr_DMRs, lengths, Scores")
        modes =["DMR based"]#,"CpG based"
        self.complete_analysis_folder = self.analysis_folder
        self.compare_Nr_DMRs()
        self.compare_lengths()
        self.compare_scores()
        #self.get_performance_all("CpG based", self.CpG_performance)
        print("STEP 3: EVALUATION")
        fractions_methods = {}
        base_refs = ["det","gsb"]
        percs = [0,0.1,0.2,0.5]
        if self.mouse:
            percs = [0,0.1]
        for b in base_refs:
            for p in percs:
                self.base_ref = b
                self.filter_perc = p
                if self.base_ref == "det":
                    self.complete_analysis_folder = os.path.join(self.analysis_folder,
                                                                 "Detection-based" + str(self.filter_perc))
                elif self.base_ref == "gsb":
                    self.complete_analysis_folder = os.path.join(self.analysis_folder,
                                                                 "Gold-Standard-Based" + str(self.filter_perc))
                if not os.path.exists(self.complete_analysis_folder):
                    os.makedirs(self.complete_analysis_folder)
                print("Evaluation in {} mode".format("DMR based"))
                print("Base ref: {}".format(self.base_ref))
                print("Filter: {}".format(self.filter_perc))
                for method in self.results_all.keys():
                    #print("Evaluating {} method".format(method))
                    fractions_methods[method] = self.get_performance_all("DMR based", method)
                    print("Done evaluating {} method".format(method))
                self.compare_conf_matrix("DMR based")
                self.compare_performance("DMR based")
                self.compare_fractions(fractions_methods)
                print("STEP 4: OVERLAP ALL AGAINST ALL")
                all_against_all_overlaps = []
                list_of_all_tabs = {method: self.results_all[method].processed_table for method in self.results_all.keys()}
                list_of_all_tabs["Gold-Standard"] = self.simulated.processed_table
                all_against_all_ranks = []
                for method1 in list_of_all_tabs.keys():
                    for method2 in list_of_all_tabs.keys():
                        all_against_all_overlaps.append(self.find_overlaps_4_methods(method1, method2,list_of_all_tabs[method2],list_of_all_tabs[method1]))
                        if method1 != "Gold-Standard" and method2 != "Gold-Standard":
                            all_against_all_ranks.append(self.find_overlaps_4_ranks(method1, method2,list_of_all_tabs[method2],list_of_all_tabs[method1]))
                all_against_all_overlaps_df = pd.DataFrame(all_against_all_overlaps)
                all_against_all_overlaps_df.to_csv(os.path.join(self.complete_analysis_folder,  "all_against_all_overlaps_df.csv"))
                self.plot_heatmap(all_against_all_overlaps_df)
                all_against_all_ranks_df = pd.DataFrame(all_against_all_ranks)
                all_against_all_ranks_df.to_csv(os.path.join(self.complete_analysis_folder, "all_against_all_ranks_df.csv"))
                self.plot_rank_heatmaps(all_against_all_ranks_df)

        #TODO rank
    def plot_roc_curve(self,tp_df, fp_df, method, mode):
        # Concatenate TPs and FPs into a single dataframe
        tp_df['label'] = [1] * len(tp_df)
        fp_df['label'] = [0] * len(fp_df)
        combined_df = pd.concat([tp_df, fp_df], ignore_index=True)

        # Extract the p-value column as scores
        y_true = combined_df['label']

        score_type = self.results_all[method].score_type

        if score_type != 'meanDiff':
            y_scores = combined_df['score']
        else:
            y_scores = combined_df['score'].abs()

        # Calculate ROC curve
        fpr, tpr, thresholds = roc_curve(y_true, y_scores)
        roc_auc = auc(fpr, tpr)

        # Plot ROC curve
        plt.figure(figsize=(8, 6))
        plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
        plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')  # Diagonal line (no skill classifier)
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver Operating Characteristic (ROC) Curve (method: '+method+"; score_type: "+ score_type +")")
        plt.legend(loc="lower right")

        if mode == "DMR based":
            savepath = os.path.join(self.complete_analysis_folder,  method +"_dmr_based"+'_ROC.png')
        elif mode == "CpG based":
            savepath = os.path.join(self.complete_analysis_folder, method +"_cpg_based"+'_ROC.png')
        #savepath = os.path.join(self.complete_analysis_folder, method +'_ROC.png')
        plt.savefig(savepath, dpi=300)
        #plt.show()

    def plot_precision_recall_curve(self, tp_df, fp_df, fn_df, method, mode):
        # Concatenate TPs, FPs, and FNs into a single dataframe
        tp_df['label'] = 1  # True Positives are labeled as 1
        fp_df['label'] = 0  # False Positives are labeled as 0

        # Assign an approximate p-value for the False Negatives
        score_type = self.results_all[method].score_type
        score_threshold = self.results_all[method].score_threshold
        if score_type != 'meanDiff':
            fn_df['score'] = np.random.uniform(score_threshold, 1.0, len(fn_df))
        else:
            fn_df['score'] = np.random.uniform(0, score_threshold, len(fn_df))
        fn_df['label'] = 1  # False Negatives are labeled as 1 (actual positives)

        # Combine all dataframes
        combined_df = pd.concat([tp_df, fp_df, fn_df], ignore_index=True)

        # Extract labels and p-value scores
        y_true = combined_df['label']
        if score_type != 'meanDiff':
            y_scores = combined_df['score']
        else:
            y_scores = combined_df['score'].abs()

          # Use the actual p-values for TPs, FPs, and assigned for FNs

        # Calculate Precision-Recall curve
        precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
        pr_auc = auc(recall, precision)

        # Plot Precision-Recall curve
        plt.figure(figsize=(8, 6))
        plt.plot(recall, precision, color='darkorange', lw=2, label=f'Precision-Recall curve (AUC = {pr_auc:.2f})')
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title('Precision-Recall Curve (method: '+method+"; score_type: "+ score_type +")")
        plt.legend(loc="lower left")

        if mode == "DMR based":
            savepath = os.path.join(self.complete_analysis_folder,  method +"_dmr_based"+'_PRC.png')
        elif mode == "CpG based":
            savepath = os.path.join(self.complete_analysis_folder, method +"_cpg_based"+'_PRC.png')
        plt.savefig(savepath, dpi=300)

    def getAUC(self): #tonight
        pass

    def getVenn(self):
        #overlaps 1,5,10,50,100
        pass

    def rankComparison(self):
        #TODO rank comparison
        pass

    def compare_fractions(self, fracs): #"all_overlaps" : NEW_TP, "overlap_per_sim" : NEW_per_sim
        fracts_TP = {method: fracs[method]["all_overlaps"]["fraction_overlap"]
                    for method in fracs.keys()}
        fracts_TP_per_sim = {method: fracs[method]["overlap_per_sim"]["fraction_overlap"]
                     for method in fracs.keys()}
        fracts_TP_df = pd.DataFrame(fracts_TP)
        #print("Fract Overlap Plot")

        #print(fracts_TP_df)
        fracts_TP_per_sim_df = pd.DataFrame(fracts_TP_per_sim)
        fracts_TP_melted = pd.melt(fracts_TP_df, var_name='Method', value_name='Fraction_Overlap')
        fracts_TP_melted.dropna(subset=['Fraction_Overlap'], inplace=True)
        fracts_TP_per_sim_melted = pd.melt(fracts_TP_per_sim_df, var_name='Method', value_name='Fraction_Overlap')
        fracts_TP_per_sim_melted.dropna(subset=['Fraction_Overlap'], inplace=True)

        DMRs_per_sim = {method: fracs[method]["overlap_per_sim"]["ID"]
                     for method in fracs.keys()}
        DMRs_per_sim_df = pd.DataFrame(DMRs_per_sim)
        DMRs_per_sim_melted = pd.melt(DMRs_per_sim_df, var_name='Method', value_name='Num_of_DMRs')
        DMRs_per_sim_melted.dropna(subset=['Num_of_DMRs'], inplace=True)

        comparison_TP = {method: [len(fracs[method]["all_overlaps"]),len(fracs[method]["overlap_per_sim"])]
                     for method in fracs.keys()}
        comparison_TP_df = pd.DataFrame(comparison_TP)
        comparison_TP_df_t = comparison_TP_df.T.reset_index()
        comparison_TP_df_t.columns = ['Method', 'Overlapping_DMR_Count','Overlap_per_True_Count']
        comparison_TP_df_melted =  pd.melt(comparison_TP_df_t, id_vars='Method', var_name='Overlap_Count', value_name='Count') #pd.melt(comparison_TP_df_t, var_name='Method', value_name='Overlap_Count')

        # Plotting a Box Plot

        plt.figure(figsize=(12, 6))
        sns.boxplot(x='Method', y='Fraction_Overlap', data=fracts_TP_melted, palette='pastel')
        plt.title('Distribution of overlap (as fraction) by Method')
        plt.xlabel('Method')
        plt.ylabel('Overlap (fraction)')
        savepath = os.path.join(self.complete_analysis_folder,  'Compare_fraction_overlap_box.png')
        plt.savefig(savepath, dpi=300, bbox_inches='tight')
        fracts_TP_df.to_csv(os.path.join(self.complete_analysis_folder,  'Compare_fraction_overlap.csv'))

        # Plotting a Box Plot
        plt.figure(figsize=(12, 6))
        sns.boxplot(x='Method', y='Fraction_Overlap', data=fracts_TP_per_sim_melted, palette='muted')
        plt.title('Distribution of overlap (as fraction) by Method (gold-standard-based)')
        plt.xlabel('Method')
        plt.ylabel('Overlap (fraction)')
        savepath2 = os.path.join(self.complete_analysis_folder,  'Compare_fraction_overlap_single_hit.png')
        plt.savefig(savepath2, dpi=300, bbox_inches='tight')
        fracts_TP_per_sim_df.to_csv(os.path.join(self.complete_analysis_folder, 'Compare_fraction_overlap_per_sim.csv'))
        comparison_TP_df.to_csv(os.path.join(self.complete_analysis_folder, 'Compare_comparison_overlaps.csv'))

        # Create the histogram fraction overlaps
        # Initialize 2x2 charts
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))
        # Flatten the axes array (makes it easier to iterate over)
        axes = axes.flatten()
        # Loop through each column and plot a histogram
        for i, column in enumerate(fracts_TP_df.columns):
            # Add the histogram
            fracts_TP_df[column].hist(ax=axes[i],  # Define on which ax we're working on
                            edgecolor='black',  # Color of the border
                            color='dodgerblue'  # Color of the bins
                            )

            # Add title and axis label
            axes[i].set_title(f'{column}: Overlap distribution')
            axes[i].set_xlabel('Overlap fraction')
            axes[i].set_ylabel('Frequency')

        plt.tight_layout()
        # Save the plot
        savepath3 = os.path.join(self.complete_analysis_folder, 'Compare_fraction_overlap_hist.png')
        plt.savefig(savepath3, dpi=300, bbox_inches='tight')

        #plot histograms single hit
        # Create the histogram fraction overlaps
        # Initialize 2x2 charts
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))
        # Flatten the axes array (makes it easier to iterate over)
        axes = axes.flatten()
        # Loop through each column and plot a histogram
        for i, column in enumerate(fracts_TP_per_sim_df.columns):
            # Add the histogram
            fracts_TP_per_sim_df[column].hist(ax=axes[i],  # Define on which ax we're working on
                            edgecolor='black',  # Color of the border
                            color='dodgerblue'  # Color of the bins
                            )

            # Add title and axis label
            axes[i].set_title(f'{column}: Overlap distribution')
            axes[i].set_xlabel('Overlap fraction')
            axes[i].set_ylabel('Frequency')

        plt.tight_layout()
        # Save the plot
        savepath3_2 = os.path.join(self.complete_analysis_folder, 'Compare_fraction_overlap_hist_single_hit.png')
        plt.savefig(savepath3_2, dpi=300, bbox_inches='tight')

        # Create the histograms per sim
        # Initialize 2x2 charts
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))
        # Flatten the axes array (makes it easier to iterate over)
        axes = axes.flatten()
        # Loop through each column and plot a histogram
        for i, column in enumerate(DMRs_per_sim_df.columns):
            # Add the histogram
            DMRs_per_sim_df[column].hist(ax=axes[i],  # Define on which ax we're working on
                            edgecolor='black',  # Color of the border
                            color='mediumseagreen'  # Color of the bins
                            )

            # Add title and axis label
            axes[i].set_title(f'{column}: Detected DMRs per true DMR')
            axes[i].set_xlabel('number of detected DMRs per true DMR')
            axes[i].set_ylabel('Frequency')

        plt.tight_layout()
        # Save the plot
        savepath4 = os.path.join(self.complete_analysis_folder, 'Compare_NrDMR_per_true_hist.png')
        plt.savefig(savepath4, dpi=300, bbox_inches='tight')

        # Plotting a Box Plot
        plt.figure(figsize=(12, 6))
        sns.boxplot(x='Method', y='Num_of_DMRs', data=DMRs_per_sim_melted, palette='muted')
        plt.title('Distribution of Detected DMRs per true DMR by Method')
        plt.xlabel('Method')
        plt.ylabel('Detected DMRs per true DMR')
        savepath5 = os.path.join(self.complete_analysis_folder, 'Compare_NrDMR_per_true_box.png')
        plt.savefig(savepath5, dpi=300, bbox_inches='tight')

        # Plotting the horizontal bar plot
        plt.figure(figsize=(12, 8))
        sns.barplot(y='Method', x='Count', hue='Overlap_Count', data=comparison_TP_df_melted, palette='muted')

        # Adding labels and title
        plt.ylabel('Methods')
        plt.xlabel('Count')
        plt.title("Comparison of overlapping DMR count vs overlapping DMR count (single hit)")
        plt.legend(title='Metric', loc='lower right')
        savepath6 = os.path.join(self.complete_analysis_folder, 'Compare_NrDMR_with_true_bar.png')
        plt.savefig(savepath6, dpi=300, bbox_inches='tight')  # Save as a PNG file with 300 DPI

        # Create subplots for a 2x2 grid of pie charts
        if not self.mouse:
            fig, axes = plt.subplots(2, 2, figsize=(12, 12))
            # Adjust the layout
            plt.subplots_adjust(hspace=0.4, wspace=0.4)
            # Plot pie charts for each method
            CompleteAnalysis.plot_pie_chart(axes[0, 0], fracs['DiMmer']["overlap_per_sim"], 'DiMmer', fracs['DiMmer']["all_sim"])
            CompleteAnalysis.plot_pie_chart(axes[0, 1], fracs['dmrseq']["overlap_per_sim"], 'dmrseq', fracs['dmrseq']["all_sim"])
            CompleteAnalysis.plot_pie_chart(axes[1, 0], fracs['DMRcate']["overlap_per_sim"], 'DMRcate', fracs['DMRcate']["all_sim"])
            CompleteAnalysis.plot_pie_chart(axes[1, 1], fracs['BSmooth']["overlap_per_sim"], 'BSmooth', fracs['BSmooth']["all_sim"])
            # Add a main title for the composed plot
            fig.suptitle('Detected DMRs Across Methods', fontsize=16)
            savepath7 = os.path.join(self.complete_analysis_folder, 'Compare_detected_perc.png')
            plt.savefig(savepath7, dpi=300, bbox_inches='tight')
        elif self.mouse:
            # Create subplots for a 1x2 grid of pie charts
            fig, axes = plt.subplots(1, 2, figsize=(12, 6))
            # Adjust the layout
            plt.subplots_adjust(wspace=0.4)
            CompleteAnalysis.plot_pie_chart(axes[0], fracs['DiMmer']["overlap_per_sim"], 'DiMmer',
                                            fracs['DiMmer']["all_sim"])
            CompleteAnalysis.plot_pie_chart(axes[1], fracs['dmrseq']["overlap_per_sim"], 'dmrseq',
                                            fracs['dmrseq']["all_sim"])
            savepath7 = os.path.join(self.complete_analysis_folder, 'Compare_detected_perc.png')
            plt.savefig(savepath7, dpi=300, bbox_inches='tight')

    def compare_Nr_DMRs(self): # tonight
        nrs_dmrs = {method: [self.individual_analysis_results[method].Nr_DMRs]
                    for method in self.individual_analysis_results.keys()}
        #nrs_dmrs["Gold-Standard"] =
        nrs_dmrs_df = pd.DataFrame(nrs_dmrs)
        #print(nrs_dmrs_df)
        df_t = nrs_dmrs_df.T.reset_index()
        df_t.columns = ['Method', 'DMR_Count']

        # Plotting with Seaborn
        plt.figure(figsize=(10, 6))
        sns.barplot(x='DMR_Count', y='Method', data=df_t, palette='muted')

        # Adding labels and title
        plt.xlabel('Number of DMRs Found')
        plt.ylabel('Methods')
        plt.title('Comparison of DMRs Found by Different Methods')

        # Display the plot
        savepath = os.path.join(self.complete_analysis_folder,  'Compare_Nr_DMRs.png')
        plt.savefig(savepath, dpi=300, bbox_inches='tight')  # Save as a PNG file with 300 DPI
        nrs_dmrs_df.to_csv(os.path.join(self.complete_analysis_folder, 'Compare_Nr_DMRs.csv'))


    def compare_lengths(self): # tonight
        lengths = {method: self.individual_analysis_results[method].tab['width']
                    for method in self.individual_analysis_results.keys()}
        lengths_df = pd.DataFrame(lengths)
        # Melting the dataframe to convert it into long format for Seaborn

        df_melted = pd.melt(lengths_df, var_name='Method', value_name='DMR_Length')
        df_melted.dropna(subset=['DMR_Length'], inplace=True)
        #print(df_melted.head())
        df_melted['DMR_Length'] = df_melted['DMR_Length'].astype(dtype="float64")#pd.to_numeric(df_melted['DMR_Length'], errors='coerce')#
        # Plotting a Violin Plot
        plt.figure(figsize=(12, 6))
        sns.violinplot(x='Method', y='DMR_Length', data=df_melted, palette='muted')
        plt.title('Distribution of DMR Widths by Method')
        plt.xlabel('Method')
        plt.ylabel('DMR Length')
        savepath1 = os.path.join(self.complete_analysis_folder,  'Compare_Lengths_Violin.png')
        plt.savefig(savepath1, dpi=300, bbox_inches='tight')  # Save as a PNG file with 300 DPI
        #plt.show()
        # Plotting a Box Plot
        plt.figure(figsize=(12, 6))
        sns.boxplot(x='Method', y='DMR_Length', data=df_melted, palette='pastel')
        plt.title('Distribution of DMR Widths by Method')
        plt.xlabel('Method')
        plt.ylabel('DMR Length')
        savepath2 = os.path.join(self.complete_analysis_folder,  'Compare_Lengths_Box.png')
        plt.savefig(savepath2, dpi=300, bbox_inches='tight')  # Save as a PNG file with 300 DPI
        lengths_df.to_csv(os.path.join(self.complete_analysis_folder, 'Compare_lengths.csv'))

    # Function to normalize scores to [0, 1]
    def normalize_scores(self, scores, score_type): #min-max normalization, alternatives: robust scaling, logarithmic scaling, max abs scaling, z-score normalization
        if score_type == 'meanDiff':
            scores = abs(scores)
            min_val, max_val = 0.1, max(scores)
        else:
            min_val, max_val = 0, 0.05
        return (scores - min_val) / (max_val - min_val)


    def compare_scores(self):
        scores = {method: self.normalize_scores(self.individual_analysis_results[method].tab['score'],
                                                self.individual_analysis_results[method].score)
                    for method in self.individual_analysis_results.keys() if method != "Gold-Standard"}
        scores_df = pd.DataFrame(scores)
        df_melted = pd.melt(scores_df, var_name='Method', value_name='Normalized_Score')
        # Set a common style for all plots
        sns.set(style="whitegrid")
        # 1. Violin Plot
        plt.figure(figsize=(10, 6))
        sns.violinplot(x='Method', y='Normalized_Score', data=df_melted, palette='muted')
        plt.title('Distribution of Normalized Scores by Method')
        plt.xlabel('Method')
        plt.ylabel('Normalized Score')
        savepath1 = os.path.join(self.complete_analysis_folder, 'Compare_Scores_Violin.png')
        plt.savefig(savepath1, dpi=300, bbox_inches='tight')  # Save as a PNG file with 300 DPI


        # 2. Box Plot
        plt.figure(figsize=(10, 6))
        sns.boxplot(x='Method', y='Normalized_Score', data=df_melted, palette='pastel')
        plt.title('Distribution of Normalized Scores by Method')
        plt.xlabel('Method')
        plt.ylabel('Normalized Score')
        savepath2 = os.path.join(self.complete_analysis_folder, 'Compare_Scores_Box.png')
        plt.savefig(savepath2, dpi=300, bbox_inches='tight')  # Save as a PNG file with 300 DPI

        # 3. Line Plot to show distribution
        plt.figure(figsize=(10, 6))
        sns.lineplot(data=df_melted, dashes=False, palette='dark')
        plt.title('Distribution of Normalized Scores by Method (Line Plot)')
        plt.xlabel('Index')
        plt.ylabel('Normalized Score')
        plt.legend(title='Method', loc='upper right')
        savepath3 = os.path.join(self.complete_analysis_folder, 'Compare_Scores_Line.png')
        plt.savefig(savepath3, dpi=300, bbox_inches='tight')  # Save as a PNG file with 300 DPI


    def time_series_comparison(self):
        pass

    def compare_conf_matrix(self, mode):
        if mode == "DMR based":
            TPs = pd.DataFrame.from_dict(self.TP, orient = 'index')
            TPs = TPs.transpose()
            FPs = pd.DataFrame.from_dict(self.FP, orient = 'index')
            FPs = FPs.transpose()
            FNs = pd.DataFrame.from_dict(self.FN, orient = 'index')
            FNs = FNs.transpose()
            title_plot = 'Comparison of TP, FP, and FN Across Different Methods'
            savepath = os.path.join(self.complete_analysis_folder, 'Compare_conf_matrix_dmr_based.png')
        TPs["type"] = "TP"
        FPs["type"] = "FP"
        FNs["type"] = "FN"
        if mode == "DMR based":
            all_conf_matrix = pd.concat([TPs, FPs, FNs], ignore_index=True)

        # Melting the dataframe to convert it into long format for Seaborn
        df_melted = pd.melt(all_conf_matrix, id_vars='type', var_name='Method', value_name='Count')
        #print("compare_conf_matrix")
        #print(df_melted)
        # Plotting the grouped bar plot
        plt.figure(figsize=(12, 8))
        sns.barplot(x='Method', y='Count', hue='type', data=df_melted, palette='muted')

        # Adding labels and title
        plt.xlabel('Methods')
        plt.ylabel('Count')
        plt.title(title_plot)
        plt.legend(title='Type', loc='upper right')

        # Display the plot
        savepath = os.path.join(self.complete_analysis_folder, 'Compare_conf_matrix.png')
        plt.savefig(savepath, dpi=300, bbox_inches='tight')  # Save as a PNG file with 300 DPI
        all_conf_matrix.to_csv(os.path.join(self.complete_analysis_folder, 'Compare_conf_matrix.csv'))



    def compare_performance(self, mode):
        if mode == "DMR based":
            precision = pd.DataFrame.from_dict(self.precision, orient='index')
            precision = precision.transpose()
            recall = pd.DataFrame.from_dict(self.recall, orient ='index')
            recall = recall.transpose()
            f1 = pd.DataFrame.from_dict(self.F1, orient ='index')
            f1 = f1.transpose()
            title_plot = 'Comparison of Precision, Recall, and F1 Scores Across Methods'
            title_plot2 = 'Radar Chart: Precision, Recall, and F1 Scores by Method'
            savepath1 = os.path.join(self.complete_analysis_folder, 'Compare_performances_bar_dmr_based.png')
            savepath2 = os.path.join(self.complete_analysis_folder, 'Compare_radar_chart_dmr_based.png')
        precision["type"] = "precision"
        recall["type"] = "recall"
        f1["type"] = "f1"

        if mode == "DMR based":
            all_performance_1 = pd.concat([precision, recall, f1], ignore_index=True)

        all_performance = all_performance_1.copy()
        # Melting the dataframe to convert it into long format for Seaborn
        df_melted = pd.melt(all_performance_1, id_vars='type', var_name='Method', value_name='Score')

        # Plotting the horizontal bar plot
        plt.figure(figsize=(12, 8))
        sns.barplot(y='Method', x='Score', hue='type', data=df_melted, palette='muted')

        # Adding labels and title
        plt.ylabel('Methods')
        plt.xlabel('Score')
        plt.title(title_plot)
        plt.legend(title='Metric', loc='lower right')
        #savepath1 = os.path.join(self.complete_analysis_folder, 'Compare_performances_bar.png')
        plt.savefig(savepath1, dpi=300, bbox_inches='tight')  # Save as a PNG file with 300 DPI

        # Sample data: Metrics for each method
        labels = all_performance['type']
        num_metrics = len(labels)

        # Method scores
        DiMmer_scores = all_performance['DiMmer']
        #print(DiMmer_scores)
        dmrseq_scores = all_performance['dmrseq']
        if not self.mouse:
            dmrcate_scores = all_performance['DMRcate']
            BSmooth_scores = all_performance['BSmooth']

        # Create angles for each axis
        angles = [n / float(num_metrics) * 2 * pi for n in range(num_metrics)]
        angles += angles[:1]
        #print(angles)
        # Function to plot a single radar chart
        def add_radar(ax, scores, color, label):
            #print(scores[1:])
            #scores += scores[:1]
            scores_extended = np.concatenate([scores, [scores[0]]])
            #print(scores_extended)
            ax.fill(angles, scores_extended, color=color, alpha=0.25)
            ax.plot(angles, scores_extended, color=color, linewidth=2, label=label)

        # Initialize radar chart
        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))

        # Adding radar plots for each method
        add_radar(ax, DiMmer_scores, '#2E8B57', 'DiMmer')
        add_radar(ax, dmrseq_scores, '#4682B4', 'dmrseq')
        if not self.mouse:
            add_radar(ax, dmrcate_scores, '#D2691E', 'DMRcate')
            add_radar(ax, BSmooth_scores, '#DA70D6', 'BSmooth')

        # Add labels to the axes
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(labels)

        # Add a legend and title
        plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.1))
        plt.title(title_plot2)
        #savepath2 = os.path.join(self.complete_analysis_folder, 'Compare_radar_chart.png')
        plt.savefig(savepath2, dpi=300, bbox_inches='tight')  # Save as a PNG file with 300 DPI
        all_performance.to_csv(os.path.join(self.complete_analysis_folder, 'Compare_performance.csv'))


    def find_overlaps_optimized(self,testtab, simtab):
        # Rename columns
        #testtab.columns = ['chr', 'start', 'end']
        testtab["ID"] = ['DMR_' + str(i + 1) for i in range(len(testtab))]
        simtab["ID"] = ['EVAL_' + str(i + 1) for i in range(len(simtab))]

        # Merge testtab and simtab by chromosome (seqnames)
        merged = pd.merge(testtab, simtab, on='chr', suffixes=('', '_sim'))

        # Calculate overlaps
        merged['overlap'] = np.maximum(0, np.minimum(merged['end'], merged['end_sim']) - np.maximum(merged['start'],
                                                                                                    merged[
                                                                                                        'start_sim']))

        filtered = merged[merged['overlap'] > 0].copy()

        # Calculate width and fraction_overlap
        filtered['width'] = filtered['end'] - filtered['start']
        filtered['sim_width'] = filtered['end_sim'] - filtered['start_sim']
        filtered['fraction_overlap'] = filtered['overlap'] / filtered['sim_width']
       # print("filtered:")
       # print(filtered.columns)
       # print(filtered.head())
       # print(len(filtered))
        grouped_filtered = filtered.copy()
        grouped_filtered = grouped_filtered[["ID","ID_sim","overlap","sim_width","fraction_overlap", "width"]]#"num_CpGs_sim",'num_CpGs_sim',
        grouped_filtered["ID_name"] = grouped_filtered["ID"]
        grouped_filtered = grouped_filtered.groupby(['ID_sim', "sim_width"]).aggregate({'ID': 'count',
                             'overlap': 'mean', 'fraction_overlap': 'sum', 'ID_name': '; '.join, 'width': 'sum'})

       # print("grouped_filtered:")
       # print(grouped_filtered.columns)
       # print(grouped_filtered.head())
       # print(len(grouped_filtered))

        ##chosen TP has to be a table like filtered
        chosen_TP = None
        if self.base_ref == "gsb":
            grouped_filtered_new = grouped_filtered[grouped_filtered['fraction_overlap'] >= self.filter_perc]
            id_names = grouped_filtered_new['ID_name'].str.split('; ').apply(lambda x: [name.strip() for name in x])
            # Flatten the list of lists and extract unique names
            unique_id_names = list(set([name for sublist in id_names for name in sublist]))
            chosen_TP = filtered[filtered['ID'].isin(unique_id_names)]
        elif self.base_ref == "det":
            chosen_TP = filtered[filtered['fraction_overlap'] >= self.filter_perc]



        # Identify rows with no overlap
        no_overlap = testtab[~testtab[['chr', 'start', 'end']].apply(tuple, axis=1).isin(
            chosen_TP[['chr', 'start', 'end']].apply(tuple, axis=1))]

        return {'TP': filtered, 'FP': no_overlap, 'TP_per_sim': grouped_filtered[grouped_filtered['fraction_overlap'] >= self.filter_perc], 'all_sim': simtab, 'chosen_TP': chosen_TP}



    def find_overlaps_4_methods(self, base_ref_name, tested_name,testtab1, testtab2): #testtab1 = to test; testtab2 = benchmark
        #output percentage per combination
        testtab1["ID"] = ['DMR_' + str(i + 1) for i in range(len(testtab1))]
        testtab2["ID"] = ['EVAL_' + str(i + 1) for i in range(len(testtab2))]

        # Merge testtab and simtab by chromosome (seqnames)
        merged = pd.merge(testtab1, testtab2, on='chr', suffixes=('', '_sim'))

        # Calculate overlaps
        merged['overlap'] = np.maximum(0, np.minimum(merged['end'], merged['end_sim']) - np.maximum(merged['start'],
                                                                                                    merged[
                                                                                                        'start_sim']))
        filtered = merged[merged['overlap'] > 0].copy()
        # Calculate width and fraction_overlap
        filtered['width'] = filtered['end'] - filtered['start']
        filtered['sim_width'] = filtered['end_sim'] - filtered['start_sim']
        filtered['fraction_overlap'] = filtered['overlap'] / filtered['sim_width']
        # print("filtered:")
        # print(filtered.columns)
        # print(filtered.head())
        # print(len(filtered))
        grouped_filtered = filtered.copy()
        grouped_filtered = grouped_filtered[
            ["ID", "ID_sim", "overlap", "sim_width", "fraction_overlap", "width"]]#"num_CpGs_sim", 'num_CpGs_sim',
        grouped_filtered["ID_name"] = grouped_filtered["ID"]
        grouped_filtered = grouped_filtered.groupby(['ID_sim',  "sim_width"]).aggregate({'ID': 'count',
                                                                                                        'overlap': 'mean',
                                                                                                        'fraction_overlap': 'sum',
                                                                                                        'ID_name': '; '.join,
                                                                                                        'width': 'sum'})
        #print("Overlap:")
        #print(len(grouped_filtered))
        grouped_filtered_new = grouped_filtered[grouped_filtered['fraction_overlap'] >= self.filter_perc]
        fract = len(grouped_filtered_new)/ len(testtab2)
        if base_ref_name == tested_name:
            fract = 1

        #print(len(testtab2))
        return {"Overlap":len(grouped_filtered_new), "Fraction":fract, "Base_Ref": len(testtab2), "Tested": len(testtab1),
                "Name_Base_Ref": base_ref_name, "Name_Tested": tested_name}


    def plot_heatmap(self, overlaps_df):
        # Pivot the DataFrame for the heatmap
        pivot_df = overlaps_df.pivot("Name_Tested", "Name_Base_Ref", "Fraction")

        # Create a heatmap
        # Convert fractions to percentages
        pivot_df_percentage = pivot_df * 100

        # Create a heatmap
        plt.figure(figsize=(8, 6))
        sns.heatmap(pivot_df_percentage, annot=True, fmt=".1f", cmap='coolwarm', cbar_kws={'label': 'Percentage (%)'}, linewidths = 0.5, linecolor = 'darkgrey')

        # Add labels and title
        plt.title("Heatmap of Overlapping Detections", fontsize=16)
        plt.xlabel("Methods (Base reference)", fontsize=12)
        plt.ylabel("Methods", fontsize=12)

        # Show the plot
        savepath1 = os.path.join(self.complete_analysis_folder, 'Compare_all_against_all_heat.png')
        plt.savefig(savepath1, dpi=300, bbox_inches='tight')

    def plot_rank_heatmaps(self, ranks_df):
        cp_ranks = ranks_df.copy()
        pivot_sp_df = cp_ranks.pivot("Name_Tested", "Name_Base_Ref", "Spearman_corr")#Kendall_corr
        pivot_ken_df = cp_ranks.pivot("Name_Tested", "Name_Base_Ref", "Kendall_corr")

        # Create a heatmap
        plt.figure(figsize=(8, 6))
        sns.heatmap(pivot_sp_df, annot=True, cmap='PiYG', vmin=-1, vmax=1, cbar=True)

        # Labeling axes and plot
        plt.title("Spearman's Rank Correlation Matrix", fontsize=16)
        plt.xlabel("Methods (Base reference)", fontsize=12)
        plt.ylabel("Methods", fontsize=12)
        savepath1 = os.path.join(self.complete_analysis_folder, 'Compare_all_against_all_rank_spear_heat.png')
        plt.savefig(savepath1, dpi=300, bbox_inches='tight')

        # Create a heatmap
        plt.figure(figsize=(8, 6))
        sns.heatmap(pivot_ken_df, annot=True, cmap='RdBu', vmin=-1, vmax=1, cbar=True)

        # Labeling axes and plot
        plt.title("Kendall’s Tau Matrix", fontsize=16)
        plt.xlabel("Methods (Base reference)", fontsize=12)
        plt.ylabel("Methods", fontsize=12)
        savepath2 = os.path.join(self.complete_analysis_folder, 'Compare_all_against_all_rank_ken_heat.png')
        plt.savefig(savepath2, dpi=300, bbox_inches='tight')



    @staticmethod
    def plot_pie_chart(ax, ov_df, method_name, all_df):
        #print("plot_pie_chart")
        ov_df = ov_df.reset_index()
        #print(ov_df['ID_sim'])
        #print(all_df['ID'])
        matched = ov_df['ID_sim'].isin(all_df['ID']).sum()
        #print(matched)
        total = len(all_df)
        not_matched = total - matched

        # Prepare data for the pie chart
        labels = ['Detected DMRs', 'Not detected DMRs']
        sizes = [matched, not_matched]

        # Create a pie chart
        ax.pie(sizes, labels=labels, autopct='%1.1f%%', colors=['#66b3ff', '#ff9999'], startangle=140, explode=(0.1, 0))
        ax.set_title(f'{method_name}')

    def getFN(self,overlaps, sim): # found_overlaps,simtab_dmr
        # Select relevant columns from found_overlaps
        #print("columns")
        #print(overlaps.columns)
        sim_overlap = overlaps[['chr', 'start_sim', 'end_sim']].copy()
        sim_overlap.columns = ['chr', 'start', 'end']  # Rename columns to match simtab_dmr

        # Perform anti-join: keep rows in simtab_dmr that are not in sim_overlap
        #sim = self.simulated.processed_table
        sim_no_overlap = pd.merge(sim, sim_overlap, on=['chr', 'start', 'end'], how='outer',
                                      indicator=True)
        sim_no_overlap = sim_no_overlap[sim_no_overlap['_merge'] == 'left_only'].drop(columns=['_merge'])
        return sim_no_overlap


    def find_overlaps_4_ranks(self, base_ref_name, tested_name,testtab1, testtab2): #testtab1 = to test; testtab2 = benchmark
        #output percentage per combination
        testtab1["ID"] = ['DMR_' + str(i + 1) for i in range(len(testtab1))]
        testtab2["ID"] = ['EVAL_' + str(i + 1) for i in range(len(testtab2))]

        testtab1 = self.add_rank(testtab1,tested_name)
        testtab2 = self.add_rank(testtab2,base_ref_name)

        # Merge testtab and simtab by chromosome (seqnames)
        merged = pd.merge(testtab1, testtab2, on='chr', suffixes=('', '_sim'))

        # Calculate overlaps
        merged['overlap'] = np.maximum(0, np.minimum(merged['end'], merged['end_sim']) - np.maximum(merged['start'],
                                                                                                    merged[
                                                                                                        'start_sim']))
        filtered = merged[merged['overlap'] > 0].copy()
        #print(filtered.columns)
        # Calculate width and fraction_overlap
        filtered['width'] = filtered['end'] - filtered['start']
        filtered['sim_width'] = filtered['end_sim'] - filtered['start_sim']
        filtered['fraction_overlap'] = filtered['overlap'] / filtered['sim_width']
        # print("filtered:")
        # print(filtered.columns)
        # print(filtered.head())
        # print(len(filtered))
        grouped_filtered = filtered.copy()
        grouped_filtered = grouped_filtered[
            ["ID", "ID_sim",  "overlap", "sim_width", "fraction_overlap","rank","rank_sim", "width"]] #"num_CpGs_sim",'num_CpGs_sim',
        grouped_filtered["ID_name"] = grouped_filtered["ID"]
        grouped_filtered = grouped_filtered.groupby(['ID_sim', "sim_width","rank_sim"], as_index=False).aggregate({'ID': 'count',
                                                                                                        'overlap': 'mean',
                                                                                                        'fraction_overlap': 'sum',
                                                                                                        'ID_name': '; '.join,
                                                                                                        'rank': 'min',
                                                                                                        'width': 'sum'
                                                                                                        })
        #activate if doesn't work
        #print(grouped_filtered.columns)
        grouped_filtered_new = grouped_filtered[grouped_filtered['fraction_overlap'] >= self.filter_perc]
        #grouped_filtered_new["rank"] = grouped_filtered_new["rank"].round()
        # Spearman’s Rank Correlation
        spearman_corr, spearman_p_value = spearmanr(grouped_filtered_new['rank_sim'], grouped_filtered_new['rank'])
        #print(f"Spearman’s Rank Correlation: {spearman_corr:.3f}, p-value: {spearman_p_value:.3f}")

        # Kendall’s Tau
        kendall_corr, kendall_p_value = kendalltau(grouped_filtered_new['rank_sim'], grouped_filtered_new['rank'])
        #print(f"Kendall’s Tau: {kendall_corr:.3f}, p-value: {kendall_p_value:.3f}")
        #print("Ranking:")
       # print(grouped_filtered)
        #print(len(testtab2))
        if base_ref_name == tested_name:
            spearman_corr, kendall_corr = 1,1
        return {"Overlap":len(grouped_filtered_new), "Fraction":len(grouped_filtered_new)/ len(testtab2), "Base_Ref": len(testtab2), "Tested": len(testtab1),
                "Name_Base_Ref": base_ref_name, "Name_Tested": tested_name, "Spearman_corr": spearman_corr, "Spearman_p_value": spearman_p_value,
                "Kendall_corr": kendall_corr, "Kendall_p_value": kendall_p_value}


    def add_rank(self, tab, method):
        if method == "DiMmer":
            tab['rank'] = tab['score'].rank(ascending=True, method='min')
        if method == "dmrseq":
            tab['rank'] = tab['score'].rank(ascending=True, method='min')
        if method == "DMRcate":
            tab['rank'] = tab['score_rank'].rank(ascending=True, method='min')
        if method == "BSmooth":
            #tab = tab.sort_values(by=['score', 'region'], ascending=[True, True])  # Sort by score, then by region
            tab['rank'] = tab['score'].abs().rank(ascending=False, method='min')
        return tab




