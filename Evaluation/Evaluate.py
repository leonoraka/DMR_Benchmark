#from .ReadResults import ReadResults
#from .Individual_Analysis import IndividualAnalysis
#from .CompleteAnalysis import CompleteAnalysis
#import snakemake
import sys
import os

#import snakemake

# Add the parent directory to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
#import snakemake

from Evaluation.CompleteAnalysis import CompleteAnalysis

#file_inputs = {"DiMmer": "/Users/leonoraka/Desktop/MA-project/Snakemake/Dimmer_results/240724015239739_DMRSearch/merged_table.csv" ,
#                "dmrseq": "/Users/leonoraka/Desktop/MA-project/Snakemake/Results/DMRseq-output-significant.tsv",
#                "DMRcate": "/Users/leonoraka/Desktop/MA-project/Snakemake/Results/DMRcate-output.tsv",
#                "BSmooth": "/Users/leonoraka/Desktop/MA-project/Snakemake/Results/BSseq-output-significant.tsv",
#                "Simulated": "/Users/leonoraka/Desktop/MA-project/Dataset/Simulated_DMR_dt.tsv",
#                "input_simulated": "/Users/leonoraka/Desktop/MA-project/Dataset/Simulated_filtered_all_CpGs.tsv"
#}

file_inputs = {"DiMmer": snakemake.input["Dimmer"],
                "dmrseq": snakemake.input["DMRseq"],
                "DMRcate": snakemake.input["DMRcate"],
                "BSmooth": snakemake.input["BSmooth"],
                "Simulated": snakemake.input["Simulated"],
                "input_simulated": snakemake.input["input_simulated"]
}
print("This is the input file Dimmer: ", snakemake.input["Dimmer"])
print("This is the input file DMRseq: ", snakemake.input["DMRseq"])

outfolder = snakemake.output["outfolder"] #"/Users/leonoraka/Desktop/MA-project/Evaluation_Results"


complete_analysis = CompleteAnalysis(file_inputs, outfolder,"sim")

complete_analysis.analysis()



































