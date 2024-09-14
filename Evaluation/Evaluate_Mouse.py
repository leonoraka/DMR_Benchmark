import snakemake
from CompleteAnalysis import CompleteAnalysis

file_inputs_Mouse = {"Mouse": snakemake.input.Mouse,
                "dmrseq": snakemake.input.DMRseq,
                "DiMmer": snakemake.input.Dimmer,
}

outfolder_Mouse = "/Users/leonoraka/Desktop/MA-project/Evaluation_Results_Mouse"
complete_analysis_mouse = CompleteAnalysis(file_inputs_Mouse, outfolder_Mouse,"Mouse_True")
complete_analysis_mouse.analysis()

outfolder_Mouse_Known = "/Users/leonoraka/Desktop/MA-project/Evaluation_Results_Mouse"
complete_analysis_mouse_Known = CompleteAnalysis(file_inputs_Mouse, outfolder_Mouse_Known,"Mouse_Known")
complete_analysis_mouse_Known.analysis()
