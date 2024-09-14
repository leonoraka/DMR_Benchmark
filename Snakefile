import glob
import os
wd_path = '/Users/leonoraka/Desktop/MA-project/Snakemake/Analysis/'
input_path = wd_path + '01_Input/'
dataset_path = wd_path + '02_Dataset/'
simulation_path = wd_path + '03_Simulation/'
result_path = wd_path + '04_Result/'
evaluation_path = wd_path + '05_Evaluation/'
evaluation_mouse_path = wd_path + '05_Evaluation_Mouse/'
benchmarking_path = wd_path + '06_Keep_Track/'

numDMRs_arr = [1500,3000,4500,6000] #6000, 4500, 6000
#unknown_wildcards = glob_wildcards(result_path + 'Res_{dmrs}/dimmer_results/{unknown}_DMRSearch/merged_table.csv')
folders = glob.glob(result_path + 'Res_{dmrs}/dimmer_results/*_DMRSearch/merged_table.csv')
#print(unknown_wildcards)
print(result_path + 'Res_{dmrs}/dimmer_results/*_DMRSearch/merged_table.csv')

def get_reads(wildcards):
    res = result_path + 'Res_{dmrs}/dimmer_results/'
    if os.path.isdir(f'/Users/leonoraka/Desktop/MA-project/Snakemake/Analysis/04_Result/Res_{wildcards.dmrs}/dimmer_results/'):
        readNames =  [f"/Users/leonoraka/Desktop/MA-project/Snakemake/Analysis/04_Result/Res_{wildcards.dmrs}/dimmer_results/{read}/merged_table.csv" for read in os.listdir(f'/Users/leonoraka/Desktop/MA-project/Snakemake/Analysis/04_Result/Res_{wildcards.dmrs}/dimmer_results/') if "_DMRSearch" in read]
        res = readNames[0]
    return res
    
#current_res = result_path + 'Res_{dmrs}/'
#print(current_res)
#rule all:
#    input:
#        expand(simulation_path + 'Sim_{dmrs}/Simulated_BSseq_obj.rds', dmrs=numDMRs_arr),
#        expand(result_path + 'Res_{dmrs}/dimmer_results/{unknown_part}_DMRSearch/merged_table.csv',
#               dmrs=unknown_wildcards.dmrs, unknown_part=unknown_wildcards.unknown_part),
#        expand(evaluation_path + 'Eval_{dmrs}', dmrs=numDMRs_arr)
rule all:
    input:
        #print(unknown_wildcards),
        expand(simulation_path + 'Sim_{dmrs}/Simulated_BSseq_obj.rds', dmrs=numDMRs_arr),
        expand(result_path + 'Res_{dmrs}/dimmer_results/dimmer_project.csv', dmrs=numDMRs_arr),
        expand(result_path + 'Res_{dmrs}/DMRseq-output.tsv', dmrs=numDMRs_arr),
        expand(result_path + 'Res_{dmrs}/DMRcate-output.tsv', dmrs=numDMRs_arr),
        expand(result_path + 'Res_{dmrs}/BSseq-output.tsv', dmrs=numDMRs_arr),
        #expand(result_path + 'Res_{dmrs}/dimmer_results/{unknown}_DMRSearch/merged_table.csv', dmrs=numDMRs_arr,unkown =unknown_wildcards.unknown),
        expand(evaluation_path + 'Eval_{dmrs}/CompleteAnalysis/', dmrs=numDMRs_arr),
        evaluation_mouse_path + "All",
        evaluation_mouse_path + "Known"


rule Get_Dataset:
    input: 
        prostate_cancer = input_path + 'BigTable.tsv',
        prostate_cancer_meta = input_path + 'GSE158927_series_matrix.xlsx'
    output:
        whole = dataset_path + 'Whole_BSseq_obj.rds',
        filtered = dataset_path + 'Filtered_BSseq_obj.rds',
        smoothed = dataset_path + 'Smoothed_BSseq_obj.rds',
        dimmer_folder = dataset_path + 'Dimmer_Dataset/'
    script:
        '/Users/leonoraka/Desktop/MA-project/Snakemake/Get_Dataset.R'   


rule Simulation:
    input: 
        rds = dataset_path + 'Filtered_BSseq_obj.rds'
    output:
        out_folder = simulation_path + 'Sim_{dmrs}/Simulated_dimmer/',
        out_obj = simulation_path + 'Sim_{dmrs}/Simulated_BSseq_obj.rds',
        filtered_out_obj = simulation_path + 'Sim_{dmrs}/Simulated_filtered_BSseq_obj.rds',
        dmr_obj = simulation_path + 'Sim_{dmrs}/Simulated_DMR_obj.rds',
        out_annot = simulation_path + 'Sim_{dmrs}/Simulated_dimmer/sample_annotation_simulated.csv',
        dmr_tab = simulation_path + 'Sim_{dmrs}/Simulated_DMR_dt.tsv',
        all_CpGs = simulation_path + 'Sim_{dmrs}/Simulated_filtered_all_CpGs.tsv'
    params:
        Covariate = 'Patient outcome',
        numDMRs = lambda wildcards: wildcards.dmrs,
        only_normal = "TRUE"
    script:
        '/Users/leonoraka/Desktop/MA-project/Snakemake/Simulation.R'

rule DMRSeq:
    input: 
        rds = simulation_path + 'Sim_{dmrs}/Simulated_filtered_BSseq_obj.rds'
    output:
        out = result_path + 'Res_{dmrs}/DMRseq-output.tsv',
        sig = result_path + 'Res_{dmrs}/DMRseq-output-significant.tsv'
    benchmark:
        benchmarking_path + "Bench_{dmrs}/DMRSeq_benchmark.txt"
    params:
        Covariate = 'condition'
    script:
        '/Users/leonoraka/Desktop/MA-project/Snakemake/Pipeline-DMR-seq.R'

rule Dimmer:
    input: 
        annot = simulation_path + 'Sim_{dmrs}/Simulated_dimmer/sample_annotation_simulated.csv',
        config = input_path + 'dimmer.config'
    output:
        out = result_path + 'Res_{dmrs}/dimmer_results/',
        out_proj = result_path + 'Res_{dmrs}/dimmer_results/dimmer_project.csv'
    benchmark:
        benchmarking_path + "Bench_{dmrs}/Dimmer_benchmark.txt"
    params:
        threads = "6",
        variable = 'condition'
    shell:
        'java -Xmx8192M -jar /Users/leonoraka/Downloads/dimmer.jar {input.config} --output_path {output.out} --variable {params.variable} --annotation_path {input.annot} --threads {params.threads} --input_type bisulfite --model T-test --save_dmr_permu_plots 1'

rule DMRcate:
    input: 
        rds = simulation_path + 'Sim_{dmrs}/Simulated_filtered_BSseq_obj.rds'
    output:
        out = result_path + 'Res_{dmrs}/DMRcate-output.tsv'
    benchmark:
        benchmarking_path + "Bench_{dmrs}/DMRcate_benchmark.txt"
    params:
        Covariate = 'condition' #patient outcome
    script:
        '/Users/leonoraka/Desktop/MA-project/Snakemake/New_DMRcate.R'

rule BSSeq:
    input: 
        rds = simulation_path + 'Sim_{dmrs}/Simulated_filtered_BSseq_obj.rds'
    output:
        out = result_path + 'Res_{dmrs}/BSseq-output.tsv',
        sig = result_path + 'Res_{dmrs}/BSseq-output-significant.tsv'
    benchmark:
        benchmarking_path + "Bench_{dmrs}/BSSeq_benchmark.txt"
    params:
        Covariate = 'condition' #patient outcome
    script:
        '/Users/leonoraka/Desktop/MA-project/Snakemake/BSseq.R' 


rule Evaluation:
    input: 
        Dimmer = get_reads,#glob.glob(dynamic(result_path + 'Res_{dmrs}/dimmer_results/*_DMRSearch/merged_table.csv')),#glob.glob(result_path + 'Res_{dmrs}/dimmer_results/*_DMRSearch/merged_table.csv'),#,
        DMRseq =  result_path + 'Res_{dmrs}/DMRseq-output-significant.tsv',
        DMRcate =  result_path + 'Res_{dmrs}/DMRcate-output.tsv',
        BSmooth =  result_path + 'Res_{dmrs}/BSseq-output-significant.tsv',
        Simulated =  simulation_path + 'Sim_{dmrs}/Simulated_DMR_dt.tsv',
        input_simulated =  simulation_path + 'Sim_{dmrs}/Simulated_filtered_all_CpGs.tsv',
        for_graph = result_path + 'Res_{dmrs}/dimmer_results/dimmer_project.csv'
    output:
        outfolder = evaluation_path + 'Eval_{dmrs}',
        check_out1 = evaluation_path + 'Eval_{dmrs}/CompleteAnalysis/'
    log:
        wd_path + 'evaluation_log.txt'
    script:
        '/Users/leonoraka/Desktop/MA-project/Snakemake/Evaluation/Evaluate.py'
    #shell:
    #    'cd /Users/leonoraka/Desktop/MA-project && python -m Evaluation.Evaluate'

#rule Evaluation_Mouse:
#    input: 
#        Mouse = 'Mouse_Results/Mouse_eval_new.xlsx',
#        DMRseq =  input_path + 'Mouse_Results/DMRseq-output.tsv',
#        Dimmer =  input_path + "Mouse_Results/merged_table.csv"
#    output:
#        outfolder = evaluation_mouse_path + "All",
#        outfolder_known = evaluation_mouse_path + "Known"
#    shell:
#        'cd /Users/leonoraka/Desktop/MA-project && python -m Evaluation.Evaluate_Mouse'

#rule end:
#    input:
#        evaluation_path + 'Eval_{dmrs}',
#        outfolder = evaluation_mouse_path + "All",
#        outfolder_known = evaluation_mouse_path + "Known"

rule summarize_benchmarks:
    input:
        expand(benchmarking_path + "Bench_{dmrs}/{rule}_benchmark.txt", rule=["DMRSeq", "Dimmer","DMRcate","BSSeq"], dmrs=numDMRs_arr)  # Add all rule names
    output:
        "benchmark_summary.txt"
    run:
        total_time = 0
        with open(output[0], 'w') as out:
            for benchmark_file in input:
                with open(benchmark_file) as f:
                    lines = f.readlines()
                    runtime_line = lines[1]  # This line has the runtime in seconds
                    runtime = float(runtime_line.split()[1])  # Extract the runtime value
                    out.write(f"{benchmark_file}: {runtime} seconds\n")
                    total_time += runtime
            out.write(f"Total runtime: {total_time} seconds\n")

