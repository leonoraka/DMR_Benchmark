import pandas as pd

#TODO: nachschauen wie wird width bei allen berechnet
#TODO: Strand abchecken
#TODO: Filter significance
class ReadResults:
    needed_colnames = ["chr", "start","end","num_CpGs","width","score"]

    def __init__(self, in_path, type):
        self.in_path = in_path
        self.method_type = type
        self.processed_table = None
        self.score_type = None
        self.score_threshold = None
        if type == "DiMmer":
            self.read_DiMmer()
        elif type == "dmrseq":
            self.read_dmrseq()
        elif type == "DMRcate":
            self.read_DMRcate()
        elif type == "BSmooth":
            self.read_BSmooth()
        elif type == "Simulated":
            self.read_Simulated()
        elif type == "input_simulated":
            self.read_input_simulated()
        elif type == "Mouse_True":
            self.read_Mouse("DMRs")
        elif type == "Mouse_Known":
            self.read_Mouse("DMRs_Known")

    def read_DiMmer(self): #gefiltert?
        """a DiMmer results table looks like:
        -> Chr, Begin, End, begin.CpG, end.CpG, score, #CpG, Num.DMRs, Average.DMRs, p-value, log.ratio, link"""
        table = pd.read_csv(self.in_path,skipinitialspace = True)
        table["width"] = table["End"].astype("Int64") - table["Begin"].astype("Int64")
        new_tab = table.copy()
        new_tab = new_tab[["Chr","Begin","End","#CpG","width","p-value"]]
        new_tab.columns = ReadResults.needed_colnames
        self.processed_table = new_tab
        self.score_type = "p-value"
        self.score_threshold = 0.05
        #new_tab.rename(columns ={ new_tab.columns[new_name_index] : ReadResults.needed_colnames[new_name_index] for new_name_index in range(len(ReadResults.needed_colnames))})

    def read_dmrseq(self):
        """a dmrseq results table looks like:
        -> seqnames, start, end, width, strand, L, area, ß, stat, pval, qval, index.start, index.end, index.width"""
        table = pd.read_table(self.in_path, sep=",")
        new_tab = table.copy()
        new_tab = new_tab[["seqnames", "start", "end", "L", "width", "qval"]]
        new_tab.columns = ReadResults.needed_colnames
        self.processed_table = new_tab
        self.score_type = "fdr" #actually fdr
        self.score_threshold = 0.05



    def read_DMRcate(self):
        """a DMRcate results table looks like:
        -> chr, start, end, width, no.cpgs, min_smoothed_fdr, Stouffer, HMFDR, Fisher, maxdiff, meandiff, overlapping genes
        min_smoothed_fdr = Minimum FDR of the smoothed estimate
        Stouffer = Stouffer summary transform of the individual CpG FDRs.
        HMFDR = Harmonic mean of the individual CpG FDRs. => hereby in script you chose fdr cutoff <0.05
        Fisher = Fisher combined probability transform of the individual CpG FDRs.
        => schon gefiltert"""
        table = pd.read_table(self.in_path, sep=";")
        new_tab = table.copy()
        new_tab = new_tab[["chr","start","end","no.cpgs","width","min_smoothed_fdr","Fisher"]] #Fisher shows rank; min_smoothed_fdr = significant
        new_tab.columns = ReadResults.needed_colnames + ["score_rank"]
        self.processed_table = new_tab
        self.score_type = "fdr"#"meanDiff"#Fisher" #HMFDR
        self.score_threshold = 0.05 #idk?


    def read_BSmooth(self):
        """a BSmooth results table looks like:
        -> chr, start, end, idxStart, idxEnd, cluster, n, width, invdensity, areaStat, maxStat, meanDiff, group1.mean, group2.mean, tstat.sd, direction"""
        table = pd.read_table(self.in_path, sep=" ")
        new_tab = table.copy()
        new_tab = new_tab[["chr", "start", "end", "n", "width", "meanDiff"]]  # or stouffer?
        new_tab.columns = ReadResults.needed_colnames
        self.processed_table = new_tab
        self.score_type = "meanDiff" #n >= 3 & abs(meanDiff) >= 0.1
        self.score_threshold = 0.1 # aber GROESSER GLEICH


    def read_Simulated(self):
        table = pd.read_table(self.in_path, sep="\t")
        new_tab = table.copy()
        new_tab = new_tab[["seqnames", "start", "end",  "L", "width", "delta", "mncov"]]# "strand" ["chr", "start","end","num_CpGs","width","score"]
        new_tab.columns = ReadResults.needed_colnames[:5] + ["delta", "mncov"]
        self.processed_table = new_tab


    def read_input_simulated(self):
        """a input simulated results table looks like:
        -> seqnames	start	end	width	strand"""
        table = pd.read_table(self.in_path, sep="\t")
        new_tab = table.copy()
        new_tab.columns = ["chr", "start", "end", "width", "strand"]
        self.processed_table = new_tab

    def read_Mouse(self, sheet):
        table = pd.read_excel(self.in_path, sheet_name=sheet, engine='openpyxl')
        new_tab = table.copy()
        new_tab = new_tab[["Chr", "Start", "End", "Length (bp)"]]# "strand"
        new_tab.columns = ["chr", "start","end","width"]
        self.processed_table = new_tab