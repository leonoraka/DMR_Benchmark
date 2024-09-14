import pandas as pd
import os


class FindData:
    def __init__(self,eval_path,expl_csv, mode):
        self.path = eval_path
        self.GSB = []
        self.DET = []
        self.GEN = []
        self.expl_csv = expl_csv
        self.GSB_df = None
        self.DET_df = None
        self.GEN_df = None
        if mode == "deep":
            self.get_complete_DFs()
        if mode == "general":
            self.getDF()


    def read_data(self,csv_path):
        df = pd.read_csv(csv_path)
        return df


    def get_eval_nr(self, p):
        #last_part = os.path.basename(os.path.normpath(p))#'/folderA/folderB/folderC/folderD/'
        _, number = str.split(p, sep='_')
        return number


    def get_based_perc(self,p):
        #last_part = os.path.basename(os.path.normpath(p))#'/folderA/folderB/folderC/folderD/'
        based, perc = str.split(p, sep='ased')
        end_based = ""
        if based.startswith('Det'):
            end_based = "DET"
        elif based.startswith('Gold'):
            end_based = "GSB"
        return {"BASE":end_based, "PERC":perc}


    def setUpList(self):
        for folder in os.listdir(self.path):
            print(folder)
            DMR_NR = self.get_eval_nr(folder)
            if os.path.isdir(os.path.join(self.path, folder)):
                into_CA_p = os.path.join(self.path, folder, "CompleteAnalysis")
                for other_folder in os.listdir(into_CA_p):
                    if os.path.isdir(os.path.join(into_CA_p, other_folder)):
                        get_info = self.get_based_perc(other_folder)
                        based, perc = get_info["BASE"], get_info["PERC"]
                        target_csv = os.path.join(into_CA_p, other_folder, self.expl_csv)
                        df = self.read_data(target_csv)
                        df["BASE"] = based
                        df["PERC"] = float(perc)
                        df["NR_DMR"] = int(DMR_NR)
                        if based == "DET":
                            self.DET.append(df)
                        if based == "GSB":
                            self.GSB.append(df)

    def generalList(self):
        for folder in os.listdir(self.path):
            print(folder)
            DMR_NR = self.get_eval_nr(folder)
            if os.path.isdir(os.path.join(self.path, folder)):
                needed_csv =  "CompleteAnalysis/" + self.expl_csv
                into_CA_p = os.path.join(self.path, folder, needed_csv)
                df = self.read_data(into_CA_p)
                df["NR_DMR"] = int(DMR_NR)
                self.GEN.append(df)

    def get_complete_DFs(self):
        self.setUpList()
        self.GSB_df = pd.concat(self.GSB)
        self.DET_df = pd.concat(self.DET)

    def getDF(self):
        self.generalList()
        self.GEN_df = pd.concat(self.GEN)








































