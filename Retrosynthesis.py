import os

def Preprocessing(DatasetFolder):
    import pandas as pd
    Datasets = [pd.read_csv(datafile) for datafile in os.listdir(DatasetFolder)]
    Datasets = pd.concat(Datasets)
    Datasets = Datasets.drop_duplicates(subset=['ReactionSmiles'])
    DatasetSplit = Datasets['ReactionSmiles'].str.split('>', expand=True)
    Datasets['Reactants']  = DatasetSplit[0]
    Datasets['Spectators'] = DatasetSplit[1]
    Datasets['Products']   = DatasetSplit[2]
    Datasets = Datasets.reset_index().rename(columns={'Index':'ID'})
    Datasets = Datasets[['ID', 'Reactants', 'Spectators', 'Products']]
    return Datasets