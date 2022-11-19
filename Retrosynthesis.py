import os

def Preprocessing(DatasetFolder):
    import pandas as pd
    from joblib import Parallel, delayed
    from ReactionExtractor import Extractor
    Datasets  = [pd.read_csv(datafile) for datafile in os.listdir(DatasetFolder)]
    Datasets  = pd.concat(Datasets)
    Datasets  = Datasets.drop_duplicates(subset=['ReactionSmiles'])
    DatasetSplit = Datasets['ReactionSmiles'].str.split('>', expand=True)
    Datasets['Reactants']  = DatasetSplit[0]
    Datasets['Spectators'] = DatasetSplit[1]
    Datasets['Products']   = DatasetSplit[2]
    Datasets  = Datasets.reset_index().rename(columns={'Index':'ID'})
    Datasets  = Datasets[['ID', 'Reactants', 'Spectators', 'Products']]
    Templates = Datasets.to_json(orient='records')
    Templates = Parallel(n_jobs=-1)(delayed(Extractor)(reaction) for reaction in Templates)
    Templates = pd.read_json(Templates, orient='records')
    Templates = Templates.dropna(subset=['Reaction_SMARTS'])
    Rule      = Templates[['Reaction_SMARTS']].groupby('Reaction_SMARTS', as_index=False).size()
    Rule      = Templates.reset_index().rename(columns={'index':'ID'})
    Rule      = Templates[['Reaction_SMARTS', 'size']]
    Rule['Rule_Index'] = [i for i in range(len(Rule))]
    Templates = [['Reactants', 'Products', 'Reaction_SMARTS']]
    Templates['Reaction_Label'] = Templates['Reaction_SMARTS'].astype('category')
    Templates['Reaction_Label'] = Templates['Reaction_Label'].cat.set_categories(Rule['Reaction_SMARTS'])
    Templates = Templates.dropna(subset=['Reaction_Label'])
    Templates['Reaction_Label'].replace({rule:idx for idx, rule in enumerate(Rule['Reaction_SMARTS'])}, inplace=True)
    Templates = Templates.to_json(orient='records')
    Rule      = Rule.to_json(orient='records')
    return Templates, Rule