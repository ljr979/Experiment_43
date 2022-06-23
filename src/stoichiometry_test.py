import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger


output_folder = 'python_results/stoichiometry_test/'

if not os.path.exists(output_folder):
        os.makedirs(output_folder)

input_folder='python_results\stoichiometries\second_repeat\complexes/'
stoich_files =[[f'{root}/{filename}' for filename in files if 'molecule_counts.csv' in filename] for root, dirs, files in os.walk(f'{input_folder}')]
stoich_files=[item for sublist in stoich_files for item in sublist]


molecule_counts=[]

for filepath in stoich_files:
    data= pd.read_csv(filepath)
    treatments = data.treatment
    timepoints=[]
    for item in treatments:
        timepoint=item.split('-')[0]
        timepoints.append(timepoint)
    data['timepoint']=timepoints

    molecule_counts.append(data)

molecule_counts=pd.concat(molecule_counts)
molecule_counts.drop([col for col in molecule_counts.columns.tolist() if ' ' in col], axis=1, inplace=True)
molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'all_small_mol_count' in col], axis=1, inplace=True)
molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'single_step_mol_count' in col], axis=1, inplace=True)
molecule_counts.drop([col for col in molecule_counts.columns.tolist() if 'max_fluorescence' in col], axis=1, inplace=True)


#proteins = [folder for folder in os.listdir(f'{input_folder}')]

#now need to find the number of CLIC molecules that are present at each time point
timepoints = molecule_counts['timepoint'].unique().tolist()



coloc_molecules=molecule_counts[molecule_counts['colocalisation']=="Coloc"]
non_coloc_molecules=molecule_counts[molecule_counts['colocalisation']=="Non-coloc"]



shsp=coloc_molecules[coloc_molecules["protein"] == "aB-c"]
client=coloc_molecules[coloc_molecules['protein']=='CLIC']
#find
zero=[]
twenty=[]
forty=[]
sixty=[]
four_hour=[]
seven_hour=[]
for timepoint in timepoints: 
    timepoint
    protein = 'shsp'
    timepoint_shsp= shsp[shsp["timepoint"] == timepoint]
    total_shsp_molecules=sum(timepoint_shsp['last_step_mol_count'].values)
    #total_shsp_timepoint=len(timepoint_shsp)
    #coloc_timepoint_shsp=len(timepoint_shsp[timepoint_shsp["colocalisation"] == "coloc"])
    #percent_colocal_shsp_timepoint=coloc_timepoint_shsp/total_shsp_timepoint*100
    listo=[timepoint, protein, total_shsp_molecules]
    if 'zero' in timepoint:
        zero.append(listo)
    elif '20min' in timepoint:
        twenty.append(listo)
    elif '40min' in timepoint:
        forty.append(listo)
    elif '60min' in timepoint:
        sixty.append(listo)
    elif '4h' in timepoint:
        four_hour.append(listo)
    elif '7h' in timepoint:
        seven_hour.append(listo)


#find % colocalisation for CLIC over time
for timepoint in timepoints: 
    timepoint
    protein = 'client'
    timepoint_client= client[client["timepoint"] == timepoint]
    total_client_molecules = sum(timepoint_client['last_step_mol_count'].values)
    listo=[timepoint, protein, total_client_molecules]
    if 'zero' in timepoint:
        zero.append(listo)
    elif '20min' in timepoint:
        twenty.append(listo)
    elif '40min' in timepoint:
        forty.append(listo)
    elif '60min' in timepoint:
        sixty.append(listo)
    elif '4h' in timepoint:
        four_hour.append(listo)
    elif '7h' in timepoint:
        seven_hour.append(listo)

all_molecule_numbers=[zero,twenty,forty,sixty,four_hour,seven_hour]
column_names=['timepoint', 'protein','total number of proteins']
all_molecule_numbers=pd.DataFrame([item for sublist in all_molecule_numbers for item in sublist])
all_molecule_numbers.columns=column_names

for timepoint in timepoints: 
    test=all_molecule_numbers[all_molecule_numbers['timepoint']==timepoint]
    test_ratio=test['total number of proteins'].values[0]/test['total number of proteins'].values[1]
    
    listo=[timepoint, test_ratio]
    if 'zero' in timepoint:
        zero.append(listo)
    elif '20min' in timepoint:
        twenty.append(listo)
    elif '40min' in timepoint:
        forty.append(listo)
    elif '60min' in timepoint:
        sixty.append(listo)
    elif '4h' in timepoint:
        four_hour.append(listo)
    elif '7h' in timepoint:
        seven_hour.append(listo)


all_molecule_numbers=[zero,twenty,forty,sixty,four_hour,seven_hour]
column_names=['timepoint', 'protein','total number of proteins', 'stoichiometry sHsp:client']
all_molecule_numbers=pd.DataFrame([item for sublist in all_molecule_numbers for item in sublist])
all_molecule_numbers.columns=column_names




ax=plt.scatter(x=test_client['timepoint'], y=test_client['value'], c='darkviolet')
ax=plt.plot(test_client['timepoint'], test_client['value'], c='plum')
ax=plt.gca()
ax.set_ylim([0,80])
plt.ylabel('percent colocalisation')
plt.title('Percentage of Rhodanese colocalised with hsp27')
plt.xlabel('timepoint')
ax=ax.get_figure()