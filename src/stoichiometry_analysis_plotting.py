import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import random 
input_folder = 'raw_data/Experiment_43_2/43_clicabc temp/20220201_Experiment43_CLIC_aBc_timepoints/'
output_folder = 'python_results/experiment_43_2/stoichiometry_test/'

python_input = 'python_results/stoichiometries/complexes/'
if not os.path.exists(output_folder):
        os.makedirs(output_folder)

#first, need to go into the trajectories folders, and grab the COLOCALISED trajectories. maybe need to do this file by file before doing all of them
#assign a variable to client and sHsp trajectory dataframes
#have to take these from the actual raw data folders because somehow when they moved from here to the 'imageJ' results folders they stopped matching to the same image? 
# client_trajectories=pd.read_csv(f'{input_folder}20220311_exp49_comlpexes_clic_ab_t20/channel488 and 647 dual ex_seq0000 (2)_#1/Clic_colocal_traj.csv')
# client_trajectories.drop([col for col in client_trajectories.columns.tolist() if ' ' in col], axis=1, inplace=True)
# hsp_trajectories=pd.read_csv(f'{input_folder}20220311_exp49_comlpexes_clic_ab_t20/channel488 and 647 dual ex_seq0000 (2)_#1/Hsp_colocal_traj.csv')
# hsp_trajectories.drop([col for col in hsp_trajectories.columns.tolist() if ' ' in col], axis=1, inplace=True)
# #read in co-ordinates data ('client colocalisation' file)
# coords_data= pd.read_csv(f'raw_data/Experiment_49/20220311_experiment49_clic_aBc/20220311_exp49_comlpexes_clic_ab_t20/channel488 and 647 dual ex_seq0000 (2)_#1/Clic_colocalisation.csv')
# coords_data.drop([col for col in coords_data.columns.tolist() if ' ' in col], axis=1, inplace=True)
client='Clic'
hsp='Hsp'

#import the step sizes data that was gathered during trajectory analysis script
step_sizes_hsp=pd.read_csv(f'{python_input}/aB-c/fitting_changepoints/median_steps.csv')
step_sizes_hsp=step_sizes_hsp.drop([col for col in step_sizes_hsp.columns.tolist() if 'Unnamed: 0' in col],axis=1)
step_sizes_client=pd.read_csv(f'{python_input}/CLIC/fitting_changepoints/median_steps.csv')
step_sizes_client=step_sizes_client.drop([col for col in step_sizes_client.columns.tolist() if 'Unnamed: 0' in col],axis=1)

#input_folder ='raw_data/20220530_exp56_flucaBc/'
#specify these at the start

timepoint_folders=[folder for folder in os.listdir(f'{input_folder}')if 'Beads' not in folder]
image_folders=[]

def grab_trajectories_paths_only(input_folder, timepoint_folders):
    hsp_trajectories_paths=[]
    client_trajectories_pahts=[]
    for tp_folder in timepoint_folders:
        subfolders=[folder for folder in os.listdir(f'{input_folder}{tp_folder}/')]
        if 'Trajectories' not in subfolders:
            new_input=f'{input_folder}{tp_folder}/647/'
        else: 
            new_input=f'{input_folder}{tp_folder}'    
                
        filenames = next(os.walk(f'{new_input}'), (None, None, []))[2]
        image_folder_list=[folder for folder in os.listdir(f'{new_input}') if folder not in filenames]
        image_folder_list=[folder for folder in image_folder_list if 'Trajectories' not in folder]

        image_folders.append(image_folder_list)
        #files=[]
        yikes=[[f'{root}/{name}' for name in files if '_colocal_traj.csv' in name]for root, dirs, files in os.walk(f'{new_input}/')]
        yikes=[item for sublist in yikes for item in sublist]
        yikes=[item for item in yikes if 'Trajectories' not in item]
        yikes=[item for item in yikes if item not in filenames]
        clients=[path for path in yikes if client in path]

        hsp_trajs=[path for path in yikes if hsp in path]
        hsp_trajectories_paths.append(hsp_trajs)
        client_trajectories_pahts.append(clients)



    return hsp_trajectories_paths, client_trajectories_pahts



def read_trajectories(hsp_trajs, clients):
    new=[]
    for trajectory in hsp_trajs:
        hsp_trajectories=pd.read_csv(trajectory)
        hsp_trajectories.columns = [str(col) + '_hsp' for col in hsp_trajectories.columns]
        for client_traj in clients:
            if client_traj.split('/')[-2]==trajectory.split('/')[-2]:
                client_trajectories=pd.read_csv(client_traj)
                client_trajectories.columns = [str(col) + '_client' for col in client_trajectories.columns]
        
        
        for item in hsp_trajectories:
            trajectory_hsp=hsp_trajectories[[item]]
            hspnumber=item.split('_')[1]
            for item1 in client_trajectories:
                client_trajectory=client_trajectories[[item1]]
                clientnumber=item1.split('_')[1]
                if hspnumber == clientnumber:
                    newtraj=trajectory_hsp.join(client_trajectory)
                    newtraj=newtraj.T
                    newtraj['path_info']=trajectory
                    newtraj['identifier']=random.randint(1,100000)
                    if '' in trajectory.split('/'):
                        newtraj['timepoint (min)']=trajectory.split('/')[-5].split('_')[-1]
                    else: 
                         newtraj['timepoint (min)']=trajectory.split('/')[-3].split('_')[-1]
                    newtraj=newtraj.reset_index().rename(columns={'index':'molecule_number'})
                    new.append(newtraj)
    new=pd.concat(new)
    new['molecule_number']=[f'{molecule_number}_{x}' for x, molecule_number in enumerate(new['molecule_number'])]
    return new



def matched_counts(new, step_sizes_hsp, step_sizes_client):
    molecule_counts = []
    for (molecule,timepoint), df in new.groupby(['molecule_number', 'timepoint (min)']):
        if 'hsp' in molecule:
            # coords_hsp=df['coords'].values[0]    
            hsp_max_fluorescence_value = df[timepoint_columns].iloc[:,0:3].values.mean()
            # Calculate average number of molecules by mean fluorescence / step size
            molecule_count_hsp=hsp_max_fluorescence_value/int(step_sizes_hsp['step_size'].values[step_sizes_hsp['step_type']=='last_step'])
            ID=df['identifier']
            molecule_counts.append(pd.DataFrame([molecule, hsp_max_fluorescence_value, molecule_count_hsp, int(timepoint), df['path_info'].values, ID[0]]).T)
            #max_fluorescence_value = np.max(sorted(df[timepoint_columns].values[0], reverse=True))

        if 'client' in molecule:
            # coords_client=df['coords'].values[0]        
            client_max_fluorescence_value = df[timepoint_columns].iloc[:,0:3].values.mean()
            molecule_count_client=client_max_fluorescence_value/int(step_sizes_client['step_size'].values[step_sizes_client['step_type']=='last_step'])
            ID=df['identifier'].values.tolist()
            ID=ID[0]
            # Calculate average number of molecules by mean fluorescence / step size
            molecule_counts.append(pd.DataFrame([molecule, client_max_fluorescence_value, molecule_count_client, int(timepoint), df['path_info'].values, ID]).T)


    molecule_counts = pd.concat(molecule_counts)
    molecule_counts.columns = ['molecule_number', 'max_fluorescence', 'molecule_count', 'timepoint (min)', 'path_info', 'UNIQUE_ID']
    return molecule_counts



hsp_trajs, clients = grab_trajectories_paths_only(input_folder=input_folder, timepoint_folders=timepoint_folders)


hsp_trajs=[item for sublist in hsp_trajs for item in sublist]
clients=[item for sublist in clients for item in sublist]

new=read_trajectories(hsp_trajs,clients)

timepoint_columns = [col for col in new.columns.tolist() if col not in ['molecule_number','path_info','timepoint (min)', 'identifier']]
timepoint_map={'t0':0, 't20':20, 't40':40, 't4h':240, 't7h': 420}
new['timepoint (min)']=new['timepoint (min)'].map(timepoint_map)

molecule_counts=matched_counts(new, step_sizes_hsp, step_sizes_client)


#pivot the molecule counts data frame so that everything is grouped correctly for plotting (and in longform for seaborn)
molecule_counts['molecule_type']= molecule_counts['molecule_number'].str.split('_').str[-2]
#molecule_counts['molecule']= molecule_counts['molecule_number'].str.split('_').str[-3]
molecule_counts.to_csv(f'{output_folder}molecule_counts_sorted.csv')
#molecule_counts['molecule_count']=molecule_counts['molecule_count'].astype(int)


molecule_pivot=pd.pivot(molecule_counts,index=['timepoint (min)', 'UNIQUE_ID'], columns='molecule_type', values='molecule_count').reset_index()
molecule_pivot['client']=molecule_pivot['client'].astype(float)
molecule_pivot['hsp']=molecule_pivot['hsp'].astype(float)
#when I turned the counts into integers, the smaller step ones turned into zeros, so I need to make them ones in order to plot this
#molecule_pivot['hsp']=molecule_pivot['hsp'].mask(molecule_pivot['hsp']==0).fillna(1)
molecule_pivot.to_csv(f'{output_folder}molecule_counts_for_plotting.csv')

#plotting using a hex bin plot for each timepoint
fig, axes = plt.subplots(ncols=5, sharey=True, figsize=(18, 3))
for x, (timepoint, df) in enumerate(molecule_pivot.groupby('timepoint (min)')):

    hb = axes[x].hexbin(x=df['client'], y=df['hsp'], gridsize=int(molecule_pivot['client'].max()),cmap='crest', vmin=0, vmax=2)
    axes[x].set(xlim=(0,30), ylim=(0,20))
    axes[x].set_title(f'{timepoint} (min)')


#molecule_pivot['hsp'].max()
#molecule_pivot['client'].max()
    #sns.jointplot(data=molecule_pivot, x='client', y='hsp', kind="hex", color="#4CB391")
cb = fig.colorbar(hb, ax=axes[4], label='counts')
plt.show()

#plotting a kernel density plot
fig, axes = plt.subplots(ncols=5, sharey=True, figsize=(18, 3))
for x, (timepoint, df) in enumerate(molecule_pivot.groupby('timepoint (min)')):
    sns.kdeplot(ax=axes[x], x=df["client"], y=df["hsp"], shade=True, palette='Greens')
    #hb = axes[x].hexbin(x=df['client'], y=df['hsp'], gridsize=int(molecule_pivot['client'].max()),cmap='Blues', vmin=0, vmax=3)
    axes[x].set(xlim=(0,molecule_pivot['client'].max()), ylim=(0,molecule_pivot['hsp'].max()))
    axes[x].set_title(f'{timepoint} (min)')
    #axes[x].get_legend().remove()



    #sns.jointplot(data=molecule_pivot, x='client', y='hsp', kind="hex", color="#4CB391")
#cb = fig.colorbar(hb, ax=axes[5], label='counts')
#axes[5].legend(molecule_pivot['timepoint (min)'].unique().tolist(), loc='upper right', ncol=2)
plt.show()
plt.savefig(f'{output_folder}CLIC_aB-c_kde_plots_timepoints.png')
kde=sns.kdeplot(x=molecule_pivot["client"], y=molecule_pivot["hsp"], hue=molecule_pivot['timepoint (min)'], palette='Blues')
figure=kde.get_figure()
figure.savefig(f'{output_folder}combined_kde_greens.png')
sns.kdeplot(x=molecule_pivot["client"], y=molecule_pivot["hsp"], shade=True)
#can change palette of the above using palette not colourmap