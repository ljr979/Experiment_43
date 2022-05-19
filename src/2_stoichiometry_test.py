import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

input_folder = 'raw_data/Experiment_49/20220311_experiment49_clic_aBc/'
output_folder = 'python_results/stoichiometry_test/'

if not os.path.exists(output_folder):
        os.makedirs(output_folder)

#first, need to go into the trajectories folders, and grab the COLOCALISED trajectories. maybe need to do this file by file before doing all of them
#assign a variable to client and sHsp trajectory dataframes
#have to take these from the actual raw data folders because somehow when they moved from here to the 'imageJ' results folders they stopped matching to the same image? 
client_trajectories=pd.read_csv(f'{input_folder}20220311_exp49_comlpexes_clic_ab_t20/channel488 and 647 dual ex_seq0000 (2)_#1/Clic_colocal_traj.csv')
client_trajectories.drop([col for col in client_trajectories.columns.tolist() if ' ' in col], axis=1, inplace=True)
hsp_trajectories=pd.read_csv(f'{input_folder}20220311_exp49_comlpexes_clic_ab_t20/channel488 and 647 dual ex_seq0000 (2)_#1/Hsp_colocal_traj.csv')
hsp_trajectories.drop([col for col in hsp_trajectories.columns.tolist() if ' ' in col], axis=1, inplace=True)
#read in co-ordinates data ('client colocalisation' file)
coords_data= pd.read_csv(f'raw_data/Experiment_49/20220311_experiment49_clic_aBc/20220311_exp49_comlpexes_clic_ab_t20/channel488 and 647 dual ex_seq0000 (2)_#1/Clic_colocalisation.csv')
coords_data.drop([col for col in coords_data.columns.tolist() if ' ' in col], axis=1, inplace=True)


#import the step sizes data that was gathered during trajectory analysis script
step_sizes_hsp=pd.read_csv('python_results/experiment_49/stoichiometries/complexes/aB-c/fitting_changepoints/median_steps.csv')
step_sizes_hsp=step_sizes_hsp.drop([col for col in step_sizes_hsp.columns.tolist() if 'Unnamed: 0' in col],axis=1)
step_sizes_client=pd.read_csv('python_results/experiment_49/stoichiometries/complexes/CLIC/fitting_changepoints/median_steps.csv')
step_sizes_client=step_sizes_client.drop([col for col in step_sizes_client.columns.tolist() if 'Unnamed: 0' in col],axis=1)

input_folder ='raw_data/Experiment_49/20220311_experiment49_clic_aBc/'
#specify these at the start
client='Clic'
hsp='Hsp'
timepoint_folders=[folder for folder in os.listdir(f'{input_folder}')if 'beads' not in folder]
image_folders=[]

for tp_folder in timepoint_folders:
    filenames = next(os.walk(f'{input_folder}{tp_folder}/'), (None, None, []))[2]
    image_folder_list=[folder for folder in os.listdir(f'{input_folder}{tp_folder}/') if folder not in filenames]
    image_folder_list=[folder for folder in image_folder_list if 'Trajectories' not in folder]
    image_folders.append(image_folder_list)
    #files=[]
    yikes=[[f'{root}/{name}' for name in files if '_colocal_traj.csv' in name]for root, dirs, files in os.walk(f'{input_folder}{tp_folder}/')]
    yikes=[item for sublist in yikes for item in sublist]
    yikes=[item for item in yikes if 'Trajectories' not in item]
    yikes=[item for item in yikes if item not in filenames]
    clients=[path for path in yikes if client in path]
    hsp_trajs=[path for path in yikes if hsp in path]
    #if I want them all in one big list , put files list outside for loop, and append here, then flatten in the loop
    #files.append(yikes)
    #files=[item for sublist in files for item in sublist]
    trajectories=[]
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
                    newtraj['path_info']=trajec_test=trajectory.split('/')[-3:-1]
                    newtraj=newtraj.reset_index().rename(columns={'index':'molecule_number'})
                    new.append(newtraj)
                    
        
        timepoint_columns = [col for col in newtraj.columns.tolist() if col not in ['molecule_number','path_info']]

        
        molecule_counts = []
        for molecule, df in newtraj.groupby('molecule_number'):
            if 'hsp' in molecule:
                # coords_hsp=df['coords'].values[0]    
                hsp_max_fluorescence_value = df[timepoint_columns].iloc[:,0:3].values.mean()
                # Calculate average number of molecules by mean fluorescence / step size
                molecule_count_hsp=hsp_max_fluorescence_value/int(step_sizes_hsp['step_size'].values[step_sizes_hsp['step_type']=='last_step'])
                molecule_counts.append(pd.DataFrame([molecule, hsp_max_fluorescence_value, molecule_count_hsp]).T)
                #max_fluorescence_value = np.max(sorted(df[timepoint_columns].values[0], reverse=True))

            if 'client' in molecule:
                # coords_client=df['coords'].values[0]        
                client_max_fluorescence_value = df[timepoint_columns].iloc[:,0:3].values.mean()
                molecule_count_client=client_max_fluorescence_value/int(step_sizes_client['step_size'].values[step_sizes_client['step_type']=='last_step'])
                # Calculate average number of molecules by mean fluorescence / step size
                molecule_counts.append([molecule, client_max_fluorescence_value, molecule_count_client])
                
        molecule_counts = pd.concat(molecule_counts)
        molecule_counts.columns = ['molecule_number', 'max_fluorescence', 'molecule_count']

                


    
    # newtraj=pd.concat(newtraj).reset_index().rename(columns={'index':'molecule_number'})


#grab the first column in each trajectories dataframe, and assign them one molecule name/number (perhaps their xy co-ordinates? + metadata). 


hsp_trajectories.columns = [str(col) + '_hsp' for col in hsp_trajectories.columns]
client_trajectories.columns = [str(col) + '_client' for col in client_trajectories.columns]
hsp_cols = hsp_trajectories.columns
client_cols = client_trajectories.columns.tolist()


#for traj_name in hsp_trajectories:
new=[]
for item in hsp_trajectories:
    trajectory_hsp=hsp_trajectories[[item]]
    hspnumber=item.split('_')[1]
    coords_hsp = coords_data.loc[int(hspnumber),['X2', 'Y2']]
    coords_hsp=(coords_hsp[0], coords_hsp[1])
    for item1 in client_trajectories:
        client_trajectory=client_trajectories[[item1]]
        clientnumber=item1.split('_')[1]
        coords_client = coords_data.loc[int(clientnumber),['X', 'Y']]
        coords_client=(coords_client[0], coords_client[1])
        if hspnumber == clientnumber:
            newtraj=trajectory_hsp.join(client_trajectory)
            newtraj=newtraj.T
            newtraj['coords']=[coords_hsp, coords_client]
            new.append(newtraj)
            
new=pd.concat(new).reset_index().rename(columns={'index':'molecule_number'})

# combine these two trajectories (hsp and client, matched) into one dataframe

# new=hsp_trajectories[['Mean_0_hsp']]
# new['Mean_0_client']=client_trajectories[['Mean_0_client']]

# new=new.T.reset_index().rename(columns={'index':'molecule_number'})
# coords=[]
# coords_client=pd.DataFrame(coords_data.loc[0,['X', 'Y']]).T
# coords_client=(coords_client.X[0], coords_client.Y[0])



# coords_hsp=pd.DataFrame(coords_data.loc[0,['X2', 'Y2']]).T
# coords_hsp=(coords_hsp.X2[0], coords_hsp.Y2[0])
# coords.append(coords_hsp)
# coords.append(coords_client)
# new['co_ords']=coords


#run analysis function 'calculating stoichiometries' on them to find the molecule size. 
timepoint_columns = [col for col in new.columns.tolist() if col not in ['molecule_number','co_ords']]


molecule_counts = []
for molecule, df in new.groupby('molecule_number'):
    if 'hsp' in molecule:
        coords_hsp=df['coords'].values[0]    
        hsp_max_fluorescence_value = df[timepoint_columns].iloc[:,0:3].values.mean()
        # Calculate average number of molecules by mean fluorescence / step size
        molecule_count_hsp=hsp_max_fluorescence_value/int(step_sizes_hsp['step_size'].values[step_sizes_hsp['step_type']=='last_step'])
        molecule_counts.append(pd.DataFrame([molecule, hsp_max_fluorescence_value, coords_hsp, molecule_count_hsp]).T)
        #max_fluorescence_value = np.max(sorted(df[timepoint_columns].values[0], reverse=True))

    if 'client' in molecule:
        coords_client=df['coords'].values[0]        
        client_max_fluorescence_value = df[timepoint_columns].iloc[:,0:3].values.mean()
        molecule_count_client=client_max_fluorescence_value/int(step_sizes_client['step_size'].values[step_sizes_client['step_type']=='last_step'])
        # Calculate average number of molecules by mean fluorescence / step size
        molecule_counts.append(pd.DataFrame([molecule, client_max_fluorescence_value, coords_client, molecule_count_client]).T)


molecule_counts = pd.concat(molecule_counts)
molecule_counts.columns = ['molecule_number', 'max_fluorescence', 'co_ords', 'molecule_count']

#plot the above dataframe, 
