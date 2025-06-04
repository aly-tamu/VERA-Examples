import h5py
import numpy as np
import matplotlib.pyplot as plt



def hdf5_to_dict(filepath, with_attributes=False):
    '''
    Get any HDF5 and turn it into a python dictionary.
    
    Parameters
    ----------
    filepath : str
        Path to HDF5 file.
    with_attributes : bool, optional
        Include HDF5 attributes in the python diction. The default is False.

    Returns
    -------
    h5_dict : dict
        HDF5 converted into python dictionary.
    '''
   
    with h5py.File(str(filepath), 'r') as h5f:
        h5_dict = {}
        stack = [(h5f, h5_dict)] # Stack holds tuples of (HDF5 object, dictionary)
        while stack:
            current_group, current_dict = stack.pop()   # Remove values from bottom of stack
        
            for key, value in current_group.items():    # If it's a group, iterate over its items
            
                if with_attributes == True:             # If you need the h5 attributes
                    current_dict['attributes'] = {attr: current_group.attrs[attr] for attr in current_group.attrs}
                
                if isinstance(value, h5py.Group):       # If it's a subgroup, add to the stack to process later
                    current_dict[key] = {}
                    stack.append((value, current_dict[key]))
                    
                elif isinstance(value, h5py.Dataset):   # If it's a dataset, add to the dictionary
                    current_dict[key] = value[()]
    return h5_dict
    
def get_mgxs_from_h5(mgxs_file, groups=None): 
    '''
    Takes mgxs.h5 file from OpenMC and turns it into a python dictionary, 
    then reformats the cross section data to be plottable.

    Parameters
    ----------
    mgxs_file : str
        Path to mgxs.h5 file.
    groups : list, optional
        Specify which sets to include in the dictionary.
        Formatted as a list of set.
        i.e ['desired set']

    Returns
    -------
    mgxs_dict : dict
        Reformatted cross section dictionary.

    '''
    attr_dict = hdf5_to_dict(mgxs_file,True)    # attributes to get the scattering order from
    h5_dict = hdf5_to_dict(mgxs_file,False)     # get mgxs hdf5 into dictionary form
    mgxs_dict = {}                              # Will become our neater, reformatted dictionary
    
    # Search and destroy
    if groups:
        for key in list(h5_dict.keys()):
            if key not in groups: del h5_dict[key]                    
    
    for set_n in h5_dict:
        mgxs_dict[set_n] = {} 
        for temp in h5_dict[set_n]:
            if temp != 'kTs':
                mgxs_dict[set_n][temp] = {}
                for XS_type in h5_dict[set_n][temp]:    # Loop through all sets, temps, and cross section types
                    
                    if isinstance(h5_dict[set_n][temp][XS_type],dict):    # If XS is matrix then it needs reformatting
                        g_max = h5_dict[set_n][temp][XS_type]['g_max']
                        g_min = h5_dict[set_n][temp][XS_type]['g_min']
                        
                        for XS_mat in h5_dict[set_n][temp][XS_type]:    # Loop through all XS matrices under label
                            if XS_mat != 'g_max' and XS_mat != 'g_min': # skip g_min and g_max data
                                
                                xy_dim = g_max.shape[0]
                                # matrix size for scatter matrices is dependent on the legendre order
                                z_dim = attr_dict[set_n]['attributes']['order']+1 if XS_mat != 'multiplicity_matrix' else 1
                                mgxs_dict[set_n][temp][XS_mat] = np.zeros([xy_dim,xy_dim,z_dim])
                                # Build matrix in neater dictionary
                                list_ind = 0
                                for i in range(0,xy_dim):
                                    for j in range(g_min[i]-1,g_max[i]):
                                        for k in range(0,z_dim):
                                            mgxs_dict[set_n][temp][XS_mat][i,j,k] =  h5_dict[set_n][temp][XS_type][XS_mat][list_ind]
                                            list_ind += 1
                                
                    else:   # Else just copy cross section vector to neater matrix 
                        mgxs_dict[set_n][temp][XS_type] = h5_dict[set_n][temp][XS_type]
                        
    return mgxs_dict         
                    
def plot_mgxs(mgxs_file, groups=None):
    
    # Get energy groups from attributes
    attr_dict = hdf5_to_dict(mgxs_file, True)
    energy_groups = attr_dict['attributes']['group structure']
    energy_mids = np.zeros(energy_groups.shape[0]-1)
    for i in range(1,energy_mids.shape[0]):
        energy_mids[i] = (energy_groups[i+1] + energy_groups[i])/2
        
    # Get clean cross section data dictionary
    mgxs_dict = get_mgxs_from_h5(mgxs_file, groups)
    for set_n in mgxs_dict:
        for temp in mgxs_dict[set_n]:
            for XS_type in mgxs_dict[set_n][temp]:  # loop through all sets, temps, and cross sections to plot
                XS = mgxs_dict[set_n][temp][XS_type]
                
                if len(XS.shape) > 1: # If 2-D array
                    plt.figure(figsize=(15,10))
                    plt.imshow(XS[::-1,::-1,0], cmap='turbo', interpolation='nearest')
                    plt.colorbar()
                    plt.title(set_n+' '+XS_type)
                    plt.xlabel('from group')
                    plt.ylabel('to group')
                    plt.show()
                    
                else: # If 1-D array
                    plt.figure(figsize=(15,10))
                    plt.loglog(energy_mids, XS[::-1])
                    plt.title(set_n+' '+XS_type)
                    plt.xlabel('Energy, eV')
                    plt.ylabel('Cross Section, 1/cm')
                    plt.grid(True)
                    plt.show()