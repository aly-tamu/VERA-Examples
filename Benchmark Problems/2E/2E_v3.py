#!/usr/bin/env python
# coding: utf-8

# # Fuel assembly: 2E

# ## Import modules

# In[1]:


import sys
import os
import pathlib
import re

import openmc
import openmc.mgxs as mgxs

import numpy as np
import matplotlib.pyplot as plt
import csv
from collections import Counter
from prettytable import PrettyTable


# ## Add location of OpenMC XS executable, setup ENDF xs path

# In[2]:


# Add path to OpenMC binary
# os.environ['PATH'] += r':/path/to/openmc/bin'
#os.environ['PATH'] += r':/Users/jean.ragusa/repo/openmc/local/bin'

# Add location of OpenMC xs data
#%env OPENMC_CROSS_SECTIONS=/Users/jean.ragusa/repo/endfb-vii.1-hdf5/cross_sections.xml

# Are you ragusa?
your_files = os.getcwd()
print(your_files)
if 'ragusa' in your_files:
    os.environ['OPENMC_CROSS_SECTIONS'] = '/home/ragusa/xs/endfb-viii.0-hdf5/cross_sections.xml'


# # Start model

# In[3]:


model=openmc.Model()


# ## Define Materials
# ***

# In[4]:


uo2 = openmc.Material(name='uo2')

uo2.add_nuclide('U234', 6.11864E-06, 'ao')
uo2.add_nuclide('U235', 7.18132E-04, 'ao')
uo2.add_nuclide('U236', 3.29861E-06, 'ao')
uo2.add_nuclide('U238', 2.21546E-02, 'ao')
uo2.add_nuclide('O16', 4.57642E-02, 'ao')

uo2.set_density('g/cm3', 10.257 )

uo2.temperature = 600


# In[5]:


zirconium = openmc.Material(name="zirconium")

zirconium.add_nuclide('Zr90',  2.18865E-02, 'ao')
zirconium.add_nuclide('Zr91',  4.77292E-03, 'ao')
zirconium.add_nuclide('Zr92',  7.29551E-03, 'ao')
zirconium.add_nuclide('Zr94',  7.39335E-03, 'ao')
zirconium.add_nuclide('Zr96',  1.19110E-03, 'ao')
zirconium.add_nuclide('Sn112', 4.68066E-06, 'ao')
zirconium.add_nuclide('Sn114', 3.18478E-06, 'ao')
zirconium.add_nuclide('Sn115', 1.64064E-06, 'ao')
zirconium.add_nuclide('Sn116', 7.01616E-05, 'ao')
zirconium.add_nuclide('Sn117', 3.70592E-05, 'ao')
zirconium.add_nuclide('Sn118', 1.16872E-04, 'ao')
zirconium.add_nuclide('Sn119', 4.14504E-05, 'ao')
zirconium.add_nuclide('Sn120', 1.57212E-04, 'ao')
zirconium.add_nuclide('Sn122', 2.23417E-05, 'ao')
zirconium.add_nuclide('Sn124', 2.79392E-05, 'ao')
zirconium.add_nuclide('Fe54',  8.68307E-06, 'ao')
zirconium.add_nuclide('Fe56',  1.36306E-04, 'ao')
zirconium.add_nuclide('Fe57',  3.14789E-06, 'ao')
zirconium.add_nuclide('Fe58',  4.18926E-07, 'ao')
zirconium.add_nuclide('Cr50',  3.30121E-06, 'ao')
zirconium.add_nuclide('Cr52',  6.36606E-05, 'ao')
zirconium.add_nuclide('Cr53',  7.21860E-06, 'ao')
zirconium.add_nuclide('Cr54',  1.79686E-06, 'ao')
zirconium.add_nuclide('Hf174', 3.54138E-09, 'ao')
zirconium.add_nuclide('Hf176', 1.16423E-07, 'ao')
zirconium.add_nuclide('Hf177', 4.11686E-07, 'ao')
zirconium.add_nuclide('Hf178', 6.03806E-07, 'ao')
zirconium.add_nuclide('Hf179', 3.01460E-07, 'ao')
zirconium.add_nuclide('Hf180', 7.76449E-07, 'ao')

zirconium.set_density('g/cm3',  6.56)

zirconium.temperature = 600


# In[6]:


water = openmc.Material(name="water")

water.add_nuclide('H1',  4.96224E-02, 'ao')
water.add_nuclide('O16', 2.48112E-02, 'ao')
water.add_nuclide('B10', 1.07070E-05, 'ao')
water.add_nuclide('B11', 4.30971E-05, 'ao')

water.add_s_alpha_beta('c_H_in_H2O')

water.set_density('g/cm3', 0.743)

water.temperature=600


# In[7]:


helium = openmc.Material(name="helium")

helium.add_nuclide('He4',1, 'ao')

helium.set_density('g/cm3', 0.178E-03 )

helium.temperature=600


# In[8]:


pyrex_mat = openmc.Material(name="pyrex_mat")

pyrex_mat.add_nuclide('B10',9.63266E-04, 'ao')
pyrex_mat.add_nuclide('B11',3.90172E-03, 'ao')
pyrex_mat.add_nuclide('O16',4.67761E-02, 'ao')
pyrex_mat.add_nuclide('Si28',1.81980E-02, 'ao')
pyrex_mat.add_nuclide('Si29',9.24474E-04, 'ao')
pyrex_mat.add_nuclide('Si30',6.10133E-04, 'ao')

pyrex_mat.set_density('g/cm3', 2.25 )

pyrex_mat.temperature=600


# In[9]:


SS304 = openmc.Material(name="SS304")

SS304.add_nuclide('C0',3.20895E-04, 'ao')
SS304.add_nuclide('Si28',1.58197E-03, 'ao')
SS304.add_nuclide('Si29',8.03653E-05, 'ao')
SS304.add_nuclide('Si30',5.30394E-05, 'ao')
SS304.add_nuclide('P31',6.99938E-05, 'ao')
SS304.add_nuclide('Cr50',7.64915E-04, 'ao')
SS304.add_nuclide('Cr52',1.47506E-02, 'ao')
SS304.add_nuclide('Cr53',1.67260E-03, 'ao')
SS304.add_nuclide('Cr54',4.16346E-04, 'ao')
SS304.add_nuclide('Mn55',1.75387E-03, 'ao')
SS304.add_nuclide('Fe54',3.44776E-03, 'ao')
SS304.add_nuclide('Fe56',5.41225E-02, 'ao')
SS304.add_nuclide('Fe57',1.24992E-03, 'ao')
SS304.add_nuclide('Fe58',1.66342E-04, 'ao')
SS304.add_nuclide('Ni58',5.30854E-03, 'ao')
SS304.add_nuclide('Ni60',2.04484E-03, 'ao')
SS304.add_nuclide('Ni61',8.88879E-05, 'ao')
SS304.add_nuclide('Ni62',2.83413E-04, 'ao')
SS304.add_nuclide('Ni64',7.21770E-05, 'ao')

SS304.set_density('g/cm3', 8.00 )

SS304.temperature=600


# In[10]:


model.materials = openmc.Materials([uo2, zirconium, water,helium,pyrex_mat,SS304]) 
#print(model.materials)


# ## Define individiual cells
# ***

# In[11]:


# Global surfaces
pitch = 1.26

# Fuel radii
fuel_outer_radius = openmc.ZCylinder(r=0.4096)
clad_inner_radius = openmc.ZCylinder(r=0.418)
clad_outer_radius = openmc.ZCylinder(r=0.475)

# Guide tube radii
gt_inner_radius = openmc.ZCylinder(r=0.561)
gt_outer_radius = openmc.ZCylinder(r=0.602)

# Instrumentation radii
it_inner_radius = openmc.ZCylinder(r=0.559)
it_outer_radius = openmc.ZCylinder(r=0.605)

#Pyrex
Inner_Tube_Inner_Radius= openmc.ZCylinder(r=0.214)
Inner_Tube_Outer_Radius= openmc.ZCylinder(r= 0.231)
Pyrex_Inner_Radius= openmc.ZCylinder(r= 0.241)
Pyrex_Outer_Radius= openmc.ZCylinder(r= 0.427)
Outer_Tube_Inner_Radius= openmc.ZCylinder(r= 0.437)
Outer_Tube_Outer_Radius= openmc.ZCylinder(r= 0.484)

# Universe bounds
left   = openmc.XPlane(-pitch/2, boundary_type='transmission')
right  = openmc.XPlane( pitch/2, boundary_type='transmission')
bottom = openmc.YPlane(-pitch/2, boundary_type='transmission')
top    = openmc.YPlane( pitch/2, boundary_type='transmission')


# ### Define pincell universe

# In[12]:


def pincell():

    fuel_region  = -fuel_outer_radius
    gap_region   = +fuel_outer_radius & -clad_inner_radius
    clad_region  = +clad_inner_radius & -clad_outer_radius
    water_region = +left & -right & +bottom & -top & +clad_outer_radius

    fuel = openmc.Cell(name='fuel')
    fuel.region = fuel_region
    fuel.fill = uo2
    #new_fuel = uo2.clone()
    #fuel.fill = new_fuel

    gap = openmc.Cell(name='gap')
    gap.region = gap_region
    gap.fill = helium     

    clad = openmc.Cell(name='clad')
    clad.region = clad_region
    clad.fill = zirconium

    moderator = openmc.Cell(name='moderator')
    moderator.region = water_region
    moderator.fill = water

    return openmc.Universe(name='pin-univ', cells=(fuel, clad, moderator,gap))


# ### Create guide tube universe

# In[13]:


def guide():

    gt_water_in_region  = -gt_inner_radius
    gt_clad_region      = +gt_inner_radius & -gt_outer_radius
    gt_water_out_region = +gt_outer_radius & +left & -right & +bottom & -top

    gt_water_in = openmc.Cell(name='gt-water-in')
    gt_water_in.region = gt_water_in_region
    gt_water_in.fill = water

    gt_clad = openmc.Cell(name='gt-clad')
    gt_clad.region = gt_clad_region
    gt_clad.fill = zirconium

    gt_water_out = openmc.Cell(name='gt-water-out')
    gt_water_out.region = gt_water_out_region
    gt_water_out.fill = water

    return openmc.Universe(name='gt-univ', cells=(gt_water_in, gt_clad, gt_water_out))


# ### Create instrumentation tube universe

# In[14]:


def instrument():

    it_water_in_region  = -it_inner_radius
    it_clad_region      = +it_inner_radius & -it_outer_radius
    it_water_out_region = +it_outer_radius & +left & -right & +bottom & -top

    it_water_in = openmc.Cell(name='it-water-in')
    it_water_in.region = it_water_in_region
    it_water_in.fill = water

    it_clad = openmc.Cell(name='it-clad')
    it_clad.region = it_clad_region
    it_clad.fill = zirconium

    it_water_out = openmc.Cell(name='it-water-out')
    it_water_out.region = it_water_out_region
    it_water_out.fill = water

    return openmc.Universe(name='it-univ', cells=(it_water_in, it_clad, it_water_out))


# In[15]:


def pyrex_pin():

    pyrex_region= +Pyrex_Inner_Radius & -Pyrex_Outer_Radius
    water_region = +left & -right & +bottom & -top & +gt_outer_radius | +Outer_Tube_Outer_Radius & -gt_inner_radius 
    guide_region = +gt_inner_radius & -gt_outer_radius
    steel_region = +Outer_Tube_Inner_Radius & -Outer_Tube_Outer_Radius | +Inner_Tube_Inner_Radius & -Inner_Tube_Outer_Radius
    gap_region = +Inner_Tube_Outer_Radius & -Pyrex_Inner_Radius |  +Pyrex_Outer_Radius & - Outer_Tube_Inner_Radius | -Inner_Tube_Inner_Radius

    steel_cell = openmc.Cell(name='steel_cell')
    steel_cell.fill = SS304
    steel_cell.region = steel_region  

    gap_cell = openmc.Cell(name='gap_cell')
    gap_cell.region = gap_region 
    gap_cell.fill = helium

    pyrex_cell = openmc.Cell(name='pyrex_cell')
    pyrex_cell.fill = pyrex_mat
    pyrex_cell.region = pyrex_region  

    water_cell = openmc.Cell(name='water_cell')
    water_cell.fill = water
    water_cell.region = water_region

    guide_cell = openmc.Cell(name='guide_cell')
    guide_cell.fill = zirconium
    guide_cell.region = guide_region

    return openmc.Universe(name='py-univ', cells=(steel_cell,gap_cell,pyrex_cell,water_cell,guide_cell))


# ### Generate universes

# In[16]:


pc_univ  = pincell()
gt_univ  = guide()
it_univ  = instrument()
py_univ = pyrex_pin()


# ## Generate assembly
# ***

# In[17]:


height = 10 # height in z-dir
dr = 0.04   # cm of water that is outside assembly
one_eighth = True # either 1/8 or 1/4 of FA


# In[18]:


def read_csv_to_2d_array(file_path):
    with open(file_path, newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        data = [row for row in reader]
    return np.asarray(data)

def count_frequencies(data):
    flattened_data = [item for row in data for item in row]  # Flatten 2D array into a 1D list
    cell_frequencies = Counter(flattened_data)
    print("cell name frequency:")
    total = 0
    for key, value in cell_frequencies.items():
        print(f'"{key}": {value}')
        total += value
    print("total: ",total)

csv_filepath = 'FA_cell_names_1_family.csv'

lattice_csv = read_csv_to_2d_array(csv_filepath)
count_frequencies(lattice_csv)

if lattice_csv.shape[0] != lattice_csv.shape[1]:
    raise Exception('CSV array of cell names is not square.')

size = lattice_csv.shape[0] #size of the assembly 


# ### Define assembly lattice

# In[19]:


assembly = openmc.RectLattice()
assembly.pitch = (pitch,pitch)
assembly.lower_left = (-size*pitch/2, -size*pitch/2)

lattice_array = np.empty((size, size), dtype=openmc.universe.Universe)

for row in range(0,lattice_csv.shape[0]):
    for col in range(0,lattice_csv.shape[1]):
        if lattice_csv[row][col] == 'gt': 
            lattice_array[row][col] = gt_univ
        elif lattice_csv[row][col] == 'it':
            lattice_array[row][col] = it_univ
        elif lattice_csv[row][col] == 'fu':
            lattice_array[row][col] = pc_univ
        elif lattice_csv[row][col] == 'p':
            lattice_array[row][col] = py_univ
        else:
            mesg = 'i={},j={},cell-name {} not recognized'.format(col, row, lattice_csv[row][col])
            raise ValueError(mesg)

assembly.universes = lattice_array


# ### Define moderator outside of the assembly

# In[20]:


# create cell that will contain the lattice
moderator_outside_cell = openmc.Cell( name = 'water_outside', fill = water )
assembly.outer = openmc.Universe( name='outer', cells = [moderator_outside_cell] )


# ### Define 1/8 or 1/4 of full assembly

# In[21]:


min_x = openmc.XPlane(x0= 0.               , boundary_type='reflective')
max_x = openmc.XPlane(x0= size*pitch/2+dr  , boundary_type='reflective')
min_y = openmc.YPlane(y0= 0.               , boundary_type='reflective')
max_y = openmc.YPlane(y0= size*pitch/2+dr  , boundary_type='reflective')
min_z = openmc.ZPlane(z0=-height/2         , boundary_type='reflective')
max_z = openmc.ZPlane(z0= height/2         , boundary_type='reflective')
xy    = openmc.Plane(-1., 1.               , boundary_type='reflective')

# define root universe
root_cell = openmc.Cell( name = 'root cell', fill = assembly )

if one_eighth:
    root_cell.region = -max_x & +min_y & -xy
else:
    root_cell.region = +min_x & -max_x & +min_y & -max_y & +min_z & -max_z

model.geometry.root_universe = openmc.Universe(name = 'root universe', cells=[root_cell])


# ## Plotting

# In[22]:


color_opt = 'material'
root_cell.plot(origin=(5,5,0), pixels=(1500,1500), width=(15,15),color_by=color_opt,outline=True)


# ##  Settings
# ***

# In[23]:


settings = openmc.Settings()


# ### Source

# In[24]:


bbox = openmc.BoundingBox(lower_left  = [0., 0., -height/2], 
                          upper_right = [size*pitch/2, size*pitch/2, height/2])
uniform_dist = openmc.stats.Box(bbox.lower_left, bbox.upper_right, only_fissionable=True)
source = openmc.IndependentSource(space=uniform_dist)
settings.source = source


# ### Destination path

# In[25]:


my_case = '2e'
my_path = './mgxs_casl_' + my_case

# check if folder exists
path = pathlib.Path(my_path)
path.mkdir(parents = True, exist_ok = True)


# ### Batching

# In[26]:


test_mode = False


# In[27]:


# add additional parameters
settings.batches =  150
settings.inactive = 20
settings.particles = 50000
settings.keff_trigger = {'type':'std_dev','threshold':0.00010}
settings.trigger_active = True
settings.trigger_max_batches = 50000
settings.output = {'tallies': True, 'path':my_path}
settings.temperature['method'] = 'interpolation'

if test_mode == True:
    settings.batches = 200
    settings.inactive = 50
    settings.particles = 1000
    settings.keff_trigger = {'type':'std_dev','threshold':0.001}

model.settings = settings


# ## Tallies
# ***

# In[28]:


# egroup_name = 'XMAS-172'
egroup_name = 'SHEM-361'
egroup = openmc.mgxs.GROUP_STRUCTURES[egroup_name]


# ### Power Tally
# 

# In[29]:


tally_power = openmc.Tally(name='power')

# Instantiate a tally Mesh
mesh = openmc.RegularMesh()
mesh._dimension = [size, size]
mesh._lower_left = [-size*pitch/2, -size*pitch/2]
mesh._upper_right = [size*pitch/2,  size*pitch/2]

# Instantiate tally Filter
mesh_filter = openmc.MeshFilter(mesh)
tally_power.filters = [mesh_filter]

tally_power.scores = ['fission', 'nu-fission', 'kappa-fission']

tallies = openmc.Tallies([tally_power])
# model.tallies=tallies


# ### MGXS Tally

# In[30]:


leg = 7
mgxs_domain = 'cell'


# In[31]:


# Setup MGXS lib
mgxs_lib = openmc.mgxs.Library(model.geometry)
mgxs_lib.energy_groups = openmc.mgxs.EnergyGroups(egroup)
mgxs_lib.scatter_format = "legendre"
mgxs_lib.mgxs_types = ['total', 'absorption', 'nu-fission', 'fission' ,'chi',
                       'consistent nu-scatter matrix', 'multiplicity matrix','kappa-fission']

# Legendre order
mgxs_lib.legendre_order = leg
if leg == leg:
    mgxs_lib.correction = None

# No by_nuclide
mgxs_lib.by_nuclide = False

# MGXS domain type
mgxs_lib.domain_type = mgxs_domain
if mgxs_domain == 'cell':
    mgxs_lib.domains = model.geometry.get_all_material_cells().values()
elif mgxs_domain == 'universe':
    mgxs_lib.domains = model.geometry.get_all_universes().values()

# Construct all tallies needed for the multi-group cross section library
mgxs_lib.build_library()
mgxs_lib.check_library_for_openmc_mgxs()

mgxs_lib.add_to_tallies_file(tallies, merge=True)

model.tallies = tallies


# ### Run OpenMC model

# In[32]:


# trick to make several runs work with jupyter
try:
    sp 
    print('sp found')
    sp.close()
except NameError:
    print('sp NOT found')


# In[33]:


statepoint_filename = model.run()


# In[34]:


# Load the last statepoint file
sp = openmc.StatePoint(statepoint_filename)


# ## MGXS Tally Outputs

# In[35]:


xs_names = [] 
for set in (mgxs_lib.domains):
    print(set.name)
    xs_names.append(set.name)

#print(mgxs_lib.domains)


# In[36]:


if one_eighth:
    txt = 'one_eighth'
else:
    txt = 'one_quarter'


# In[37]:


mgxs_lib.load_from_statepoint(sp)

h5_file_path = my_path + f'/mgxs_{my_case}_{txt}_{egroup_name}.h5'
print(h5_file_path)

# below, no need for xs_type = 'macro' as it is the default
mgxs_lib.create_mg_library(xsdata_names=xs_names).export_to_hdf5(h5_file_path)


# ## Power Tally Outputs

# In[38]:


computed_power_tallies = sp.get_tally()
power_tally_values = computed_power_tallies.get_values()

pin_power_file_path = my_path + f'/pinpow_{my_case}_{txt}_{egroup_name}.npy'
np.save(pin_power_file_path, power_tally_values)


# In[39]:


for score_id in range(0,len(computed_power_tallies.scores)):
     # get name of the score
    score = computed_power_tallies.scores[score_id]
    print('score = ',score)
    # extract + shape the data
    vals = power_tally_values[:,0,score_id].reshape(17,17)
    vals = vals[8:,8:]

    dd = np.copy(vals)

    # multiply axis of symeetry by 2 for power display
    dd[0,:] *= 2
    dd[np.diag_indices_from(dd)] *= 2

    dd = np.flip(dd,axis=0)

    # normalize
    idx_fuel_cells = np.argwhere(dd>0)
    nnz_fuel_cells = np.asarray(idx_fuel_cells).shape[0]
    dd /= np.sum(dd)
    dd *= nnz_fuel_cells

    table = PrettyTable()

    n_list = [str(n) for n in range(1, 10)]
    table.field_names = n_list

    # Convert ndarray to list of lists
    data_list = dd.tolist()

    for row in data_list:
        # Format each value to 2 decimal places
        formatted_row = [f"{x:.6g}" for x in row]
        table.add_row(formatted_row)

    print(table)


# In[40]:


plt.figure()
plt.imshow(dd, cmap='viridis')  # You can change the colormap if desired
plt.title('2D Array Plot')
plt.colorbar()  # Adds a colorbar to show the scale
plt.clim([0.9,1.1])
plt.show()


# ## Clean up by deleting unwanted files

# In[41]:


def delete_runtime_files(directory='.'):
    """
    Deletes all files with a .xml extension, files named statepoint.NNNNN.h5 where N is a digit,
    and the file named summary.h5 in the specified directory.

    Parameters:
    directory (str): The directory to search for files. Defaults to the current directory.
    """

    # Regular expression pattern to match files named statepoint.NNNNN.h5
    pattern = re.compile(r'statepoint\.\d{5}\.h5')

    # Iterate through files in the directory
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)

        # Check if the file has a .xml extension, matches the pattern, or is named summary.h5
        if filename.endswith('.xml') or pattern.match(filename) or filename == 'summary.h5' \
           or filename == 'tallies.out':
            # Delete the file
            try:
                os.remove(file_path)
                print(f"Deleted: {file_path}")
            except Exception as e:
                print(f"Error deleting {file_path}: {e}")


# Example usage:
# delete_runtime_files('/path/to/directory')  # specify the directory path if needed


# In[42]:


delete_runtime_files('./')
delete_runtime_files(my_path)


# In[ ]:




