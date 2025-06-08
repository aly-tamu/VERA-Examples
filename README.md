# VERA-Examples
VERA CASL benchmark problem 2 modeled with OpenSn and OpenMC for validation and comparison.\
Performed as a final project for NUEN 689 under Dr. Jean Ragusa, based on Pablo Garcia's work.\
Pablo Garcia's repository: https://github.com/pablogarcia44/Garcia_OpenSn.git


## REPOSITORY CONTENTS:
### Benchmark Problems:
Folder with OpenMC and OpenSn work for the CASL benchmark, along with some utility files.
<details>
<summary>Benchmark Problems/ Folder Contents:</summary>
  <hr>
  <ul>
  <li>2[case]:</li>
    <dd>Folder with the OpenMC and OpenSn inputs and outputs for each case,<br/>
    as well as the cross sections from the OpenMC run and the lattice csv file.</dd><br/>
    
  <details>
  <summary>2[case]/ Folder Contents:</summary>
  <hr>
  <ul>
    <li>mgxs_casl_2[case]:</li>
        <dd>Folder containing OpenMC outputs: the multi-group cross section .h5 file and power distribution .npy file.</dd>
        <br/>
    <li>2[case]_v3.ipynb:</li>
        <dd>Jupyter notebook of the OpenMC script for the given case, produces.</dd>
        <br/>
    <li>2[case]_v3.py:</li>
        <dd>Python script version of the 2[case]_v3.ipynb file.<br/>
        Does not remove plotter functions so it may freeze and not finish running.<br/>
        Use the python scripts in the 'OpenMC Python Files' folder instead of these.</dd>
        <br/>
    <li>FA_cell_names_1_family.csv:</li>
        <dd>CSV file of the lattice structure using abbreviations for each input, such as 'fu', 'it' or 'gt'.<br/>
        This file is read by the OpenMC files, OpenSn files, and spydermesh driver when they need the case lattice composition.</dd>
        <br/>
    <li>lattice_2[case].obj:</li>
        <dd>Mesh object output by the spydermesh drivers, which OpenSn uses to run.</dd>
        <br/>
    <li>openmc_2[case]_keff.txt:</li>
        <dd>keff output from OpenMC python file run.</dd>
        <br/>
    <li>opensn-2[case]-keff.py:</li>
        <dd>keff output from OpenSn python file run.</dd>
  </ul>
  <hr>
  </details>
  
<li>make_all_mesh_objs.ipynb:</li>
  <dd>Jupyter notebook to create all OpenSn mesh objects for every case.</dd>
<br/>
<li>OpenMC_h5_reader.py:</li>
  <dd>Python script to read OpenMC cross sections from an hdf5 and plot them.</dd>
<br/>
<li>Opensn_xs_estimate.ipynb:</li>
  <dd>Plots cross sections from OpenSn after it has imported OpenMC MGXS.<br/>
  Was used as a test to understand negative cross section warning when running OpenSn.</dd>
<br/>
<li>power_plotter.ipynb:</li>
  <dd>Python script to take a power.txt output file from an OpenSn benchmark run and plot it.</dd>
<br/>
<li>opensn_sample_input.py:</li>
  <dd>Sample input for running OpenSn from the command line.</dd>
<br/>
<li>spydermesh.py:</li>
  <dd>Mesh object creater designed by Dr. Ragusa.</dd>
<br/>
<li>spydermesh_driver.ipynb:</li>
  <dd>Jupyter notebook to create one case's mesh object at a time based on Pablo Garcia's 'spydermesh_driver_gap.ipynb.'<br/>
  Includes option to cut mesh for symmetry.</dd>
<br/>
<li>spyermesh_driver_gap.ipynb:</li>
  <dd>Jupyter notebook to create a single case's mesh object written by Pablo Garcia.</dd>
  
</ul>
<hr>
</details>

### LUA Problems:               
Folder with Pablo Garcia's LUA OpenSn inputs.\
These were converted into the opensn python inputs in 'Benchmark Problems.'

### OpenMC Python Files:        
OpenMC benchmark inputs and csv files as .py files to easily run from the command line.

### NUEN 689 Presentation.pdf:  
Final presentation for NUEN 689.
