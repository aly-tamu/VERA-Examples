# VERA-Examples
VERA CASL benchmark problem 2 modeled with OpenSn and OpenMC for validation and comparison.
Performed as a final project for NUEN 689 under Dr. Jean Ragusa, based on Pablo Garcia's work.
Pablo Garcia's repository: https://github.com/pablogarcia44/Garcia_OpenSn.git

REPOSITORY CONTENTS:
Benchmark Problems :         Folder with OpenMC and OpenSn work for the CASL benchmark,
                             along with some utility files.
LUA Problems :               Folder with Pablo Garcia's LUA OpenSn inputs. These were 
                             converted into the opensn python inputs in 'Benchmark Problems.'
OpenMC Python Files :        OpenMC benchmark inputs and csv files as .py files to easily
                             run from the command line.
NUEN 689 Presentation.pdf :  Final presentation NUEN 689.
                      
'Benchmark Problems/' Folder Contents:
2[case] :                    Folder with the OpenMC and OpenSn inputs and outputs for each case,
                             as well as the cross sections from the OpenMC run and the lattice 
                             csv file.
make_all_mesh_objs.ipynb :   Jupyter notebook to create all OpenSn mesh objects for every case.
OpenMC_h5_reader.py :        Python script to read OpenMC cross sections from an hdf5 and plot them.
Opensn_xs_estimate.ipynb :   Plots cross sections from OpenSn after it has imported OpenMC MGXS.
                             Was used as a test to understand negative cross section warning when 
                             running OpenSn.
power_plotter.ipynb :        Python script to take a power.txt output file from an OpenSn benchmark
                             run and plot it.
opensn_sample_input.py :     Sample input for running OpenSn from the command line.
spydermesh.py :              Mesh object creater designed by Dr. Ragusa.
spydermesh_driver.ipynb :    Jupyter notebook to create one case's mesh object at a time based on Pablo
                             Garcia's 'spydermesh_driver_gap.ipynb.' Includes option to cut mesh for symmetry.
spyermesh_driver_gap.ipynb : Jupyter notebook to create a single case's mesh object written by Pablo Garcia.

'2[case]/' Folder Contents:
mgxs_casl_2[case] :          Folder containing OpenMC outputs: the multi-group cross section .h5 file and power 
                             distribution .npy file.
2[case]_v3.ipynb :           Jupyter notebook of the OpenMC script for the given case, produces.   
2[case]_v3.py :              Python script version of the 2[case]_v3.ipynb file. Does not remove plotter
                             functions so it may freeze and not finish running. Use the python scripts 
                             in the 'OpenMC Python Files' folder instead of these.
FA_cell_names_1_family.csv : CSV file of the lattice structure using abbreviations for each input, such as 'fu',
                             'it' or 'gt'. This file is read by the OpenMC files, OpenSn files, and spydermesh 
                             driver when they need the case lattice composition.
lattice_2[case].obj :        Mesh object output by the spydermesh drivers, which OpenSn uses to run.
openmc_2[case]_keff.txt :    keff output from OpenMC python file run.
opensn-2[case]-keff.py :     keff output from OpenSn python file run.    