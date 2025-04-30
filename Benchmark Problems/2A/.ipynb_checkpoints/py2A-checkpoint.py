import os
import sys
import numpy as np

sys.path.append("../../..")

from pyopensn.mesh import FromFileMeshGenerator, PETScGraphPartitioner
from pyopensn.xs import MultiGroupXS
from pyopensn.aquad import GLCProductQuadrature2DXY
from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
from pyopensn.context import UseColor, Finalize

UseColor(False)

casename = '2A'
h5_name = '2a'

mesh_filepath = 'lattice_'+casename+'.obj'
meshgen = FromFileMeshGenerator(
    filename=mesh_filepath,
    partitioner=PETScGraphPartitioner(type='parmetis')
)

grid = meshgen.Execute()
grid.ExportToPVTU('mesh_2A')

xs_filepath = 'mgxs_casl_'+h5_name+'/mgxs_'+h5_name+'_one_eighth_SHEM-361.h5'
xs_dict = {}
xs_list = []

h5_mat_names = ['fuel', 'fuel clad', 'fuel gap', 
                'gt-clad', 'gt-water-in', 'gt-water-out', 
                'it-clad', 'it-water-in', 'it-water-out', 
                'moderator', 'water_outside']

for name in h5_mat_names:
    xs_dict[name] = MultiGroupXS()
    xs_dict[name].LoadFromOpenMC(xs_filepath, name, 294.0)
    xs_list = np.append(xs_list,xs_dict[name])

block_ids = [i for i in range(0,len(xs_list))]

scat_order = 3 #xs_list[0].scattering_order

pquad = GLCProductQuadrature2DXY(32, 4)

num_groups = 361
#num_groups = 361

group_sets = [{
            "groups_from_to": (0, num_groups - 1),
            "angular_quadrature": pquad,
            "angle_aggregation_num_subsets": 1,
            "inner_linear_method": "petsc_gmres",
            "l_abs_tol": 1.0e-6,
            "l_max_its": 300,
            #"gmres_restart_interval": 30
            }]

group_sets = [{
            "groups_from_to": (0, num_groups - 1),
            "angular_quadrature": pquad,
            "angle_aggregation_num_subsets": 1,
            "inner_linear_method": "classic_richardson",
            "l_abs_tol": 1.0e-1,
            "l_max_its": 30,
            #"gmres_restart_interval": 30
            }]

# fix this when automating it stops breaking
bound_conditions = [
                        { 'name' : "xmin", 'type' : "reflecting" },
                        { 'name' : "xmax", 'type' : "reflecting" },
                        { 'name' : "ymin", 'type' : "reflecting" },
                        { 'name' : "ymax", 'type' : "reflecting" },
                        { 'name' : "zmin", 'type' : "reflecting" },
                        { 'name' : "zmax", 'type' : "reflecting" }
                        ]
# fix this when automating it stops breaking
xs_mapping = [
            {'block_ids' : [0],'xs' : xs_list[0]},
            {'block_ids' : [1],'xs' : xs_list[1]},
            {'block_ids' : [2],'xs' : xs_list[2]},
            {'block_ids' : [3],'xs' : xs_list[3]},
            {'block_ids' : [4],'xs' : xs_list[4]},
            {'block_ids' : [5],'xs' : xs_list[5]},
            {'block_ids' : [6],'xs' : xs_list[6]},
            {'block_ids' : [7],'xs' : xs_list[7]},
            {'block_ids' : [8],'xs' : xs_list[8]},
            {'block_ids' : [9],'xs' : xs_list[9]},
            {'block_ids' : [10],'xs' : xs_list[10]}
            ]
#xs_mapping = [
#            {'block_ids' : [0],'xs' : xs_list[0]},
#            {'block_ids' : [1],'xs' : xs_list[1]},
#            {'block_ids' : [2],'xs' : xs_list[2]},
#            {'block_ids' : [3],'xs' : xs_list[3]}
#            ]
#xs_mapping = [{'block_ids' : [0],'xs' : xs_list[0]}]

phys = DiscreteOrdinatesProblem(mesh=grid,
                                num_groups=num_groups,
                                groupsets= group_sets,
                                xs_map=xs_mapping
                                        )
phys.SetOptions(scattering_order=scat_order,
                verbose_inner_iterations=True,
                verbose_outer_iterations=True,
                power_default_kappa=1.0,
                power_normalization=1.0,
                save_angular_flux=False,
                #write_restart_time_interval = 3660,
                #write_restart_path = "2A_restart/2A",
                boundary_conditions = bound_conditions
                                    )

k_solver = PowerIterationKEigenSolver(lbs_problem = phys,
                                      k_tol = 1.0e-8
                                     )
k_solver = PowerIterationKEigenSolver(lbs_problem = phys,
                                      k_tol = 1.0e-1
                                     )
k_solver.Initialize()

k_solver.Execute()