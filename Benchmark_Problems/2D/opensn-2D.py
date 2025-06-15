import os
import sys
import numpy as np
import csv
from collections import Counter
import matplotlib.pyplot as plt

path = os.getcwd()

sys.path.append("../../..")

casename = "2D"
h5_name = "2d"

if (casename not in path) and ("Benchmark_Problems/" not in path):
    path = path + "/" + casename

mesh_filepath = path + "/" + "lattice_" + casename + ".obj"
meshgen = FromFileMeshGenerator(
    filename=mesh_filepath, partitioner=PETScGraphPartitioner(type="parmetis")
)

grid = meshgen.Execute()
grid.ExportToPVTU("mesh_"+casename)

xs_filepath = path + "/" + "mgxs_casl_" + h5_name + "/mgxs_" + h5_name + "_one_eighth_SHEM-361.h5"
xs_dict = {}
xs_list = []

h5_mat_names = [
    "fuel",
    "clad",
    "gap",
    "gt-clad",
    "gt-water-in",
    "gt-water-out",
    "it-clad",
    "it-water-in",
    "it-water-out",
    "moderator",
    "water_outside",
]

for name in h5_mat_names:
    xs_dict[name] = MultiGroupXS()
    xs_dict[name].LoadFromOpenMC(xs_filepath, name, 294.0)
    xs_list = np.append(xs_list, xs_dict[name])

block_ids = [i for i in range(0, len(xs_list))]

scat_order = 3  # xs_list[0].scattering_order

pquad = GLCProductQuadrature2DXY(4, 32)

num_groups = 361

group_sets = [
    {
        "groups_from_to": (0, num_groups - 1),
        "angular_quadrature": pquad,
        "angle_aggregation_type": "polar",
        "angle_aggregation_num_subsets": 1,
        "inner_linear_method": "classic_richardson",
        "l_abs_tol": 1.0e-4,
        "l_max_its": 300,
    }
]

# fix this when automating it stops breaking
bound_conditions = [
    {"name": "xmin", "type": "reflecting"},
    {"name": "xmax", "type": "reflecting"},
    {"name": "ymin", "type": "reflecting"},
    {"name": "ymax", "type": "reflecting"},
    {"name": "zmin", "type": "reflecting"},
    {"name": "zmax", "type": "reflecting"},
]

# fix this when automating it stops breaking
xs_mapping = [
    {"block_ids": [0], "xs": xs_list[0]},
    {"block_ids": [1], "xs": xs_list[1]},
    {"block_ids": [2], "xs": xs_list[2]},
    {"block_ids": [3], "xs": xs_list[3]},
    {"block_ids": [4], "xs": xs_list[4]},
    {"block_ids": [5], "xs": xs_list[5]},
    {"block_ids": [6], "xs": xs_list[6]},
    {"block_ids": [7], "xs": xs_list[7]},
    {"block_ids": [8], "xs": xs_list[8]},
    {"block_ids": [9], "xs": xs_list[9]},
    {"block_ids": [10], "xs": xs_list[10]},
]

phys = DiscreteOrdinatesProblem(
    mesh=grid, num_groups=num_groups, groupsets=group_sets, xs_map=xs_mapping
)
phys.SetOptions(
    scattering_order=scat_order,
    verbose_inner_iterations=True,
    verbose_outer_iterations=True,
    use_precursors=False,
    power_default_kappa=1.0,
    power_normalization=1.0,
    save_angular_flux=False,
    boundary_conditions=bound_conditions,
    restart_writes_enabled=True,
    write_delayed_psi_to_restart=True,
    write_restart_path="./2D_",
    #read_restart_path="./restart_32_4_tight/2D_",
)


k_solver = PowerIterationKEigenSolver(lbs_problem=phys, k_tol=1.0e-6)
k_solver.Initialize()
k_solver.Execute()

keff = k_solver.GetEigenvalue()
fflist = phys.GetScalarFieldFunctionList()


def compute_cell_center_old(i, j, num_cells, pitch):
    """
    i, j are in {0, 1, …, num_cells-1}.
    Returns (x_center, y_center) so that:
      - When num_cells is odd, the cell ( (num_cells-1)//2, (num_cells-1)//2 ) sits at (0,0).
      - When num_cells is even, the grid is centered between the four middle cells.
    """
    center_index = (num_cells - 1) / 2.0
    x_center = (i - center_index) * pitch
    y_center = (j - center_index) * pitch
    return x_center, y_center


def compute_cell_center(i, j, offset_x, offset_y):
    """
    i, j are in {0, 1, …, num_cells-1}.
    Returns (x_center, y_center)
    """
    x_center = i * pitch + offset_x 
    y_center = j * pitch + offset_y
    return x_center, y_center

def read_csv_to_2d_array(file_path):
    with open(file_path, newline="", encoding="utf-8") as csvfile:
        reader = csv.reader(csvfile)
        data = [row for row in reader]
    return np.asarray(data)


def count_frequencies(data):
    flattened_data = [item for row in data for item in row]  # Flatten 2D array into a 1D list
    return Counter(flattened_data)


csv_filename = "FA_cell_names_1_family.csv"
csv_filepath = path + "/" + csv_filename
lattice_csv = read_csv_to_2d_array(csv_filepath)

cell_frequencies = count_frequencies(lattice_csv)
if rank == 0:
    print("cell name frequency:")
    total = 0
    for key, value in cell_frequencies.items():
        print(f'"{key}": {value}')
        total += value
    print("total: ", total)

num_cells = lattice_csv.shape[0]
if num_cells != lattice_csv.shape[1]:
    raise Exception("CSV array of cell names is not square.")

fuel_xs = xs_dict["fuel"]
sig_f = np.array(fuel_xs.sigma_f)

pitch = 1.26

num_cells_quarter = np.ceil(num_cells/2).astype(np.int64)
quarter_lattice_center = -5.375

offset_x = -4.705
offset_y =  4.705 - 16*pitch

val_table = np.zeros([num_cells, num_cells])

for i in range(num_cells):
    for j in range(num_cells):
        if lattice_csv[i, j] == "fu":
            x_center, y_center = compute_cell_center(i, j, offset_x, offset_y)
            if rank == 0:
                print("centers=", i, j, x_center, y_center)
            my_lv = RCCLogicalVolume(r=0.4060, x0=x_center, y0=y_center, z0=-1.0, vz=2.0)

            val = 0
            for g in range(0, num_groups):
                ffi = FieldFunctionInterpolationVolume()
                ffi.SetOperationType("sum")
                ffi.SetLogicalVolume(my_lv)
                ffi.AddFieldFunction(fflist[g])
                ffi.Initialize()
                ffi.Execute()
                val_g = ffi.GetValue()
                val += val_g * sig_f[g]
            val_table[i, j] = val

val_table_ori = val_table.copy()

val_table = np.flip(val_table, axis=1)
# quarter array
A = val_table[:num_cells_quarter, :num_cells_quarter]
A[:,-1] *= 2.
A[-1,:] *= 2.
val_table_quarter = A.copy()
A_flipped = np.flip(A, axis=1)
B = np.hstack([A,A_flipped[:,1:]])
B_flipped = np.flip(B, axis=0)
val_table = np.vstack([B,B_flipped[1:,:]])

norm = np.sum(val_table) / cell_frequencies["fu"]
val_table /= norm

MPIBarrier()

if rank == 0:
    #print("sum=", np.sum(val_table))
    np.savetxt("power.txt", val_table)
    #np.savetxt("power_quarter.txt", val_table_quarter)
    #np.savetxt("power_ori.txt", val_table_ori)
    with open("keff.txt", "w") as file:
        file.write(str(keff))

vtk_basename = "flx_4_32_"
#FieldFunctionGridBased.ExportMultipleToVTK([fflist[g] for g in range(num_groups)], vtk_basename)
