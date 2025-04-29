import numpy as np
import h5py
import math

my_filename = "lattice_2A.obj"
meshgen1 = mesh.MeshGenerator.Create({
    "inputs": [
        mesh.FromFileMeshGenerator.Create({"filename": my_filename})
    ]
})
mesh.MeshGenerator.Execute(meshgen1)
mesh.ExportToPVTU("mesh_2B")

mat_names = [
    'clad_guide', 'clad_instru', 'clad_pincell', 'fuel_pincell', 'gap_pincell',
    'moderator_guide', 'moderator_instru', 'moderator_pincell',
    'water_guide', 'water_instru', 'water_outside'
]
materials = {}
my_xs = {}
xs_file = 'mgxs_2A.h5'

for name in mat_names:
    my_xs[name] = xs.Create()
    xs.Set(my_xs[name], OPENMC_XSLIB, xs_file, 294.0, name)
    mat_handle = mat.AddMaterial(name)
    mat.SetProperty(mat_handle, TRANSPORT_XSECTIONS, OPENMC_XSLIB, xs_file, 294, name)
    materials[name] = mat_handle


pquad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 32, 4)
aquad.OptimizeForPolarSymmetry(pquad, 4.0 * math.pi)

num_groups = 361

lbs_block = {
    "num_groups": num_groups,
    "groupsets": [
        {
            "groups_from_to": [0, num_groups - 1],
            "angular_quadrature_handle": pquad,
            "inner_linear_method": "krylov_richardson",
            "l_max_its": 20,
            "l_abs_tol": 1e-6,
            "angle_aggregation_type": "polar",
        },
    ]
}

lbs_options = {
    "boundary_conditions": [
        {"name": bc, "type": "reflecting"} for bc in ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]
    ],
    "scattering_order": 3,
    "verbose_inner_iterations": True,
    "verbose_outer_iterations": True,
    "power_default_kappa": 1.0,
    "power_normalization": 1.0,
    "save_angular_flux": False,
    "write_restart_time_interval": 3660,
    "write_restart_path": "2A_restart/2A"
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys1, lbs_options)

k_solver = lbs.PowerIterationKEigen.Create({
    "lbs_solver_handle": phys1,
    "k_tol": 1.0e-8,
})

solver.Initialize(k_solver)
solver.Execute(k_solver)

fflist, count = lbs.GetScalarFieldFunctionList(phys1)

if location_id == 0:
    print('Length of the field function list:', len(fflist), count)

MPIBarrier()
log.Log(LOG_0, f'Length of the field function list (v2): {len(fflist)}')


pitch = 1.26
num_cells = 17

def compute_cell_center(i, j):
    offset = (num_cells // 2) * pitch
    x_center = (i - 1) * pitch - offset
    y_center = (j - 1) * pitch - offset
    return x_center, y_center

avoid_pairs = {
    4: [14, 4],
    12: [3, 6, 9, 15, 12],
    3: [6, 12, 9],
    9: [3, 6, 12, 15, 9],
    14: [14, 4],
    15: [6, 12, 9],
    6: [3, 6, 12, 9, 15]
}

def should_avoid(i, j):
    return j in avoid_pairs.get(i, [])

val_table = np.full((num_cells, num_cells), -1.0)

xs_fuel_pincell = xs.Get(my_xs["fuel_pincell"])
sig_f = xs_fuel_pincell.get("sigma_f", None)

if sig_f is None or not isinstance(sig_f, list):
    raise ValueError("sigma_f is not a valid list or missing")

for i in range(1, num_cells + 1):
    for j in range(1, num_cells + 1):
        if should_avoid(i, j):
            continue

        x_center, y_center = compute_cell_center(i, j)
        my_lv = logvol.RCCLogicalVolume.Create({"r": 0.4060, "x0": x_center, "y0": y_center, "z0": -1.0, "vz": 2.0})

        val = 0.0
        for g in range(1, num_groups + 1):
            ffi = fieldfunc.FFInterpolationCreate(VOLUME)
            fieldfunc.SetProperty(ffi, OPERATION, OP_SUM)
            fieldfunc.SetProperty(ffi, LOGICAL_VOLUME, my_lv)
            fieldfunc.SetProperty(ffi, ADD_FIELDFUNCTION, fflist[g - 1])
            fieldfunc.Initialize(ffi)
            fieldfunc.Execute(ffi)
            val_g = fieldfunc.GetValue(ffi)
            val += val_g * sig_f[g - 1]

        val_table[i - 1][j - 1] = val


if location_id == 0:
    with open('power_2A.txt', 'w') as f:
        for i in range(num_cells):
            for j in range(num_cells):
                f.write(f"val_table[{i+1}][{j+1}] = {val_table[i][j]:.5f}\n")

MPIBarrier()
log.Log(LOG_0, 'Finished computing power distribution')

fieldfunc.ExportToVTKMulti(fflist, 'flux_2A')