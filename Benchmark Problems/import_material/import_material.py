import openmc
import os
import numpy as np
import csv

#os.getcwd()
#os.chdir(os.path.expanduser('~')+'/GitHub/VERA-Examples/Benchmark Problems')
#print(os.getcwd())
#file_path = 'Material Files/3.1_fuel.csv'


def csv_to_matrix(filepath):
    with open(filepath, newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        csv_mat = [row for row in reader]
        csv_mat = np.array(csv_mat)[1::,:]
        return csv_mat

def mat_from_csv(filepath):
    csv_mat = csv_to_matrix(filepath)
    mat_name = str(csv_mat[0,0])
    mat = openmc.Material(name=mat_name)
    mat.set_density('g/cm3',float(csv_mat[0,3]))
    for row in csv_mat:
        mat.add_nuclide(str(row[1]),float(row[2]),'ao')
    return mat