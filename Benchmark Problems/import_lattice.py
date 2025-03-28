import os
import numpy as np
import csv

os.getcwd()
#os.chdir(os.path.expanduser('~')+'/GitHub/VERA-Examples/Benchmark Problems')
#print(os.getcwd())
#file_path = 'Material Files/3.1_fuel.csv'

def csv_to_matrix(filepath):
    with open(filepath, newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        csv_mat = [row for row in reader]
        csv_mat = np.array(csv_mat).transpose()
        return csv_mat