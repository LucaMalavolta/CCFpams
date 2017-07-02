import numpy as np
import os
import argparse

parser = argparse.ArgumentParser(prog='rdb2input.py', description='read a .rdb file to create the input file required by CCFPams')
parser.add_argument('rdb_directory', type=str, nargs=1, help='config file')
parser.add_argument('star_name', type=str, nargs=1, help='config file')


args = parser.parse_args()
rdb_directory = args.rdb_directory[0]
star_name = args.star_name[0]


fileout_input = open(star_name+'_input.list','w')

data = np.genfromtxt(rdb_directory + '/' +  star_name+'_harpn_extract.rdb', delimiter='\t', dtype=str, skip_header=2)
for d_file, d_data, d_mask, d_simu in zip(data[:,0], data[:,1], data[:,52], data[:,72]):
    fileout_input.write(d_data + '   ' + d_file + '   ' + d_mask + '   A   ' + star_name + '  \n')

print  rdb_directory,  star_name, star_name+'_input.list succesfully created'

fileout_input.close()
