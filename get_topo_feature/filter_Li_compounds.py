'''
filter the Li-compounds from ICSD dataset
'''


import pymatgen as mg
import shutil
import argparse
import glob
import sys
import numpy as np
import os


def filter_Li_compounds(args):
    log = open(args.issues_log, 'w')
    filenames = glob.glob(os.path.join(args.dirname, '*.cif'))
    Li_filenames_flag = np.zeros([len(filenames), ])
    for i, file in enumerate(filenames):
        try:
            structure = mg.Structure.from_file(file)
            for j, sp in enumerate(structure.types_of_specie):
                if sp.symbol == 'Li':
                    Li_filenames_flag[i] = 1
                    break
                else:
                    continue
        except:
            print(file, end='\n', file=log)
    log.close()
    Li_filenames = np.array(filenames)[np.where(Li_filenames_flag == 1)]

    save_dirname = args.save_dirname
    for i, file in enumerate(Li_filenames):
        shutil.copy2(file, save_dirname)
    return None


def parse_args(args):
    parser = argparse.ArgumentParser(description='filter Li compounds')
    parser.add_argument('--dirname', default=None, type=str,
                        help='The path of the folder containing .cif files')
    parser.add_argument('--save_dirname', default=None, type=str,
                        help='The path of the folder saving Li-compounds .cif files')
    parser.add_argument('--issues_log', type=str, default='filter.log')
    args = parser.parse_args()
    return args


def main(args):
    filter_Li_compounds(args)


def cli_main():
    args = parse_args(sys.argv[1:])
    print(args)
    main(args)


if __name__ == "__main__":
    cli_main()
