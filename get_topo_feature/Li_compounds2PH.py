'''
    convert Li-containing compounds into PH features
    Rips complexes
'''

import numpy as np
import glob
import os
import sys
import argparse
import ripser as rp
import time


__author__ = "Dong Chen"
__date__ = "Oct. 10, 2020"


class Compounds2PH(object):
    def __init__(self, r_sphere=5, supercell_extend_dis=15):
        self.r_sphere = r_sphere
        self.supercell_extend_dis = supercell_extend_dis

    def poscar2newcell(self, filename, supercell_extend_dis=15):
        # read the poscar file
        with open(filename, 'r') as fr:
            lines = fr.readlines()
        
        # read the convert matrix
        convert_matrix = np.empty((3, 3), dtype=float)
        for i in range(0, 3):
            d1 = list(map(float, lines[i+2].strip().split()))
            convert_matrix[i][:] = d1
        convert_matrix = np.round(convert_matrix, 5)

        # read unit cell elements
        ele_type = lines[5].strip().split()
        each_ele_num = list(map(int, lines[6].strip().split()))
        unitcell_elements = np.array([])
        cont = 0
        for j in range(len(ele_type)):
            for k in range(each_ele_num[j]):
                unitcell_elements = np.append(unitcell_elements, ele_type[j])
        
        # read unit cell relative position
        if 's' == lines[7][0] or 'S' == lines[7][0]:
            print(lines[7])
            start_line = 9
        else:
            start_line = 8
        unitcell_relative_coordinates = np.empty((sum(each_ele_num), 3), dtype=float)
        for l in range(sum(each_ele_num)):
            d2 = list(map(float, lines[l+start_line].strip().split()[0:3]))
            unitcell_relative_coordinates[l][:] = d2
        unitcell_relative_coordinates = np.round(unitcell_relative_coordinates, 4)

        # creat new supercell
        # extend layer, 333 at least
        extend_layer = 0
        convert_m_dis = [np.linalg.norm(v1-[0, 0, 0]) for v1 in convert_matrix]
        convert_m_dis_min = min(convert_m_dis)
        while (convert_m_dis_min*extend_layer - supercell_extend_dis) < 0:
            extend_layer += 1
            convert_m_dis_min = convert_m_dis_min * extend_layer
        
        # unit cell elements number
        unitcell_elements_num = unitcell_elements.shape[0]

        # extend, 333 or more
        supercell_relative_coordinates = np.empty((unitcell_elements_num*(1+extend_layer*2)**3, 3), dtype=float)

        # (0 0 0) is first
        displacement = list(range(extend_layer+1)) + list(range(-extend_layer, 0))
        supercell_elements = np.array([], dtype=np.float64)
        start_p = 0
        for i in range(1+extend_layer*2):
            for j in range(1+extend_layer*2):
                for k in range(1+extend_layer*2):
                    dis_p = [displacement[i], displacement[j], displacement[k]]
                    supercell_relative_coordinates[start_p: start_p+unitcell_elements_num][:] = \
                        unitcell_relative_coordinates + dis_p
                    supercell_elements = np.append(
                        supercell_elements, list(unitcell_elements)
                    )
                    start_p += unitcell_elements_num
        supercell_relative_coordinates = np.round(supercell_relative_coordinates, 5)

        # unitcell, supercell coordinates
        unitcell_coordinates = np.round(np.dot(unitcell_relative_coordinates, convert_matrix), 4)
        supercell_coordinates = np.round(np.dot(supercell_relative_coordinates, convert_matrix), 4)

        return unitcell_elements, unitcell_coordinates, supercell_elements, supercell_coordinates

    def get_sub_environment(
            self, unitcell_elements, unitcell_coordinates,
            supercell_elements, supercell_coordinates,
            consider_distance=5
        ):
        """
        Arguments:
            unitcell_elements {array} - Unitcell's elements
            unitcell_coordinates {array} - Unitcell's absolute coordinates
            supercell_elements {array} - Supercell's elements
            supercell_coordinates {array} - Supercell's absolute coordinates
        Keyword Arguments:
            consider_distance {int} - Radius of demon for each element (default: 5)
        Returns:
            sub_environment_index {list} - Each element's neighbors' index
        """

        # number of unitcell elements and supercell elements
        unitcell_ele_num = unitcell_elements.shape[0]
        supercell_elements_num = supercell_elements.shape[0]

        # initial sub_environment
        sub_environment_index = []
        sub_environment_Li_index = []
        center_Li_coordinates = []

        for i in range(unitcell_ele_num):
            center_element = unitcell_elements[i]
            if center_element == 'Li':
                center_coordinates = unitcell_coordinates[i]
                neighbor_index = []
                neighbor_Li_index = []
                for j in range(supercell_elements_num):
                    distance_value = np.round(
                        np.linalg.norm(center_coordinates - supercell_coordinates[j]), 4)

                    # index of elements whose distance within consider distance
                    if distance_value <= consider_distance and distance_value != 0:
                        if supercell_elements[j] == 'Li':
                            neighbor_Li_index.append(j)
                        else:
                            neighbor_index.append(j)

                # all Li's neighbors index
                sub_environment_index.append(neighbor_index)
                sub_environment_Li_index.append(neighbor_Li_index)
                center_Li_coordinates.append(center_coordinates)
            else:
                continue

        return sub_environment_index, sub_environment_Li_index, center_Li_coordinates

    def get_barcode(
        self, unitcell_elements, unitcell_coordinates,
        supercell_elements, supercell_coordinates,
        sub_environment_index, sub_environment_Li_index, center_Li_coordinates
    ):
        PH_dict = {}
        for n, center_coordinates in enumerate(center_Li_coordinates):
            sphere_Li_coords = np.vstack(
                [center_coordinates, supercell_coordinates[sub_environment_Li_index[n]]])
            sphere_other_coords = supercell_coordinates[sub_environment_index[n]]

            # get PH features
            print(np.shape(sphere_Li_coords))
            Li_ph = rp.ripser(sphere_Li_coords, maxdim=2, metric='euclidean', thresh=self.r_sphere)['dgms']
            other_ph = rp.ripser(sphere_other_coords, maxdim=2, metric='euclidean', thresh=self.r_sphere)['dgms']
            PH_dict[f'Li_{n}'] = [Li_ph, other_ph]
        return PH_dict

    def from_poscar(self, filename):
        unitcell_elements, unitcell_coordinates, supercell_elements, supercell_coordinates = \
            self.poscar2newcell(filename, supercell_extend_dis=self.supercell_extend_dis)
        sub_environment_index, sub_environment_Li_index, center_Li_coordinates = \
            self.get_sub_environment(
                unitcell_elements, unitcell_coordinates, supercell_elements,
                supercell_coordinates, consider_distance=self.r_sphere)
        PH_dict = self.get_barcode(
            unitcell_elements, unitcell_coordinates, supercell_elements, supercell_coordinates,
            sub_environment_index, sub_environment_Li_index, center_Li_coordinates)
        return PH_dict

    def single_from_poscar(self, poscar_filename, save_dirname):
        PH_dict = self.from_poscar(poscar_filename)

        save_filename = os.path.join(
            save_dirname,
            os.path.split(poscar_filename)[-1]
        )
        np.save(save_filename, PH_dict)
        print(poscar_filename)
        # try:
        #     PH_dict = self.from_poscar(poscar_filename)

        #     save_filename = os.path.join(
        #         save_dirname,
        #         os.path.split(poscar_filename)[-1]
        #     )
        #     np.save(save_filename, PH_dict)
        #     print(poscar_filename)
        # except:
        #     print('failed', poscar_filename)

    def batch_from_poscar(self, poscar_dirname, save_dirname, issues_log):
        log = open(issues_log, 'w')
        poscar_filenames = glob.glob(os.path.join(poscar_dirname, '*'))
        for i, poscar_file in enumerate(poscar_filenames):
            t0 = time.time()
            try:
                PH_dict = self.from_poscar(poscar_file)

                save_filename = os.path.join(
                    save_dirname,
                    os.path.split(poscar_file)[-1].split('.cif')[0]
                )

                np.save(save_filename, PH_dict)
            except:
                print('failed', i, poscar_file, end='\n', file=log)
            print(i, time.time()-t0, poscar_file, end='\n', file=log)
        log.close()

    def from_cif(self, file):
        import pymatgen as mg
        structure = mg.Structure.from_str(open(file).read(), fmt='cif')

        all_Li_sites = []
        for i, it in enumerate(structure):
            sp = []
            sp_occ = []
            has_Li = False
            for s, v in it.species.items():
                if s.symbol == 'Li':
                    has_Li = True
                    break
                else:
                    sp.append(s.symbol)
                    sp_occ.append(v)
            if has_Li:
                structure[i] = 'Li'
            else:
                structure[i] = sp[np.argmax(sp_occ)]

            if has_Li:  # get Li sites
                all_Li_sites.append(structure.sites[i])

        PH_dict = {}
        for j, Li_site in enumerate(all_Li_sites):
            sphere_Li_coords = np.vstack([Li_site.coords, ])
            sphere_other_coords = np.vstack([Li_site.coords, ])  # first is Li atom

            # get the sphere centered with Li atom
            sphere = structure.get_all_neighbors(r=self.r_sphere, sites=[Li_site])[0]
            for k, site in enumerate(sphere):
                coords = site[0].coords
                symbol = [symbol_type.symbol for symbol_type, _ in site.species.items()]
                if 'Li' in symbol:
                    sphere_Li_coords = np.vstack([sphere_Li_coords, coords])
                else:
                    sphere_other_coords = np.vstack([sphere_other_coords, coords])

            # get PH features
            Li_ph = rp.ripser(sphere_Li_coords, maxdim=2, metric='euclidean', thresh=self.r_sphere)['dgms']
            other_ph = rp.ripser(
                sphere_other_coords[1::, :], maxdim=2, metric='euclidean', thresh=self.r_sphere)['dgms']
            PH_dict[f'Li_{j}'] = [Li_ph, other_ph]
        return PH_dict

    def single_from_cif(self, cif_filename, save_dirname):
        try:
            PH_dict = self.from_cif(cif_filename)

            save_filename = os.path.join(
                save_dirname,
                os.path.split(cif_filename)[-1].split('.cif')[0]
            )
            np.save(save_filename, PH_dict)
            print(cif_filename)
        except:
            print('failed', cif_filename)

    def batch_from_cif(self, cif_dirname, save_dirname, issues_log):
        log = open(issues_log, 'w')
        cif_filenames = glob.glob(os.path.join(cif_dirname, '*.cif'))
        for i, cif_file in enumerate(cif_filenames):
            t0 = time.time()
            try:
                PH_dict = self.from_cif(cif_file)

                save_filename = os.path.join(
                    save_dirname,
                    os.path.split(cif_file)[-1].split('.cif')[0]
                )

                np.save(save_filename, PH_dict)
            except:
                print('failed', i, cif_file, end='\n', file=log)
            print(i, time.time()-t0, cif_file, end='\n', file=log)
        log.close()


def parse_args(args):
    parser = argparse.ArgumentParser(description='convert cif to PH features')

    parser.add_argument('--cif_dirname', default=None, type=str)
    parser.add_argument('--poscar_dirname', default=None, type=str)
    parser.add_argument('--save_dirname', default=None, type=str)
    parser.add_argument('--issues_log', default='convert.log', type=str)
    parser.add_argument('--single_convert', action='store_true', default=False)
    parser.add_argument('--from_poscar', action='store_true', default=False)
    parser.add_argument('--cif_filename', default=None, type=str)
    parser.add_argument('--poscar_filename', default=None, type=str)
    parser.add_argument('--sphere_radius', default=5.0, type=float)
    args = parser.parse_args()
    return args


def main(args):
    PH = Compounds2PH(args.sphere_radius)
    if args.from_poscar:
        if args.single_convert:
            PH.single_from_poscar(args.poscar_filename, args.save_dirname)
        else:
            PH.batch_from_poscar(args.poscar_dirname, args.save_dirname, args.issues_log)
    else:
        if args.single_convert:
            PH.single_from_cif(args.cif_filename, args.save_dirname)
        else:
            PH.batch_from_cif(args.cif_dirname, args.save_dirname, args.issues_log)


def cli_main():
    args = parse_args(sys.argv[1:])
    print(args)
    main(args)


if __name__ == "__main__":
    cli_main()
    print('End!')

# from poscar
# from cif
