'''
barcode to topological feature
'''


import numpy as np
import sys
import os
import argparse
import pandas as pd


__author__ = "Dong Chen"
__date__ = "Oct. 22, 2020"


class generate_feature(object):
    """
    Introduction: split the barcode in some way
    """

    def __init__(self, bins_range=0.1, max_dis=10):
        self.bins_range = bins_range
        self.max_dis = max_dis

    # --------------------------split barcode------------------------
    def split_vb(self, betti_num, bins_range, max_dis):
        # the shape of betti_num must be (n,2)
        bins_num = int(max_dis/bins_range)
        birth_vector = np.zeros((1, bins_num), dtype=int)
        for j in range(0, bins_num):
            birth_vector[0][j] = np.shape(
                betti_num[
                    (betti_num[:, 0] >= j*bins_range) * (betti_num[:, 0] <= (j+1)*bins_range)
                ]
            )[0]

        return birth_vector

    def split_vd(self, betti_num, bins_range, max_dis):
        # the shape of betti_num must be (n,2)
        bins_num = int(max_dis/bins_range)
        death_vector = np.zeros((1, bins_num), dtype=int)
        for j in range(0, bins_num):
            death_vector[0][j] = np.shape(
                betti_num[
                    (betti_num[:, 1] >= j*bins_range) * (betti_num[:, 1] <= (j+1)*bins_range)
                ]
            )[0]
        return death_vector

    def split_vp(self, betti_num, bins_range, max_dis):
        # the shape of betti_num must be (n,2)
        bins_num = int(max_dis/bins_range)
        bt_persist = betti_num
        persist_vector = np.zeros((1, bins_num), dtype=int)
        for j in range(0, bins_num):
            persist_vector[0][j] = np.shape(
                bt_persist[
                    (bt_persist[:, 0] <= j*bins_range) * (bt_persist[:, 1] >= (j+1)*bins_range)
                ]
            )[0]
        return persist_vector

    def statistic_vp(self, betti_num):
        # only statistic the info of betti_num 0
        # just consider the death value between 0.1 and 10, filter the toplogy noise
        betti_value_consider = list(filter(lambda v: v > 0.1 and v <= 10, betti_num[:, 1]))
        
        # reflect the atom density
        atom_n = len(betti_value_consider)
        
        if atom_n > 0:
            betti_feature = [
                len(betti_value_consider),  # reflect the atom density
                np.min(betti_value_consider),
                np.max(betti_value_consider),
                np.mean(betti_value_consider),
                np.std(betti_value_consider),
                np.sum(betti_value_consider),
                np.median(betti_value_consider),
            ]
        else:
            betti_feature = [0] * 7
        return betti_feature

    def statistic_value(self, data):
        if np.max(data) > self.max_dis:
            for i, v in enumerate(data):
                data[i] = v if v <= self.max_dis else self.max_dis
        return [np.min(data), np.max(data), np.mean(data), np.std(data), np.sum(data)]

    def statistic_feature(self, betti_num):
        # statistic features from betti1 and betti2
        # min_birth, max_birth, mean_birth, std_birth, sum_birth
        # min_death, max_death, mean_death, std_death, sum_death
        # min_length, max_length, mean_length, std_length, sum_length

        statistic_feature_num = 5*3
        if len(betti_num) == 0:
            betti_feature = [0] * statistic_feature_num  # total features number
        else:
            birth_statistic = self.statistic_value(betti_num[:, 0])
            death_statistic = self.statistic_value(betti_num[:, 1])
            length_statistic = self.statistic_value(betti_num[:, 1]-betti_num[:, 0])
            betti_feature = birth_statistic + death_statistic + length_statistic
        return betti_feature

    def bar2feature(self, dgms_data):
        # betti0
        persist_vector = self.statistic_vp(dgms_data[0])
        # betti-1
        bettin_feature = np.array([])
        for i in range(1, 2):
            bettin_feature = np.append(bettin_feature, self.statistic_feature(dgms_data[i]))
        out_features = np.append(persist_vector, bettin_feature)  # shape=(max_dis/bins_range+15*1, )
        return out_features

    def bar2feature_split_all(self, dgms_data):
        # betti0
        persist_vector = self.split_vp(dgms_data[0], self.bins_range, self.max_dis)
        # betti-1
        persist_vector_1 = self.split_vp(dgms_data[1], self.bins_range, self.max_dis)
        out_features = np.append(persist_vector, persist_vector_1)  # shape=(2*max_dis/bins_range, )
        return out_features


def from_rips_alpha_complex(filenames, file_idx, args):
    feature_generator = generate_feature(bins_range=args.bins_range, max_dis=args.cutoff_dis)
    if len(filenames) > 1:
        Li_out_features = np.zeros([len(filenames), 22])
        other_out_features = np.zeros([len(filenames), 22])
        for i, filename in enumerate(filenames):
            print(i, filename)
            barcode_data = np.load(filename, allow_pickle=True).item()
            Li_out_features_temp = []
            other_out_features_temp = []
            for key, value in barcode_data.items():
                Li_out_features_temp.append(feature_generator.bar2feature(value[0]))
                other_out_features_temp.append(feature_generator.bar2feature(value[1]))
            Li_out_features[i, :] = np.mean(Li_out_features_temp, axis=0)
            other_out_features[i, :] = np.mean(other_out_features_temp, axis=0)
        Li_out_features_df = pd.DataFrame(Li_out_features, index=file_idx)
        other_out_features_df = pd.DataFrame(other_out_features, index=file_idx)
        Li_out_features_df.to_csv(f"Li_tp_from_{os.path.split(args.dirname)[-1]}.csv")
        other_out_features_df.to_csv(f"other_tp_from_{os.path.split(args.dirname)[-1]}.csv")
    else:
        barcode_data = np.load(filenames[0], allow_pickle=True).item()
        Li_out_features = []
        other_out_features = []
        for key, value in barcode_data.items():
            print(key, np.shape(value[0][0]), np.shape(value[0][1]), np.shape(value[0][1]))
            Li_out_features.append(feature_generator.bar2feature(value[0]))
            other_out_features.append(feature_generator.bar2feature(value[1]))
        Li_out_features_df = pd.DataFrame(Li_out_features, index=file_idx)
        other_out_features_df = pd.DataFrame(other_out_features, index=file_idx)
        Li_out_features_df.to_csv(f"Li_tp_from_{os.path.split(args.dirname)[-1]}.csv")
        other_out_features_df.to_csv(f"other_tp_from_{os.path.split(args.dirname)[-1]}.csv")


def parse_args(args):
    parser = argparse.ArgumentParser(description='get topological features')
    parser.add_argument('--file_type', default='Li_other', type=str)
    parser.add_argument('--single_file', default=False, action='store_true')
    parser.add_argument('--filename', default=None, type=str)
    parser.add_argument('--dirname', default=None, type=str)
    parser.add_argument('--save_dirname', default='./', type=str)
    parser.add_argument('--bins_range', default=0.1, type=float,
                        help='The length of each bin')
    parser.add_argument('--cutoff_dis', default=10, type=float,
                        help='Max distance to consider during splitting the barcode.')
    parser.add_argument('--file_ids', default=None, type=str,
                        help='record the information of all files, format=.csv')
    parser.add_argument('--all_split', default=False, action='store_true')
    args = parser.parse_args()
    return args


def main(args):
    if args.single_file:
        print(args.bins_range, args.cutoff_dis)
        if args.file_type == 'Li_other':
            from_rips_alpha_complex([args.filename], args)

    else:
        file_idx = pd.read_csv(args.file_ids, header=0, index_col=0)['cif_id']
        all_files = []
        for i, id_n in enumerate(file_idx):
            all_files.append(os.path.join(args.dirname, f"{id_n}_POSCAR.npy"))
        if args.file_type == 'Li_other':  # include anion
            from_rips_alpha_complex(all_files, file_idx, args)


def cli_main():
    args = parse_args(sys.argv[1:])
    print(args)
    main(args)


if __name__ == "__main__":
    cli_main()
    print('End!')
