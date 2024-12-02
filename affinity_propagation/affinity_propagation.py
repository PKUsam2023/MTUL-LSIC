import numpy as np
import os
import argparse
from sklearn.cluster import AffinityPropagation
import pandas as pd
# import time
import sys


__author__ = "Dong Chen"
__date__ = "Jan. 18, 2021"


def data_norm(*args):
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    scaler.fit(args[0])
    norm_args = [scaler.transform(args[i]) for i in range(len(args))]
    norm_args = norm_args if len(args) > 1 else norm_args[0]
    return norm_args


def get_labels(args, filenames, max_iteration=1000):
    rand_seed = 1234
    np.random.seed(rand_seed)

    features_data = pd.read_csv(filenames[0], header=0, index_col='cif_id')
    features = features_data.values
    # features = data_norm(features)

    clustering = AffinityPropagation(
        random_state=rand_seed,
        max_iter=max_iteration
    ).fit(features)

    cluster_label_df = pd.DataFrame(
        clustering.labels_,
        index=list(features_data.index),
        columns=['label']
    )
    print(f"number of cluster: {len(clustering.cluster_centers_indices_)}")
    cluster_label_df.to_csv(
        os.path.join(args.save_dirname, f"affinity_clusters_{len(clustering.cluster_centers_indices_)}_Li_free__.csv")
    )


def parse_args(args):
    parser = argparse.ArgumentParser(description='run k-means')
    parser.add_argument('--filenames', nargs='+', default=None, type=str)
    parser.add_argument('--max_clusters', default=10, type=int)
    parser.add_argument('--max_iteration', default=20, type=int)
    parser.add_argument('--save_dirname', default='./', type=str)
    args = parser.parse_args()
    return args


def main(args):
    get_labels(args, args.filenames, args.max_iteration)


def cli_main():
    args = parse_args(sys.argv[1:])
    print(args)
    main(args)


if __name__ == "__main__":
    cli_main()
    print('End!')
