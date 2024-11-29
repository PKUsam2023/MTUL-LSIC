import pandas as pd
import numpy as np


# label_file = r"affinity_clusters_34.csv"  # move same stru + Li only filter + Li free filter + anion tp features similarity(no norm)
# target_label = [9, 16, 18, 20, 21, 26, 27, 31]
label_file = r"affinity_clusters_32_Li_free.csv"  # move same stru + Li only filter unit cell + Li free filter unit cell + anion tp features similarity (no norm)
# target_label = [4, 5, 7, 9, 10, 15, 17, 43, 52]
target_label = [4, 5, 7, 9, 10, 15, 17, 26]  # label for affinity_clusters_32_Li_free.csv
label_data = pd.read_csv(label_file, header=0, index_col=0)
target_ids = label_data[label_data['label'].isin(target_label)].index
print(len(target_ids))  # 403

known_ids = [174202, 191542, 189530, 185539, 185540, 238693, 238692, 238690, 238689, 195182, 68251, 68252, 69076, 71212, 71213, 202768, 202769, 171168, 171169, 171170, 171171, 174612, 174615, 174617, 54867, 154400, 157628, 158372, 158373, 159426, 159731, 159732, 161342, 161343, 161386, 161387, 161388, 161389, 163860, 163861, 163914, 422259, 182034, 182035, 182036, 182312, 245641, 245934, 245935, 245936, 245937, 419624, 419625, 419626, 246817, 261302, 183684, 183685, 183686, 183687, 183873, 184230, 191524, 189529, 191529, 185602, 191525, 191527, 191526, 191528, 183607, 237202, 238687, 238686, 237145, 195435, 237199, 195437, 195436, 237200, 62629, 62244, 50422, 65025, 69345, 69347, 51333, 92242, 98361, 240269, 151919, 62300, 62301, 50420, 86457, 91113, 78503, 83502, 202539, 89992, 94029, 94030, 94031, 95973, 55751, 55752, 163297, 163299, 163300, 163301, 163302, 163303, 163304, 163305, 166862, 184087, 193036, 290805, 193037, 193032, 193038, 193039, 193033, 62333, 62878, 69763, 69764, 69765, 78504, 83832, 201935, 89456, 92250, 97658, 97659, 97660, 191891, 83831, 190657, 190656, 245527, 257190, 20208, 25816, 34221, 50058, 66576, 66577, 79427, 77095, 415976, 174385, 109092, 150918, 150919, 161792, 167251, 167252, 167253, 8280, 65643, 79537, 78513, 100167, 100169, 202631, 202632, 166547, 426050, 150920, 150921, 238600, 238601, 157654, 180319, 290831, 290832, 35018, 50578, 92200, 95649, 188886, 188887, 193755, 193947, 241439, 252036, 252037, 252038, 252039, 252040, 252041, 95785]

# target_ids = [i for i in target_ids if i not in known_ids]
target_ids = [i for i in target_ids]
print(len(target_ids))

# feature_file = r"../top_simi_stru.csv"
feature_file = r"../top_simi_stru_from_uc_list_same.csv"
feature_data = pd.read_csv(feature_file, header=0, index_col=0)
target_features = feature_data.loc[target_ids]
print(np.shape(target_features))
# target_features.to_csv('target_cluster_stru.csv')
target_features.to_csv('affinity_clusters_result.csv')  # add 26 label
