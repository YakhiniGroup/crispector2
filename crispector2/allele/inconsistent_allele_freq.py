from scipy.spatial import distance
import numpy as np
import sys


def analyze_gap_and_JS(all_site_dict):
    problematic_sites = list()
    for site_name, site_info in all_site_dict.items():
        sorted_by_mock_tx = dict()
        sorted_mock = dict(sorted(site_info['mock'].items(), key=lambda item: item[1], reverse=True))
        for allele, freq in site_info['mock'].items():
            sorted_by_mock_tx[allele] = site_info['tx'][allele]
        array_mock = np.array(list(sorted_mock.values()))
        array_tx = np.array(list(sorted_by_mock_tx.values()))
        gap_array = list()
        agg_gap = 0
        for i in range(1, len(sorted_mock)):
            cur_mock = array_mock[i-1]
            next_mock = array_mock[i]
            if next_mock == 0:
                next_mock = sys.float_info.epsilon
            cur_tx = array_tx[i - 1]
            next_tx = array_tx[i]
            if next_tx == 0:
                next_tx = sys.float_info.epsilon
            curr_gap = abs((cur_mock / next_mock) - (cur_tx / next_tx))
            gap_array.append(curr_gap)
            agg_gap += curr_gap

        j_s_d = distance.jensenshannon(array_mock, array_tx)

        if j_s_d >= 0.1:
            if (j_s_d >= 0.2) or (agg_gap >= 2.0):
                problematic_sites.append(site_name)
    return problematic_sites
