from crispector2.utils.constants_and_types import IS_EDIT, FREQ, TX_READ_NUM, MOCK_READ_NUM, ON_TARGET, CUT_SITE, AlgResult, \
    PAM_WINDOW, GRNA_WINDOW, RANDOM_NUMBER, ALLELE, SNP_PHASE, RANDOM_EDIT_READS, TX_EDIT, EDIT_PERCENT, CI_LOW, \
    CI_HIGH, REFERENCE
import numpy as np
import pandas as pd
import random
import os
import copy
from scipy.stats import norm
from crispector2.input_processing.alignment import align_random_df
from crispector2.utils.logger import LoggerWrapper
from typing import Dict
from scipy.spatial import distance
import matplotlib.pyplot as plt
from crispector2.algorithm.binomial_probability import compute_binom_p
from crispector2.algorithm.core_algorithm import CoreAlgorithm
from crispector2.modifications.modification_tables import ModificationTables
from crispector2.modifications.modification_types import ModificationTypes
from tqdm import tqdm
import sys


class Interval:
    """
    class to handle intervals of CI of results summary
    """

    def __init__(self, interval):
        self.start = interval[0]
        self.end = interval[1]

class RandomReadsHandler:
    """
    Class to handle all random reads across all sites
    """

    def __init__(self, mock_random_reads, tx_random_reads, tables_d, result_summary_d, alleles_ref_d,
                 amplicon_min_score, translocation_amplicon_min_score, enable_substitutions, confidence_interval, donor,
                 min_num_of_reads, override_noise_estimation, allele_ref_df, OUTPUT_DIR):

        self._mock_random_reads = mock_random_reads
        self._tx_random_reads = tx_random_reads
        self._tables_d = tables_d
        self._result_summary_d = result_summary_d
        self._alleles_ref_d = alleles_ref_d
        self._amplicon_min_score = amplicon_min_score
        self._translocation_amplicon_min_score = translocation_amplicon_min_score
        self._enable_substitutions = enable_substitutions
        self._confidence_interval = confidence_interval
        self._donor = donor
        self._min_num_of_reads = min_num_of_reads
        self._override_noise_estimation = override_noise_estimation
        self._allele_ref_df = allele_ref_df
        self._outdir = OUTPUT_DIR
        self._allele_site_mapping_d = dict()
        # dict to save the ratios between editing reads in alleles of the same site
        self._editing_ratios_d = dict()
        # dict to save ratios of allele's mock
        self._mock_ratios_d = dict()
        # Set logger
        logger = LoggerWrapper.get_logger()
        self._logger = logger

    def get_info_for_alleles_plots(self, alleles_details):
        info_allele_plot = dict()
        for parent_site, child_alleles in alleles_details.items():
            if len(child_alleles) == len(self._mock_ratios_d[parent_site]):
                info_allele_plot[parent_site] = dict()
                for allele_site_char, allele_site_info in child_alleles.items():
                    site_name = allele_site_info[0]
                    snvs = allele_site_info[2]
                    ref_amp = self._allele_ref_df.at[site_name, REFERENCE]
                    pam = self._allele_ref_df.at[site_name, PAM_WINDOW]
                    grna = self._allele_ref_df.at[site_name, GRNA_WINDOW]
                    freq = self._mock_ratios_d[parent_site][site_name]
                    info_allele_plot[parent_site][site_name] = {'RefAmp': ref_amp, 'frequency': freq,
                                                                'snvs': snvs, 'pam': pam, 'grna': grna}

        return info_allele_plot

    def _random_reads_statistics_and_warnings(self):
        random_reads_for_analysis_d = dict()
        for site_to_check in self._tx_random_reads.keys():
            random_reads_for_analysis_d[site_to_check] = {'total_mock': 0,
                                                          'total_tx': 0,
                                                          'random_mock': 0,
                                                          'random_tx': 0,
                                                          'random_mock_freq': 0,
                                                          'random_tx_freq': 0}
            try:
                total_site_mock_reads = sum(self._tables_d[site_to_check].mock_reads[FREQ])
                total_site_tx_reads = sum(self._tables_d[site_to_check].tx_reads[FREQ])
                random_tx_freq = sum(self._tx_random_reads[site_to_check][FREQ])
                random_mock_freq = 0
                if site_to_check in self._mock_random_reads.keys():
                    random_mock_freq = sum(self._mock_random_reads[site_to_check][FREQ])
                per_random_mock = round((random_mock_freq / total_site_mock_reads) * 100, 2)
                per_random_tx = round((random_tx_freq / total_site_tx_reads) * 100, 2)

                random_reads_for_analysis_d[site_to_check]['total_mock'] = total_site_mock_reads
                random_reads_for_analysis_d[site_to_check]['total_tx'] = total_site_tx_reads
                random_reads_for_analysis_d[site_to_check]['random_mock'] = random_mock_freq
                random_reads_for_analysis_d[site_to_check]['random_tx'] = random_tx_freq
                random_reads_for_analysis_d[site_to_check]['random_mock_freq'] = per_random_mock
                random_reads_for_analysis_d[site_to_check]['random_tx_freq'] = per_random_tx

                if per_random_mock > 1 or per_random_tx > 1:
                    self._logger.warning("Site {} - Has {}% random mock reads and {}% random tx reads.".format(
                        site_to_check, per_random_mock, per_random_tx))
            except:
                continue

        random_for_analysis_df = pd.DataFrame.from_dict(random_reads_for_analysis_d, orient='index')
        # random_for_analysis_df.to_csv(os.path.join(self._outdir, 'random_reads_stats.csv'))

    def map_alleles_to_site(self, allele_info_d):
        """
        function to create mapping between alleles to original sites
        :param allele_info_d:
        :return:
        """
        for site, alleles in allele_info_d.items():
            self._allele_site_mapping_d[site] = list()
            for snp, allele_info in alleles.items():
                self._allele_site_mapping_d[site].append(allele_info[0])

    def _calc_mock_and_editing_ratios(self):
        """
        function to calculate the editing ratios between alleles of the same site
        :return:
        """
        mock_reads_per_allele = dict()
        total_mock_per_site_d = dict()
        editing_reads_per_allele = dict()
        total_edit_per_site_d_before = dict()
        total_edit_per_site_d_after = dict()
        for site in self._alleles_ref_d.keys():
            try:
                site_name = site+'_'
                total_edit_per_site_d_before[site] = 0
                total_edit_per_site_d_after[site] = 0
                editing_reads_per_allele[site] = dict()
                total_mock_per_site_d[site] = 0
                mock_reads_per_allele[site] = dict()
                for subsite, subsite_info in self._tables_d.items():
                    if (site_name in subsite) and ('[' in subsite):
                        # for edit reads ratio
                        subsite_edit_reads = sum(subsite_info.tx_reads[subsite_info.tx_reads[IS_EDIT]==True][FREQ])
                        total_edit_per_site_d_before[site] += subsite_edit_reads
                        if subsite_edit_reads == 0:
                            subsite_edit_reads = sys.float_info.epsilon
                        total_edit_per_site_d_after[site] += subsite_edit_reads
                        editing_reads_per_allele[site][subsite] = subsite_edit_reads
                        # for mock reads ratio
                        subsite_mock_reads = sum(subsite_info.mock_reads[FREQ])
                        mock_reads_per_allele[site][subsite] = subsite_mock_reads
                        total_mock_per_site_d[site] += subsite_mock_reads
            except:
                continue

        # normalize dicts
        for site in editing_reads_per_allele.keys():
            self._editing_ratios_d[site] = dict()
            self._mock_ratios_d[site] = dict()
            edit_reads = True
            if total_edit_per_site_d_before[site] == 0:
                edit_reads = False

            if edit_reads:
                for allele_tx, edit_reads_count in editing_reads_per_allele[site].items():
                    self._editing_ratios_d[site][allele_tx] = edit_reads_count / total_edit_per_site_d_after[site]

            for allele_mock, mock_read_count in mock_reads_per_allele[site].items():
                self._mock_ratios_d[site][allele_mock] = mock_read_count / total_mock_per_site_d[site]
                if not edit_reads:
                    self._editing_ratios_d[site][allele_mock] = mock_read_count / total_mock_per_site_d[site]


    def _split_random_reads_according_to_ratios(self):
        """
        function that split the df of random reads, both mock and tx, according to ratios that calculated previously
        :return: new tables_d with random reads assignment according to ratios
        """
        new_tables_d = copy.deepcopy(self._tables_d)
        for site in self._editing_ratios_d.keys():
            random_mock = self._generate_randomness_number(self._mock_random_reads[site])
            random_tx = self._generate_randomness_number(self._tx_random_reads[site])
            # erase the snp/allele determination
            random_mock[SNP_PHASE] = ''
            random_tx[ALLELE] = ''
            # order by the random number that was generated
            random_mock = random_mock.sort_values(by=[RANDOM_NUMBER], ascending=True)
            random_tx = random_tx.sort_values(by=[RANDOM_NUMBER], ascending=True)
            # for edit ratios
            ratios_edit_list = list()
            allele_edit_list = list()
            # for mock ratios
            ratios_mock_list = list()
            allele_mock_list = list()

            # save list of edit ratios
            for i, (allele, ratio) in enumerate(self._editing_ratios_d[site].items()):
                if len(ratios_edit_list) == 0:
                    ratios_edit_list.append(ratio)
                else:
                    ratios_edit_list.append(ratio + ratios_edit_list[i-1])
                allele_edit_list.append(allele)

            # save list of mock ratios
            for i, (allele, ratio) in enumerate(self._mock_ratios_d[site].items()):
                if len(ratios_mock_list) == 0:
                    ratios_mock_list.append(ratio)
                else:
                    ratios_mock_list.append(ratio + ratios_mock_list[i-1])
                allele_mock_list.append(allele)

            # assign tx reads according to edit ratios
            new_allele_assign_tx_d = dict()
            for i, allele_name in enumerate(allele_edit_list):
                new_allele_assign_tx_d[allele_name] = dict()
                if i == 0:
                    j = 0
                else:
                    j = ratios_edit_list[i-1]
                allele_tx = random_tx[random_tx[RANDOM_NUMBER].between(j, ratios_edit_list[i])]
                allele_tx = allele_tx.astype({FREQ: np.int64})
                new_allele_assign_tx_d[allele_name] = allele_tx

            # assign mock reads according to mock ratios
            new_allele_assign_mock_d = dict()
            for i, allele_name in enumerate(allele_mock_list):
                new_allele_assign_mock_d[allele_name] = dict()
                if i == 0:
                    j = 0
                else:
                    j = ratios_mock_list[i-1]
                allele_mock = random_mock[random_mock[RANDOM_NUMBER].between(j, ratios_mock_list[i])]
                allele_mock = allele_mock.astype({FREQ: np.int64})
                new_allele_assign_mock_d[allele_name] = allele_mock

            # align mock
            aligned_mock_d_allele = align_random_df(new_allele_assign_mock_d, self._allele_ref_df,
                                                    self._amplicon_min_score, self._translocation_amplicon_min_score)
            # align tx
            aligned_tx_d_allele = align_random_df(new_allele_assign_tx_d, self._allele_ref_df,
                                                  self._amplicon_min_score, self._translocation_amplicon_min_score)
            # concat into new_tables_d (both mock and tx) and group reads
            for allele_name, allele_df in aligned_mock_d_allele.items():
                mock_w_random = pd.concat([self._tables_d[allele_name].mock_reads, allele_df])
                new_tables_d[allele_name].mock_reads = mock_w_random
                new_tables_d[allele_name].n_reads_mock = sum(new_tables_d[allele_name].mock_reads[FREQ])
            for allele_name, allele_df in aligned_tx_d_allele.items():
                tx_w_random = pd.concat([self._tables_d[allele_name].tx_reads, allele_df])
                new_tables_d[allele_name].tx_reads = tx_w_random
                new_tables_d[allele_name].n_reads_tx = sum(new_tables_d[allele_name].tx_reads[FREQ])

        return new_tables_d

    def _parse_freq_df(self, df):
        """
        function that parse the 'frequency' of read into separate rows in order to split the random reads in a correct
        manner
        :param df:
        :return: new df with line per each frequency of read
        """
        all_row_array = []
        for i, row in df.iterrows():
            frequency = row[FREQ]
            row[FREQ] = 1
            single_row = row.to_numpy()
            multiple_rows = np.tile(single_row, (frequency, 1))
            if len(all_row_array) == 0:
                all_row_array = multiple_rows
            else:
                all_row_array = np.concatenate((all_row_array, multiple_rows), axis=0)
        parsed_df = pd.DataFrame(all_row_array, columns=df.columns)

        return parsed_df

    def _generate_randomness_number(self, df):
        """
        function that generate random numbers for randomly assign the reads
        :param df:
        :return: df with additional column of random numbers range from 0 to 1
        """
        random_array = [random.random() for i in range(len(df))]
        df[RANDOM_NUMBER] = random_array

        return df

    def _compute_best_stats(self, dfs_data):
        """
        Takes all ata regarding alleles sites (random and not random reads) and compute new statistics
        :param dfs_data: tables_d. all the modifications of all sites, including the mock & tx dfs in it
        :return: The new dfs of the alleles and the new statistics scores
        """

        tx_reads_d = dict()
        mock_reads_d = dict()
        alleles_to_run = list()
        for site in dfs_data.keys():
            if '[' in site:
                alleles_to_run.append(site)
        # create tx_reads_d and mock_reads_d
        for site, site_info in dfs_data.items():
            if site in alleles_to_run:
                tx_reads_d[site] = dfs_data[site].tx_reads.reset_index()
                mock_reads_d[site] = dfs_data[site].mock_reads.reset_index()

        # Get modification types and positions priors
        modifications_allele = ModificationTypes.init_from_cfg(self._enable_substitutions)

        # Convert alignment to modification tables
        tables_d_allele: Dict[str, ModificationTables] = dict()
        for site, row in self._allele_ref_df.iterrows():
            if site in tx_reads_d.keys():
                tx_reads_num = tx_reads_d[site][FREQ].sum().astype(int)
                mock_reads_num = mock_reads_d[site][FREQ].sum().astype(int)
                if self._donor and row[ON_TARGET]:
                    pass
                elif min(tx_reads_num, mock_reads_num) < self._min_num_of_reads:
                    pass
                else:
                    tables_d_allele[site] = ModificationTables(tx_reads_d[site], mock_reads_d[site],
                                                               modifications_allele, row)
            else:
                try:
                    tables_d_allele[site] = dfs_data[site]
                except:
                    continue
                # Compute binomial coin for all modification types
        binom_p_d = compute_binom_p(tables_d_allele, modifications_allele, self._override_noise_estimation,
                                    self._allele_ref_df)
        # TBD: delete from binom the log
        # Run crispector core algorithm on all sites
        result_summary_d_allele: AlgResult = dict()  # Algorithm result dictionary
        algorithm_d: Dict[str, CoreAlgorithm] = dict()
        for site, row in self._allele_ref_df.iterrows():
            cut_site = row[CUT_SITE]
            if site in tx_reads_d.keys():
                algorithm_d[site] = CoreAlgorithm(cut_site, modifications_allele, binom_p_d[site],
                                                  self._confidence_interval, row[ON_TARGET])
                result_summary_d_allele[site] = algorithm_d[site].evaluate(tables_d_allele[site])
                result_summary_d_allele[site][ON_TARGET] = row[ON_TARGET]
                result_summary_d_allele[site][PAM_WINDOW] = self._allele_ref_df.at[site, PAM_WINDOW]
                result_summary_d_allele[site][GRNA_WINDOW] = self._allele_ref_df.at[site, GRNA_WINDOW]
            else:
                try:
                    result_summary_d_allele[site] = self._result_summary_d[site]
                except:
                    tx_reads_num = 0
                    mock_reads_num = 0
                    result_summary_d_allele[site] = {TX_READ_NUM: tx_reads_num, MOCK_READ_NUM: mock_reads_num,
                                                     ON_TARGET: row[ON_TARGET]}

        return tables_d_allele, result_summary_d_allele

    def _is_intersect(self, ci_arr):
        """
        function that returns True if some intervals are intersects and False if all not intersects.
        :param: ci_arr: array of intervals
        :return: bool: True / False intersection of intervals
        """
        n = len(ci_arr)
        ci_arr.sort(key=lambda x: x.start)
        # In the sorted array, if the start of an interval is less than end of previous one - there is an overlap
        for i in range(1, n):
            if ci_arr[i - 1].end > ci_arr[i].start:
                return True
        # Else - no overlap
        return False

    def _compute_confidence_interval(self, editing_activity, n_reads_tx, confidence_interval):
        """
        Compute confidence interval and returns low & high CI boundary
        :param editing_activity: the calculated editing activity
        :param n_reads_tx: number of treatment reads
        :param confidence_interval: the confidence interval parameter
        :return: Tuple of low & high CI boundary
        """
        confidence_inv = norm.ppf(confidence_interval + (1 - confidence_interval) / 2)
        half_len_CI = confidence_inv * np.sqrt((editing_activity * (1 - editing_activity)) / n_reads_tx)
        return max(0, editing_activity - half_len_CI), editing_activity + half_len_CI

    def _compute_new_statistics_based_on_random_assignment(self, original_n_reads_tx, original_n_reads_edited, allele_random,
                                                           sum_all, confidence_interval):
        """
        compute new high and low CI, where for high we are taking into consideration all random reads, and for low none
        :param original_n_reads_tx: the original number of treatment reads
        :param original_n_reads_edited: the original number of edited treatment reads
        :param allele_random: number of random edited reads of the specific site
        :param sum_all: sum of all random edited reads from all alleles of the same site
        :param confidence_interval: the confidence interval parameter
        :return: array of [low_CI, high_CI]
        """
        # high_confidence - consider all random edited reads from both alleles
        n_reads_tx_high_CI = original_n_reads_tx + sum_all - allele_random
        n_reads_edited_high_CI = original_n_reads_edited + sum_all - allele_random
        editing_activity_high_CI = (n_reads_edited_high_CI / n_reads_tx_high_CI)
        _, CI_high = self._compute_confidence_interval(editing_activity_high_CI, n_reads_tx_high_CI, confidence_interval)
        # low_confidence - subtract random edited reads
        n_reads_tx_low_CI = original_n_reads_tx - allele_random
        n_reads_edited_low_CI = original_n_reads_edited - allele_random
        editing_activity_low_CI = n_reads_edited_low_CI / n_reads_tx_low_CI
        CI_low, _ = self._compute_confidence_interval(editing_activity_low_CI, n_reads_tx_low_CI, confidence_interval)

        return [CI_low * 100, CI_high * 100]

    def _estimate_random_reads_editing_effect(self, dict_of_alleles, confidence_interval):
        """
        function that checks if site has alleles that CI-overlap after manipulation of:
            high_CI = all random edited reads assign to an allele
            low_CI = none random edited reads assign to an allele
        :param dict_of_alleles: dictionary of all alleles
        :param confidence_interval: the confidence interval parameter
        :return: True if the alleles of a site are overlapping, and False if not
        """
        new_alleles_CI = list()
        try:
            sum_randoms = sum(allele[RANDOM_EDIT_READS] for allele in dict_of_alleles.values())
        except:
            # read count per allele is too low that there is no statistics for it
            return False

        for allele in dict_of_alleles.values():
            n_allele_random = allele[RANDOM_EDIT_READS]
            n_reads_tx = allele[TX_READ_NUM]
            n_reads_edited = allele[TX_EDIT]
            new_allele = self._compute_new_statistics_based_on_random_assignment(n_reads_tx, n_reads_edited,
                                                                            n_allele_random, sum_randoms,
                                                                            confidence_interval)
            new_alleles_CI.append(Interval(new_allele))
        if self._is_intersect(new_alleles_CI):
            return True
        return False

    def _get_medoid(self, d, sites_linkage):
        """
        function that takes all iterations of summary results and select the medoid of the data with respect to
        each site's alleles
        :param d: the dictionary with all statistics results of all iterations
        :param sites_linkage: a key-value of all sites and their alleles
        :return: The medoid of the results. medoid per site (same alleles in the site would get same medoid)
        """

        best_index_per_allele = dict()

        for site, alleles in sites_linkage.items():
            # get all results as data points in the dimension of number of alleles
            data_points = list()
            for indx, results in enumerate(d):
                sub_point = list()
                for allele_name in alleles:
                    try:
                        editing_activity = results[allele_name][EDIT_PERCENT]
                    except:  # for read count less than threshold
                        editing_activity = 0
                    sub_point.append(editing_activity)
                data_points.append(sub_point)

            # creating an empty distances matrix
            dict_len = len(d)
            dist_matrix = np.empty(shape=(dict_len, dict_len), dtype='object')

            # calculate each pairwise euclidean distance
            for i, fix_point in enumerate(data_points):
                for j, var_point in enumerate(data_points):
                    dist_matrix[i, j] = distance.euclidean(fix_point, var_point)

            # get the smallest aggregated distance as the medoid and as the selected point.
            medoid = np.argmin(dist_matrix.sum(axis=0))

            for allele_name in alleles:
                best_index_per_allele[allele_name] = medoid

            # enable plotting medoid
            # # plotting the data points
            # if len(alleles) == 2:
            #     x = [row[0] for row in data_points]
            #     y = [row[1] for row in data_points]
            #     label = [idx for idx in range(len(data_points))]
            #     plt.scatter(x, y)

            #     for i, txt in enumerate(label):
            #         plt.annotate(txt, (x[i], y[i]))
            #     plt.title(f'The medoid of site {site} is: {medoid}')
            #     plt.savefig(os.path.join(self._outdir, f'medoid_{site}.png'), bbox_inches='tight')
            #     plt.close()
            #     plt.clf()

        return best_index_per_allele

    def _estimated_upper_lower_CI_ea(self, result_summary_d_dict):
        """
        function that retrieve the highest CI score out of all runs and lowest, per each site
        :param result_summary_d_dict:
        :return: dictionary of lowest/highest CI for each allele site
        """
        highest_ea = dict()
        lowest_ea = dict()

        for rep in result_summary_d_dict:
            for site_name, site_info in rep.items():
                if '[' in site_name:
                    try:
                        if site_name not in highest_ea:
                            highest_ea[site_name] = site_info[CI_HIGH]
                            lowest_ea[site_name] = site_info[CI_LOW]
                        else:
                            if site_info[CI_HIGH] > highest_ea[site_name]:
                                highest_ea[site_name] = site_info[CI_HIGH]

                            if site_info[CI_LOW] < lowest_ea[site_name]:
                                lowest_ea[site_name] = site_info[CI_LOW]
                    except:
                        continue

        df_low = pd.DataFrame.from_dict(lowest_ea, orient='index', columns=['low'])
        df_high = pd.DataFrame.from_dict(highest_ea, orient='index', columns=['high'])

        df_final = pd.concat([df_low, df_high], axis=1)
        df_final.to_csv(os.path.join(self._outdir, 'CI_upper_lower_estimation.csv'))

        # return highest_ea, lowest_ea

    def _re_calculate_statistics(self, outdir):
        """
        function that takes all ambiguous sites (regarding CI), re-assign the random reads to the alleles and re-compute
        the statistics. This done 11 times, taking the median score of each site as the "best" score
        :param outdir: path of out directory
        :return: The best selected statistics and dfs out of 11 random possibles
        """

        tables_d_dict = list()
        result_summary_d_dict = list()
        for i in tqdm(range(11)):
            # assign "random reads" randomly between sites
            tables_d_random_temp = self._split_random_reads_according_to_ratios()

            # compute new statistics per each site that was modified randomly
            tables_d_allele, result_summary_d_allele = self._compute_best_stats(tables_d_random_temp)

            # store the dfs and the statistics results
            tables_d_dict.append(tables_d_allele)
            result_summary_d_dict.append(result_summary_d_allele)

        # calculate upper and lower estimated CI editing activity
        self._estimated_upper_lower_CI_ea(result_summary_d_dict)
        # calculate the medoid of each site to set as the best result
        best_index_per_allele = self._get_medoid(result_summary_d_dict, self._allele_site_mapping_d)

        # assign the "best score" into new dictionaries
        best_tables_d = dict()
        best_results_summary = dict()
        for allele, medoid_index in best_index_per_allele.items():
            try:
                best_tables_d[allele] = tables_d_dict[medoid_index][allele]
                best_results_summary[allele] = result_summary_d_dict[medoid_index][allele]
            except:
                continue

        return best_tables_d, best_results_summary

    def run_initial_assigner(self):
        self._random_reads_statistics_and_warnings()
        self._calc_mock_and_editing_ratios()
        for site in self._editing_ratios_d.keys():
            random_mock = self._parse_freq_df(self._mock_random_reads[site])
            self._mock_random_reads[site] = random_mock
            random_tx = self._parse_freq_df(self._tx_random_reads[site])
            self._tx_random_reads[site] = random_tx
        tables_d_new = self._split_random_reads_according_to_ratios()
        tables_d_allele, result_summary_d_allele = self._compute_best_stats(tables_d_new)

        return tables_d_allele, result_summary_d_allele

    def run_best_assigner(self, result_summary_d):
        re_run_overlapping_sites = dict()
        # get the list of all the alleles of site
        for site, allele_names in self._allele_site_mapping_d.items():
            # get the results_summary_stat_per all sites' alleles
            alleles_results_summary = dict()
            for allele_name in allele_names:
                alleles_results_summary[allele_name] = result_summary_d[allele_name]
            # if the quality statistics is poor:
            if self._estimate_random_reads_editing_effect(alleles_results_summary, self._confidence_interval):
                re_run_overlapping_sites[site] = allele_names

        # find the best statistics for the overlapping sites
        self._logger.info("Start re-evaluating alleles...")
        new_tables_d, best_stats_per_site = None, None
        if len(re_run_overlapping_sites) > 0:
            new_tables_d, best_stats_per_site = self._re_calculate_statistics(self._outdir)

        return new_tables_d, best_stats_per_site

