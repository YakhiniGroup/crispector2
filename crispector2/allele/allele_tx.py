from crispector2.utils.constants_and_types import FREQ, IS_RANDOM, AVG_SCORE, SITE_NAME, ALLELE, ALIGN_SCORE, READ
import numpy as np
import pandas as pd
import math
import random
from crispector2.utils.configurator import Configurator
from crispector2.input_processing.alignment import LocalLooseAlignment
from crispector2.utils.logger import LoggerWrapper

# np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)


def _is_long_del(matches, pos, snp_locus):
    char = matches[pos]
    counter = 0
    while char == '-':
        counter += 1
        pos += 1
        char = matches[pos]
    if counter > len(snp_locus):
        return True, counter
    else:
        return False, counter


class AlleleForTx:
    """
    Class to handle the allele case of the treatment
    """

    def __init__(self, tx_df, new_alleles, alleles_ref_reads, alleles_ratios):
        self._tx_df = tx_df
        self.new_alleles = new_alleles
        # create list of sites to be handled
        self._sites_to_be_alleles = list(self.new_alleles.keys())
        self._new_tx_df = dict()
        # self._sites_score = pd.DataFrame(columns=[SITE_NAME, AVG_SCORE])
        self.alleles_ref_reads = alleles_ref_reads
        self.random_reads = dict()
        self.alleles_mock_tx_ratios = alleles_ratios

        # Set logger
        logger = LoggerWrapper.get_logger()
        self._logger = logger

        self._cfg = Configurator.get_cfg()
        # create alignment instances
        self._ll_aligner = LocalLooseAlignment(self._cfg["local_loose_alignment"])

    # -------------------------------#
    ######### Public methods #########
    # -------------------------------#

    def run(self, ratios_df):
        # original_sites = dict()
        for tx_site_name, tx_site_df in self._tx_df.items():
            # if tx site has an allele site
            if tx_site_name in self._sites_to_be_alleles:
                # prepare relevant information for execution
                self._new_tx_df[tx_site_name] = dict()
                # these modifications changing the tx_reads_d from the main
                self.alleles_mock_tx_ratios[tx_site_name]['tx'] = dict()
                tx_site_df[ALLELE] = None
                tx_site_df[IS_RANDOM] = False
                tx_site_df[ALIGN_SCORE] = math.nan
                # list of all relevant allele dfs for this site
                list_sites_of_allele = self.new_alleles[tx_site_name]
                allele_relevant_ref_reads = self.alleles_ref_reads[tx_site_name]

                # iterate over all rows in df and assign to relevant allele
                for i_d, row in tx_site_df.iterrows():
                    seq = row[READ]
                    temp_scores = list()
                    is_random = False

                    for allele_type, allele_ref in allele_relevant_ref_reads.items():
                        snp_loci = list_sites_of_allele[allele_type][2]
                        alignment = self._ll_aligner._align_seq_to_read(allele_ref, seq)
                        alignment_score = alignment.score
                        [ref_aligned, matches, read_aligned, _] = format(alignment[0]).split("\n")

                        # compensate the situation where the deletion is in the reference thus it won't take it
                        # into consideration in the alignment score
                        # if '-' in set(ref_aligned).intersection(set(read_aligned)):
                        #     max_len = len(matches)
                        #     j = 0
                        #     while j < max_len:
                        #         if (matches[j] == '-') and (ref_aligned[j] == '-') and (read_aligned[j] == '-'):
                        #             long_del, counter = _is_long_del(matches, j,
                        #                                                   list_sites_of_allele[allele_type][2])
                        #             if not long_del:
                        #                 alignment_score += self._ll_aligner._aligner.match_score
                        #                 j += 1
                        #             elif long_del:
                        #                 j += counter
                        #         else:
                        #             j += 1

                        if '-' in set(ref_aligned).intersection(set(read_aligned)):
                            max_len = len(matches)
                            j = 0
                            mismatch_counter = 0
                            snps_in_interval = 0
                            real_del_counter = 0
                            while j < max_len:
                                k = j

                                while matches[k] == '-':
                                    if (matches[k] == '-') and (ref_aligned[k] == '-') and (read_aligned[k] == '-'):
                                        real_del_counter += 1
                                        k += 1
                                    else:
                                        mismatch_counter += 1
                                        k += 1

                                if (mismatch_counter > 0) or (real_del_counter > 0):
                                    # deletion_interval = [k - mismatch_counter - real_del_counter, k]
                                    # for snp in snp_loci:
                                    #     if (snp >= deletion_interval[0]) and (snp < deletion_interval[1]):
                                    #         snps_in_interval += 1
                                    if (mismatch_counter == 0) and (real_del_counter > 0):
                                        alignment_score += (self._ll_aligner._aligner.match_score * real_del_counter)

                                    j = k
                                    mismatch_counter = 0
                                    snps_in_interval = 0
                                    real_del_counter = 0
                                else:
                                    j += 1

                        temp_scores.append(alignment_score)
                    max_score = np.max(temp_scores)
                    if temp_scores.count(max_score) == 1:
                        best_allele = list(allele_relevant_ref_reads.keys())[np.argmax(temp_scores)]
                    else:
                        is_random = True
                        best_allele = random.choice(list(allele_relevant_ref_reads.keys()))

                    tx_site_df.loc[i_d, ALLELE] = list_sites_of_allele[best_allele][0]
                    tx_site_df.loc[i_d, ALIGN_SCORE] = max_score
                    tx_site_df.loc[i_d, IS_RANDOM] = is_random

                # separate to random reads and not random
                self.random_reads[tx_site_name] = tx_site_df[tx_site_df[IS_RANDOM] == True]
                df_no_random = tx_site_df[tx_site_df[IS_RANDOM] == False]
                # sum reads without random reads
                site_sum_reads = sum(df_no_random[FREQ])
                # set threshold for percentage of random reads out of all reads. If more than that - do not report alleles
                if sum(tx_site_df[tx_site_df[IS_RANDOM] == True][FREQ]) / sum(tx_site_df[FREQ]) > 0.5:  # TBD: make hyperparameter (0.5)
                    self._logger.info("Site {} has to many random reads".format(tx_site_name))

                ratios_dict = dict()
                # set the new df for the tx allele
                for allele_type, allele_info in list_sites_of_allele.items():
                    # get the tx ratios for this site
                    tx_site_df_temp = df_no_random.loc[df_no_random[ALLELE] == allele_info[0]]
                    sum_reads = sum(list(tx_site_df_temp[FREQ]))
                    # append information to ratios
                    self.alleles_mock_tx_ratios[tx_site_name]['tx'][allele_type] = sum_reads / site_sum_reads
                    ratios_dict[allele_type] = sum_reads
                    # calculate the alignment average score. If less than one [TBD]
                    avg_score = tx_site_df_temp[ALIGN_SCORE].mean(axis=0, skipna=True)
                    # self._sites_score = self._sites_score.append({SITE_NAME: allele_info[0], AVG_SCORE: avg_score},
                    #                                              ignore_index=True)
                    tx_site_df_temp = tx_site_df_temp.drop(labels=[ALLELE, ALIGN_SCORE], axis=1)
                    self._new_tx_df[tx_site_name][allele_type] = [allele_info[0], tx_site_df_temp]
                # original_sites[tx_site_name] = tx_site_df  # TBD: maybe delete

                # normalize the dict
                try:
                    sum_all_reads = sum(ratios_dict.values())
                    for _i in ratios_dict.keys():
                        ratios_dict[_i] /= sum_all_reads
                except:
                    print(f'problem with {tx_site_name}')

                # append to df ratios in order to keep track the ratios with compare to mock- TBD: return later
                ratios_df.at[list(ratios_df['site_name']).index(tx_site_name), 'tx_ratios'] = ratios_dict

        return self._new_tx_df, ratios_df  # original_sites # , round(self._sites_score, 3)