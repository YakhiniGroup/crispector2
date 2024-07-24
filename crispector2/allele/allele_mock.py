from crispector2.utils.constants_and_types import HALF_WINDOW_LEN, LENGTH_RATIO, SNP_NUC_TYPE, FREQ, ALIGNMENT_W_DEL, LEN, \
    SNP_PHASE, IS_RANDOM, BASE, REFERENCE, ALIGN_CUT_SITE, N, SCORE_NORM_RATE, READ
import numpy as np
import pandas as pd
from scipy.stats import entropy
import scipy.stats
import math
import random
from crispector2.input_processing.alignment import LocalStrictAlignment, LocalLooseAlignment, GlobalStrictAlignment, MainAlignment
from crispector2.utils.configurator import Configurator
from crispector2.utils.logger import LoggerWrapper
from collections import Counter


class AlleleForMock:
    """
    Class to handle the allele case of the Mock
    """
    def __init__(self, ratios, ref_df, max_poten_snvs, max_allele_mismatches, max_len_snv_ctc, random_percentage):
        self._ref_df = ref_df
        self._ratios = ratios
        self.new_alleles = dict()
        self.alleles_ref_reads = dict()
        self._consensus_len = None  # set the consensus of the majority length of the reads
        self._number_of_alleles = 0  # determine the number of alleles
        self._length_ratio = LENGTH_RATIO
        self._half_window_len = HALF_WINDOW_LEN
        self._max_potential_snvs = max_poten_snvs
        self._max_mm = max_allele_mismatches
        self._max_len_snv_ctc = max_len_snv_ctc
        self._random_reads_threshold_percentage = random_percentage
        self.alleles_mock_tx_ratios = dict()

        # dict for random reads
        self.random_reads = dict()
        # df that saves the ratios of nb over the whole process - mock before, after and tx
        self.df_mock_tx_snp_ratios = pd.DataFrame(
            columns=['site_name', 'mock_ratios_before', 'mock_ratios_after',
                     'tx_ratios', 'number_of_alleles'])

        # set the lower decision bound for significant entropy
        self._entropy_lower_bound_score = 0
        for p in self._ratios:
            self._entropy_lower_bound_score -= p * math.log(p, len(self._ratios))

        # Set logger
        logger = LoggerWrapper.get_logger()
        self._logger = logger

        self._cfg = Configurator.get_cfg()
        # create alignment instances
        self._ls_aligner = LocalStrictAlignment(self._cfg["local_strict_alignment"])
        self._ll_aligner = LocalLooseAlignment(self._cfg["local_loose_alignment"])
        self._gs_aligner = GlobalStrictAlignment(self._cfg["global_strict_alignment"])

        # create alignment instance as the main alignment
        self._main_aligner = MainAlignment(self._cfg["alignment"], self._cfg["NHEJ_inference"]["window_size"])

    # -------------------------------#
    ######### Public methods #########
    # -------------------------------#

    def run(self, site_name, df):
        self._site_name = site_name
        self._ref_amplicon = self._ref_df.at[site_name, REFERENCE]
        self._df = df
        self._snp_locus = list()  # holds the position of all snp loci in this site
        self._cut_site = None
        self._total_num_reads = 0
        self._nuc_distribution_before = None  # holds all distributions of all snps
        self._nuc_distribution_after = None
        self._windows = dict()  # holds all window that open around each snp in this site

        # get all snps and a new df relevant to the snps
        df, SNP = self._mock_SNP_detection()

        # if the function finds SNPs
        if SNP:
            self.alleles_mock_tx_ratios[self._site_name] = dict()
            # separate to random reads and not random
            self.random_reads[self._site_name] = df[df[IS_RANDOM] == True]
            df_no_random = df[df[IS_RANDOM] == False]
            sum_no_rand_reads = sum(df_no_random[FREQ])
            # set threshold for percentage of random reads out of all reads. If more than that - do not report alleles
            if sum(df[df[IS_RANDOM] == True][FREQ]) / sum(df[FREQ]) > self._random_reads_threshold_percentage:
                self._logger.info("Site {} has to many random reads".format(self._site_name))
            else:
                self.alleles_mock_tx_ratios[self._site_name]['mock'] = dict()
                # iterate over all possible alleles
                for allele in list(self._windows.keys()):
                    # filter df for reads with the current snp and clean columns
                    df_for_curr_allele = df_no_random.loc[df_no_random[SNP_PHASE] == allele]
                    # get the ratio of this allele, without the random reads
                    self.alleles_mock_tx_ratios[self._site_name]['mock'][allele] = sum(
                        df_for_curr_allele[FREQ]) / sum_no_rand_reads
                    # get the most relevant amplicon for reference from df
                    amplicon = self.alleles_ref_reads[site_name][allele]
                    df_for_curr_allele = df_for_curr_allele.drop(labels=[SNP_PHASE, LEN, SNP_NUC_TYPE], axis=1)
                    allele_window_list = self._windows[allele]
                    # TBD: change the name to more relevant one
                    _new_name = self._site_name + '_' + str(self._snp_locus) + '_' + allele

                    # add list of:(new name, filtered df, snp positions, windows ,amplicon)
                    if self._site_name not in self.new_alleles.keys():
                        self.new_alleles[self._site_name] = dict()
                    self.new_alleles[self._site_name][allele] = [  # TBD: maybe some values unnecessary. check
                        _new_name,
                        df_for_curr_allele,
                        self._snp_locus,
                        allele_window_list,
                        amplicon
                    ]

                # insert ratios to the self variable. TBD: move to the end and normalized by number of alleles
                curr_mock_df = pd.DataFrame(data=[[self._site_name, self._nuc_distribution_before,
                                                   self._nuc_distribution_after, None, self._number_of_alleles]],
                                            columns=['site_name', 'mock_ratios_before',
                                                     'mock_ratios_after', 'tx_ratios', 'number_of_alleles'])
                self.df_mock_tx_snp_ratios = pd.concat([self.df_mock_tx_snp_ratios, curr_mock_df], ignore_index=True)

        return self.new_alleles

    # -------------------------------#
    ######### Private methods #######
    # -------------------------------#

    def _alignment_to_return_reads(self, relevant_read, window, CTC, start_from):
        """
         Return the start and end coordinates of the alignment`
         :param: relevant_read: The relevant read to align to
         :param: window: The window around the snp
         :param: CTC: Does the snp close to cut-site?
         :param: start_from: start align from this point
         :return: start: Start coordinate of the alignment;
                  end: End coordinate of the alignment
         """
        if not CTC:
            alignment = self._ls_aligner._align_seq_to_read(relevant_read[start_from:], window)
        else:  # TBD: I think that it does not matter for this part if it is CTC or not!
            alignment = self._ll_aligner._align_seq_to_read(relevant_read[start_from:], window)

        [read_aligned, matches, window_aligned, _] = format(alignment[0]).split("\n")
        # finding positions of start and end in the original read
        start1 = matches.find('|')
        start2 = matches.find('.')
        end1 = matches.rfind('|')
        end2 = matches.rfind('.')
        start = min(start1, start2) + start_from
        end = max(end1, end2) + start_from + 1

        return start, end

    def _get_specific_windows(self, df):  # reference_read
        """
         Create list of windows to be aligned to reads with different length then consensus, and create self._window
         :param: reference_read: The reference read to extract window from
         :param: cut_site: The consensus cut-site
         :return: windows_list: List of all windows for further alignment
         """

        windows_dict = dict()
        reads_per_allele = dict()
        allele_phases = list(self._nuc_distribution_before.keys())[:self._number_of_alleles]
        # the specific case - make a window with correct snp for each allele for each snp position
        for allele in allele_phases:
            df_for_curr_allele = df.loc[df[SNP_PHASE] == allele]
            new_reference_read = df_for_curr_allele.iloc[0, 3]
            reads_per_allele[allele] = new_reference_read

        # check that all nb are equal except the snps
        for i in range(len(list(reads_per_allele.values())[0])):
            if i not in self._snp_locus:
                temp_nb_l = list()
                for allele_ref in list(reads_per_allele.values()):
                    temp_nb_l.append(allele_ref[i])
                if len(set(temp_nb_l)) >= 2:
                    most_freq_nb = Counter(temp_nb_l).most_common(1)
                    if most_freq_nb[0][1] > 1:
                        selected_nb = most_freq_nb[0][0]
                    else:
                        selected_nb = temp_nb_l[0]

                    for allele in allele_phases:
                        temp_allele_ref = reads_per_allele[allele]
                        correct_allele_ref = temp_allele_ref[:i] + selected_nb + temp_allele_ref[i+1:]
                        reads_per_allele[allele] = correct_allele_ref

        consensus_read = list(reads_per_allele.values())[0]
        for allele, read in reads_per_allele.items():
            window_list_per_allele = self._prepare_snp_windows(consensus_read, allele)
            windows_dict[allele] = window_list_per_allele

        return windows_dict, reads_per_allele

    def _get_general_windows(self, reference_read):
        """
         Create list of windows to be aligned to reads with different length then consensus, and create self._window
         :param: reference_read: The reference read to extract window from
         :param: cut_site: The consensus cut-site
         :return: windows_list: List of all windows for further alignment
         """
        new_reference_read = reference_read

        # Assigning all SNPs to be N in the reference read
        for snp_locus in self._snp_locus:
            new_reference_read = new_reference_read[:snp_locus] + N + new_reference_read[snp_locus + 1:]
        windows_lists = self._prepare_general_snp_windows(new_reference_read)

        return windows_lists

    def _prepare_general_snp_windows(self, ref_read):
        """
         for each snp, create all possible windows size self._half_window_len for each type of snp
         :param: reference_read: the relevant site's amplicon
         :return: windows_list: list of lists that each list contains the window around the snp,
                                if it's CTC (close to cut), and two variables that help the length of the window to be
                                equal, even if the snp is too close to beginning or end of the amplicon
         """
        windows_list_for_snp = list()
        for i, snp_locus in enumerate(self._snp_locus):
            add_nb_before = 0
            add_nb_after = 0
            if abs(self._cut_site - snp_locus) < self._half_window_len:
                CTC = True
            else:
                CTC = False  # CTC = Close To Cut-site

            # if the read is close to the end, add from the left part some more nb
            if len(ref_read) - snp_locus < self._half_window_len + 1:
                add_nb_before = self._half_window_len - (len(ref_read) - snp_locus) + 1
                window_search = ref_read[
                                snp_locus - self._half_window_len - add_nb_before:snp_locus + 1 + self._half_window_len]
                windows_list_for_snp.append((window_search, CTC, add_nb_before, add_nb_after))

            # if the read is close to the start, add from the right part some more nb
            elif snp_locus - self._half_window_len < 0:
                add_nb_after = self._half_window_len - snp_locus
                window_search = ref_read[:snp_locus + 1 + self._half_window_len + add_nb_after]
                windows_list_for_snp.append((window_search, CTC, add_nb_before, add_nb_after))
            # if "normal"
            else:
                window_search = ref_read[snp_locus - self._half_window_len:snp_locus + 1 + self._half_window_len]
                windows_list_for_snp.append((window_search, CTC, add_nb_before, add_nb_after))

        return windows_list_for_snp

    def _prepare_snp_windows(self, ref_read, allele):
        """
         for each snp, create all possible windows size self._half_window_len for each type of snp
         :param: reference_read: the relevant site's amplicon
         :return: windows_list: list of lists that each list contains the window around the snp,
                                if it's CTC (close to cut), and two variables that help the length of the window to be
                                equal, even if the snp is too close to beginning or end of the amplicon
         """
        windows_list_for_snp = list()
        for i, snp_locus in enumerate(self._snp_locus):
            add_nb_before = 0
            add_nb_after = 0
            if abs(self._cut_site - snp_locus) < self._half_window_len:
                CTC = True
            else:
                CTC = False  # CTC = Close To Cut-site

            # if the read is close to the end, add from the left part some more nb
            if len(ref_read) - snp_locus < self._half_window_len + 1:
                add_nb_before = self._half_window_len - (len(ref_read) - snp_locus) + 1
                window_search = ref_read[snp_locus - self._half_window_len - add_nb_before:snp_locus] + \
                                allele[i] + \
                                ref_read[snp_locus+1:snp_locus + 1 + self._half_window_len]
                windows_list_for_snp.append((window_search, CTC, add_nb_before, add_nb_after))

            # if the read is close to the start, add from the right part some more nb
            elif snp_locus - self._half_window_len < 0:
                add_nb_after = self._half_window_len - snp_locus
                window_search = ref_read[:snp_locus] + \
                                allele[i] + \
                                ref_read[snp_locus + 1:snp_locus + 1 + self._half_window_len + add_nb_after]
                windows_list_for_snp.append((window_search, CTC, add_nb_before, add_nb_after))
            # if "normal"
            else:
                window_search = ref_read[snp_locus - self._half_window_len:snp_locus] + \
                                allele[i] + \
                                ref_read[snp_locus+1:snp_locus + 1 + self._half_window_len]
                windows_list_for_snp.append((window_search, CTC, add_nb_before, add_nb_after))

        return windows_list_for_snp

    def _return_reads_to_nuc_dist(self, df, filtered_df, windows_list):
        """
         Return reads with different length than the consensus to the df
         :param: df: The df with the reads with the same length
         :param: filtered_df: The df of the reads that were filtered due to different length of the consensus
         :param: windows_list: list of all windows to be aligned to each read
         :return: df: The df with some returned reads (that were able to be returned);
                  df_dropped: Reads that were not able to be returned
         """
        df_dropped_list = list()

        for i, row in filtered_df.iterrows():
            position_passed = [0]
            read = row[ALIGNMENT_W_DEL]
            for idx, (window, CTC, additional_nb_before, additional_nb_after) in enumerate(windows_list):
                # get the start and end coordinates according to the alignment
                start, end = self._alignment_to_return_reads(read, window, CTC, position_passed[-1])
                position_passed.append(max(0, end - len(window)))
                # if succeed to align properly # TBD: Check the case of CTC==True
                if end - start == len(window):
                    snp_nb = read[start:end][
                        self._half_window_len + additional_nb_before - additional_nb_after]
                    filtered_df.at[i, SNP_NUC_TYPE].append(snp_nb)
                else:
                    # only rows that all snps have been found # TBD: check with Zohar if ok
                    df_dropped_list.append(row)
                    filtered_df = filtered_df.drop(i)
                    break

        df_dropped = pd.DataFrame(data=df_dropped_list, columns=list(filtered_df.columns))
        df = pd.concat([df, filtered_df])

        return df, df_dropped

    def _find_all_potential_snp(self, df):
        """
         Find all potential snp among the filtered df
         :param df: The df with the additional information regarding snps, filtered with the same read length
         :return: df: return the same df as self._df however with additional columns regarding the snps;
                  SNP: True/False to indicate if SNP is exists or not
         """
        # obtain entropies for the site run
        local_score = []
        freq_list = df[FREQ]

        # split the DataFrame into separated columns to calculate snps
        parsed_df = df[ALIGNMENT_W_DEL].apply(lambda x: pd.Series(list(x)))  # TBD: expensive in time

        # iterate over all columns, calc entropy and append to local score
        for i in range(len(parsed_df.columns)):
            temp_df = pd.concat([parsed_df[i], freq_list], axis=1, join='inner')
            temp_score, _ = calc_entropy(temp_df)
            local_score.append(temp_score)

            # if score is higher >> it is a snp
            if temp_score >= self._entropy_lower_bound_score:
                self._snp_locus.append(i)
            else:
                continue

        # get all snps to column SNP_NUC_TYPE
        for i_row, row in df.iterrows():
            for snp_loc in self._snp_locus:
                df.at[i_row, SNP_NUC_TYPE].append(str(row[ALIGNMENT_W_DEL][snp_loc]))

        return df

    def _norm_nuc_distribution(self, df):
        """
         Normalize a given df with frequency
         :param: df: The site df with the returned reads
         :return: nuc_dict: normalized nuc dictionary distribution
         """
        # create dictionary and assign values inside
        nuc_dict = dict()
        for i, row in df.iterrows():
            phase = df.at[i, SNP_PHASE]
            freq = df.at[i, FREQ]
            if phase in nuc_dict.keys():
                nuc_dict[phase] += freq
            else:
                nuc_dict[phase] = freq
        # making sure that the dictionary is sorted
        nuc_dict = dict(sorted(nuc_dict.items(), key=lambda item: item[1], reverse=True))
        # normalized the dictionary of the nucleotide distribution
        self._total_num_reads = sum(nuc_dict.values())
        for snp in nuc_dict.keys():
            nuc_dict[snp] /= self._total_num_reads
        return nuc_dict

    def _get_num_alleles(self):
        """
         Calculate the optimal number of alleles in a given distribution
         :param: self: Get some self variables
         :return: coverage_list: The list of coverage;
                  entropy_scores: The scores of all entropies of all number of potential alleles;
                  norm_scores_list: Score after normalization;
                  number_of_alleles: The optimum number of alleles
         """
        entropy_scores = list()
        norm_scores_list = list()
        coverage_list = list(self._nuc_distribution_before.values())

        for i in range(1, len(coverage_list)):
            curr_vals = coverage_list[:i + 1]
            _entropy = scipy.stats.entropy(pd.Series(curr_vals), base=i + 1)
            entropy_scores.append(_entropy)
            norm_score = score_normalization(_entropy, sum(curr_vals))
            norm_scores_list.append(norm_score)

        number_of_alleles = np.argmax(norm_scores_list) + 2

        return coverage_list, entropy_scores, norm_scores_list, number_of_alleles

    def _get_ref_read(self, most_freq_read):
        """
         Generates the new reference read for the allele
         :param: most_freq_read: The most frequent read in the allele's df
         :return: new_ref_read: The new reference read of the allele
         """
        alignment = self._ll_aligner._align_seq_to_read(self._ref_amplicon, most_freq_read)
        [ref_w_i, signs, read_w_d, _] = format(alignment[0]).split("\n")

        return read_w_d

    def _mock_SNP_detection(self):
        """
        Find all snps and return df with the snp bases
        :return: df: return the same df as self._df however with additional columns regarding the snps,
                SNP: True/False to indicate if SNP is exists or not
        """
        SNP = False
        false_df = False
        try:
            df = self._df
            df.loc[:, SNP_NUC_TYPE] = [list() for _ in range(len(df.index))]  # create snp list for each row
            df = df.sort_values([FREQ], ascending=False)  # sort according to frequency

            # adding length of each read to each row
            len_reads = list(df[ALIGNMENT_W_DEL].str.len())
            df[LEN] = len_reads

            # compute the distribution of the length of the reads
            lengths_ratio, len_dict = compute_lengths_distribution(df, set(len_reads))

            # if the highest freq of length is smaller than (1 - self._length_ratio) than don't proceed with allele
            if lengths_ratio[0] < (1 - self._length_ratio):
                self._logger.info(
                    "Site {} has ambiguous read lengths, thus won't be analyzed as alleles".format(self._site_name))
                return false_df, SNP

            # calc the majority number of reads by length and filter Dataframe based on it
            self._consensus_len = max(len_dict, key=len_dict.get)
            mask = df[LEN] == self._consensus_len
            same_len_df = df.loc[mask]
            filtered_reads_diff_len = df.loc[~mask]

            # find all potential stand-alone snps
            same_len_df = self._find_all_potential_snp(same_len_df)  # TBD: takes long time - fix
            self._cut_site = get_most_frq_cut_site(same_len_df)

            # if there is more than 1 snv however not more than self._max_potential_snvs:
            if (len(self._snp_locus) > 0) and (len(self._snp_locus) < self._max_potential_snvs):
                SNP = True
                same_len_indexes = list(same_len_df.index)

                # handle reads that are not at length of the self._consensus length:
                # get the relevant reference read
                reference_read_w_al = same_len_df.loc[same_len_indexes[list(same_len_df[FREQ]).index(
                    max(same_len_df[FREQ]))], ALIGNMENT_W_DEL]
                # create all relevant "windows" around each potential SNP
                windows_list_for_realignment = self._get_general_windows(reference_read_w_al)
                # assign reads with different length back to the main df
                df_w_returned, df_dropped = self._return_reads_to_nuc_dist(
                    same_len_df, filtered_reads_diff_len, windows_list_for_realignment)

                # phasing the multi-snp into one phase instead of a list of snps
                df_w_returned[SNP_PHASE] = df_w_returned[SNP_NUC_TYPE].apply(lambda x: ''.join(x))
                df_dropped[SNP_PHASE] = df_dropped[SNP_NUC_TYPE].apply(lambda x: ''.join(x))
                # normalized the snp distribution
                self._nuc_distribution_before = self._norm_nuc_distribution(df_w_returned[[SNP_PHASE, FREQ]])
                # get the number of alleles
                _, _, _, self._number_of_alleles = self._get_num_alleles()
                self._windows, reads_per_allele = self._get_specific_windows(df_w_returned)
                self.alleles_ref_reads[self._site_name] = reads_per_allele

                # assign alleles that are not represented to alleles that do
                df_complete = self._determine_allele_for_unknowns(df_w_returned, df_dropped)

                # normalizing again, after determine unknown reads' allele
                self._nuc_distribution_after = self._norm_nuc_distribution(df_complete[[SNP_PHASE, FREQ]])

                return df_complete, SNP
            else:
                return false_df, SNP
        except:
            return false_df, SNP

    def _determine_allele_for_unknowns(self, df_w_returned, df_dropped):
        """
        Find all snps and return df with the snp bases
        :param: df_w_returned: The main df that contains most of the reads, that classified to alleles
        :param: df_dropped: The df with reads that couldn't be assigned with allele on the straight-forward process,
                            and different approach is needed
        :return: df: return the same df as self._df however with additional columns regarding the snps,
                SNP: True/False to indicate if SNP is exists or not
        """
        # one_snp is true if one snp only found
        one_snp = False
        if len(self._snp_locus) == 1:
            one_snp = True

        allele_phases = list(self._windows.keys())
        # reads_to_be_filtered_list = list()
        reads_with_alleles = df_w_returned[df_w_returned[SNP_PHASE].isin(allele_phases)]
        reads_with_no_alleles = df_w_returned[~df_w_returned[SNP_PHASE].isin(allele_phases)]
        reads_with_alleles = reads_with_alleles.assign(is_random=False)
        reads_with_no_alleles = reads_with_no_alleles.assign(is_random=False)
        if len(df_dropped) > 0:
            reads_with_no_alleles = pd.concat([reads_with_no_alleles, df_dropped])
        # if only one snp, then it has to be random
        if one_snp:  # TBD: check if needed. Maybe redundant (966). It's NOT!!!
            for i, row in reads_with_no_alleles.iterrows():
                reads_with_no_alleles.at[i, IS_RANDOM] = True
        else:
            windows_list = list(self._windows.values())
            for i, row in reads_with_no_alleles.iterrows():
                seq = row[ALIGNMENT_W_DEL]
                # find the best allele to assign this read. return best index of allele list
                best_index, score, is_random, _ = self._best_index_wrapper(seq, windows_list)
                # if score is None or it is random:
                if (score is None) or is_random:
                    reads_with_no_alleles.at[i, IS_RANDOM] = True
                else:
                    # if there is better alignment between the alleles - set it
                    snp_type = allele_phases[best_index]
                    reads_with_no_alleles.at[i, SNP_PHASE] = snp_type

        reads_with_alleles = pd.concat([reads_with_alleles, reads_with_no_alleles])

        return reads_with_alleles

    def _best_index_wrapper(self, curr_amplicon, SNP_window_information):
        """
        wrapper of the two functions that determine the best_index
        :param: curr_amplicon: The amplicon to be aligned to
        :param: SNP_window_information: list of SNP windows with additional information
        :return: df: the best index of allele to be set, alignment score, if it was picked by random, the whole list
        """

        # add to list the True/False for each CTC window
        CTC_list = list()  # TBD: put CTC into a function outside here and outside every row iteration
        for window_info in SNP_window_information[0]:
            CTC_list.append(window_info[1])
        scores = self._calc_alignments_scores(curr_amplicon, SNP_window_information, CTC_list)
        best_index, alignment_score, is_random, index_list = self._determine_best_idx(
            curr_amplicon, SNP_window_information, CTC_list, scores)
        return best_index, alignment_score, is_random, index_list

    def _calc_alignments_scores(self, curr_amplicon, SNP_window_information, CTC_list):
        """
        Calculates the best alignment score between each window to a given amplicon
        :param: curr_amplicon: The amplicon to be aligned to
        :param: SNP_window_information: list of SNP windows with additional information
        :param: CTC_list: list of True/False regarding the snp close to site or not
        :param: aligner: aligner to align according to
        :return: total_scores: list of lists - best score for each snp for each allele option.
                               Inner list - scores of all alleles in a certain snp
        """
        total_scores = list()
        total_i = list()

        # iterate over different snp
        for i, snp_option in enumerate(SNP_window_information[0]):
            temp_scores = list()
            smallest_i_list = list()
            relevant_read = curr_amplicon
            # iterate over all alleles
            for j in range(len(SNP_window_information)):

                SNP_window = SNP_window_information[j][i][0]
                CTC = CTC_list[i]

                smallest_gap = np.inf
                i_smallest = None

                # if the snp not close to the cut site - k-mars to find the best part of the read for alignment
                if not CTC:
                    if i == 0:
                        t = 0
                    else:
                        t = self._snp_locus[i-1] - len(SNP_window)
                    while t < len(relevant_read) - (len(SNP_window)) + 1:
                        window = relevant_read[t: t + len(SNP_window)]
                        curr_gap_score = count_diff(SNP_window, window)
                        if curr_gap_score < smallest_gap:
                            smallest_gap = curr_gap_score
                            i_smallest = t
                        t += 1

                    temp_scores.append(smallest_gap)
                    smallest_i_list.append(i_smallest)

                # if the snp is close to the cut site - the read might be very different then the window -
                # thus norm local alignment
                elif CTC:
                    alignment_for_norm = self._ll_aligner._align_seq_to_read(SNP_window, SNP_window)
                    max_CTC_score = alignment_for_norm.score
                    alignment = self._ll_aligner._align_seq_to_read(relevant_read, SNP_window)
                    alignment_score = alignment.score
                    norm_score = - alignment_score / max_CTC_score
                    temp_scores.append(norm_score)

            total_scores.append(temp_scores)
            total_i.append(smallest_i_list)

        return total_scores

    def _determine_best_idx(self, curr_amplicon, SNP_window_information, CTC_list, scores):
        """
        Determine the index with the best score out of list of lists
        :param: scores: list of lists with scores as values
        :param: CTC_list: list of Close To Cut with relative to each snp
        :param: curr_amplicon: The amplicon to be aligned to
        :param: SNP_window_information: list of SNP windows with additional information
        :return: returns the best index (or False) and the index_list
        """
        index_list = list()
        CTC_alignment = False
        is_random = False
        # iterate over all list. each list is the scores of all alleles in a given snp
        for score in scores:
            # set the minimum score in the list and find all indexes that comply the minimum score
            min_score = min(score)
            all_min_idx = [idx for idx, val in enumerate(score) if val == min_score]
            index_list.append(all_min_idx)
        # get the best indexes across all lists
        best_indexes = list(set(index_list[0]).intersection(*index_list))
        # if the best_indexes contains only one index - return it
        if len(best_indexes) == 1:
            best_index = best_indexes[0]
            alignment_score = np.mean(scores, axis=0)[best_index]
        # else - if there are several "best index":
        else:
            # if it is only one snp
            if len(scores) == 1:  # Not redundant to line (582). important for Tx analysis
                is_random = True
                best_index = random.choice(index_list[0])
                alignment_score = np.mean(scores, axis=0)[best_index]
                if CTC_list[0]:
                    CTC_alignment = True
            else:
                # flatten all lists to one list of best indexes
                flat_best_index_list = [item for sublist in index_list for item in sublist]
                # get the indexes that repeated the most out of the flatted list
                new_best_indexes = get_most_freq_val(flat_best_index_list)
                # if there is an allele that is more dominant in the sense of fitting to window - get it.
                if len(new_best_indexes) == 1:
                    best_index = new_best_indexes[0]
                    alignment_score = np.mean(scores, axis=0)[best_index]
                else:
                    # create 2 lists - one to contain all scores, and
                    # one to contain all scores but without scores of snps that close to cut-site
                    sub_scores_w_CTC = list()
                    sub_scores_wo_CTC = list()
                    # iterate over all scores
                    for j, score in enumerate(scores):
                        # append to the lists the scores in the new_best_indexes
                        sub_scores_w_CTC.append([x for i, x in enumerate(score) if i in new_best_indexes])
                        if not CTC_list[j]:
                            sub_scores_wo_CTC.append([x for i, x in enumerate(score) if i in new_best_indexes])
                    sums = np.sum(sub_scores_wo_CTC, axis=0)
                    # if all sums are equal - we need to try different approach
                    if all(CTC_list):
                        CTC_alignment = True
                        ctc_sums = np.sum(sub_scores_w_CTC, axis=0)
                        # if all sums are equal - choose randomly out of the option of new_best_indexes
                        if np.count_nonzero(np.array(ctc_sums) == min(ctc_sums)) > 1:
                            is_random = True
                            best_index = random.choice(new_best_indexes)
                            alignment_score = np.mean(scores, axis=0)[new_best_indexes.index(best_index)]
                        # else - we have a better index in a sense of scoring, thus we will choose it.
                        else:
                            best_index = new_best_indexes[np.argmin(ctc_sums)]
                            alignment_score = np.mean(scores, axis=0)[new_best_indexes.index(best_index)]
                    elif np.count_nonzero(np.array(sums) == min(sums)) > 1:
                        # if there is a CTC between snps
                        if (sub_scores_w_CTC != sub_scores_wo_CTC) and (len(sub_scores_wo_CTC) > 0):
                            CTC_alignment = True
                            # consider all snps CTC and calculate from the beginning the scores
                            new_CTC_list = [True] * len(CTC_list)
                            scores_new = self._calc_alignments_scores(curr_amplicon, SNP_window_information,
                                                                      new_CTC_list)
                            scores = np.array(scores_new)[:, new_best_indexes]
                            new_sums = np.sum(scores, axis=0)
                            # if all sums are equal - choose randomly out of the option of new_best_indexes
                            if np.count_nonzero(np.array(new_sums) == min(new_sums)) > 1:
                                is_random = True
                                best_index = random.choice(new_best_indexes)
                                alignment_score = np.mean(scores, axis=0)[new_best_indexes.index(best_index)]
                            # else - we have a better index in a sense of scoring, thus we will choose it.
                            else:
                                best_index = new_best_indexes[np.argmin(new_sums)]
                                alignment_score = np.mean(scores, axis=0)[new_best_indexes.index(best_index)]
                        # else - there is no CTC among snps,
                        # thus we don't have more calculation to make, so we will pick index randomly
                        else:
                            is_random = True
                            new_min_idx = [idx for idx, val in enumerate(sums) if val == min(sums)]
                            best_index = new_best_indexes[random.choice(new_min_idx)]
                            alignment_score = np.mean(scores, axis=0)[new_best_indexes.index(best_index)]
                    # else - we have a better index in a sense of scoring, thus we will choose it.
                    else:
                        best_index = new_best_indexes[np.argmin(sums)]
                        alignment_score = np.mean(scores, axis=0)[new_best_indexes.index(best_index)]

        if CTC_alignment:
            MAX_CTC_MM = - self._ll_aligner._aligner.match_score * (self._max_len_snv_ctc * 2 + 1 - self._max_mm) \
                         / (self._ll_aligner._aligner.match_score * (self._max_len_snv_ctc * 2 + 1))
            if alignment_score > MAX_CTC_MM:
                alignment_score = None
        else:
            if alignment_score > self._max_mm:
                alignment_score = None

        return best_index, alignment_score, is_random, index_list

# -------------------------------#
######### Global Functions #######
# -------------------------------#


def count_diff(seq_A, seq_B):
    """
    Count the differences between 2 sequences
    :param: seq_A: first sequence
    :param: seq_B: second sequence
    :return: counter: the number of mismatches at the same positions
    """
    counter = 0
    for i, j in zip(seq_A, seq_B):
        if i != j:
            counter += 1
    return counter


def calc_entropy(temp_df):
    """
     Calculates entropy of the passed `pd.Series`
     :param: df: One nucleotide base through all reads in one mock site, along with frequency of each read
     :return: entropy: The entropy for this particular nucleotide base;
              sorted_nuc_dict: A dictionary of the nucleotides' distribution
     """
    # computing unique nuc bases
    df = temp_df.set_axis([BASE, FREQ], axis=1)
    nuc_bases = set(df.iloc[:, 0])
    nuc_dict = {}
    # computing nuc bases distribution (including deletions '-')
    for nuc in nuc_bases:
        nuc_dict[nuc] = df.loc[df[BASE] == nuc, FREQ].sum()
    sorted_nuc_dict = dict(sorted(nuc_dict.items(), key=lambda item: item[1], reverse=True))
    # consider changing the base according to expected alleles, variable number_of_relevant_alleles
    # number_of_relevant_alleles = len(self._ratios)
    dict_for_entropy = {k: sorted_nuc_dict[k] for k in list(sorted_nuc_dict)[:2]}
    # calc entropy
    _entropy = scipy.stats.entropy(pd.Series(dict_for_entropy), base=2)

    return _entropy, sorted_nuc_dict


def compute_lengths_distribution(sorted_df, possible_len):
    """
     Compute the distribution of the length of all the *aligned* reads in the df
     :param sorted_df: sorted df by frequency
     :param possible_len: list of all possible lengths in the df
     :return: lengths_ratio: The ratio of the largest length percentage;
              len_dict: Dictionary of lengths percentages
     """
    len_dict = {}
    for i_len in possible_len:
        len_dict[i_len] = sorted_df.loc[sorted_df[LEN] == i_len, FREQ].sum()
    sum_len = sum(len_dict.values())
    values_len = list(len_dict.values())
    percent_len = [x / sum_len for x in values_len]
    lengths_ratio = np.sort(percent_len)[::-1]

    return lengths_ratio, len_dict


def get_most_frq_cut_site(same_len_df):
    """
     Retrieve the most frequent cut-site to be the cut-site of the site
     :param: same_len_df: df of the site
     :return: the cut-site
     """
    cut_site_dict = dict()
    freq = list(same_len_df[FREQ])
    cut_sites = list(same_len_df[ALIGN_CUT_SITE])
    for i, cut in enumerate(cut_sites):
        if cut in cut_site_dict.keys():
            cut_site_dict[cut] += freq[i]
        else:
            cut_site_dict[cut] = freq[i]

    return int(max(cut_site_dict, key=cut_site_dict.get))


def score_normalization(_entropy, _coverage):
    # TBD: optimize the ratios between coverage and entropy
    """
     Find all snps and return df with the snp bases
     :param: _entropy: The entropy of the current number of alleles
     :param: _coverage: The current coverage of these number of alleles out of 100%
     :return: Normalized calculation between coverage and entropy
     """
    return math.pow(_coverage, SCORE_NORM_RATE) * math.pow(_entropy, 1 - SCORE_NORM_RATE)


def get_most_freq_val(lst):
    """
    Get the most frequent values from a list
    :param: lst: list of values
    :return: the most frequent values (can be many)
    """
    # # using df to get the desired output
    # df = pd.DataFrame({'Number': lst})
    # df1 = pd.DataFrame(data=df['Number'].value_counts())
    # df1['Count'] = df1['Number']
    # df1['Number'] = df1.index
    # df1 = df1.reset_index(drop=True)
    #
    # return list(df1[df1['Count'] == df1.Count.max()]['Number'])

    counts = Counter(lst)
    max_count = counts.most_common(1)[0][1]
    out = [value for value, count in counts.most_common() if count == max_count]

    return out




# TBD: consider delete this function as it's not having any functionality
def _prepare_windows_for_alignment(windows_list_for_realignment, allele_phases):
    windows_list = list()
    for phase in allele_phases:
        temp_allele_list = list()
        for i, window in enumerate(windows_list_for_realignment):
            new_window = window[0]
            for j in range(len(phase)):
                n_index = new_window.find(N)
                new_window = new_window[:n_index] + phase[j] + new_window[n_index + 1:]
                # new_window = window[0].replace(N, phase[i])
            temp_allele_list.append((new_window, window[1], window[2], window[3]))
        windows_list.append(temp_allele_list)
    return windows_list

