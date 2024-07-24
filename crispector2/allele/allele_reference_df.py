from crispector2.utils.constants_and_types import SITE_NAME, CUT_SITE, PAM_WINDOW, GRNA_WINDOW, SGRNA_REVERSED, \
    REFERENCE
import pandas as pd
import ast


class ref_dfAlleleHandler:
    """
    Class to handle the creation of new ref_df with allele sites
    """

    def __init__(self):
        self._ref_df = None

    def run(self, ref_df, new_allele_sites):
        self._ref_df = ref_df
        for original_site_name, sites_list in new_allele_sites.items():
            for allele, allele_info in sites_list.items():
                new_site_name = allele_info[0]
                new_amplicon = allele_info[4]

                new_row_for_ref_df = self._ref_df[self._ref_df[SITE_NAME] == original_site_name]
                new_row_for_ref_df = new_row_for_ref_df.rename(index={original_site_name: new_site_name})
                new_row_for_ref_df.loc[new_site_name, REFERENCE] = new_amplicon
                new_row_for_ref_df.loc[new_site_name, SITE_NAME] = new_site_name
                self._ref_df = pd.concat([self._ref_df,
                                          new_row_for_ref_df])
                # TBD: maybe replace concat with append to list like I did in previous cases
        '''add PAM coordinates and gRNA coordinates to df'''
        # TBD: WHY to get the PAM?? cannot remember
        self._get_PAM_site()

        return self._ref_df

    def _get_PAM_site(self):
        """
        Set the gRNA coordinate and PAM coordinate for each site
        :param:
        :return: assign values directly to df
        """
        self._ref_df[PAM_WINDOW] = None
        self._ref_df[GRNA_WINDOW] = None

        for i, row in self._ref_df.iterrows():
            cut_site = row[CUT_SITE]
            if row[SGRNA_REVERSED]:
                PAM = [cut_site - 6, cut_site - 4]
                grna = [cut_site - 3, cut_site + 16]
            elif not row[SGRNA_REVERSED]:
                PAM = [cut_site + 3, cut_site + 5]
                grna = [cut_site - 17, cut_site + 2]

            self._ref_df.at[i, PAM_WINDOW] = PAM
            self._ref_df.at[i, GRNA_WINDOW] = grna


def is_snp_in_pam_grna(ref_df):
    """
    If the snps that were found are part of the PAM or gRNA, raise a WARNING
    :param ref_df:
    :return: None if all snps not in grna or pam
    """
    sites_done = list()
    for site_name, site_info in ref_df.iterrows():
        if '[' in site_name:
            temp_general_site_name = site_name[:site_name.index('[')]
            if temp_general_site_name not in sites_done:
                PAM_window = site_info[PAM_WINDOW]
                grna_window = site_info[GRNA_WINDOW]
                snps_locus = site_name[site_name.index('['): site_name.index(']')+1]
                snps_locus = ast.literal_eval(snps_locus)

                for snp in snps_locus:
                    if PAM_window[0] <= snp <= PAM_window[1]:
                        return f'one of the SNPs in site {temp_general_site_name[:-1]} is located inside the PAM! ' \
                               'Make sure to examine the results carefully'
                    elif grna_window[0] <= snp <= grna_window[1]:
                        return f'one of the SNPs in site {temp_general_site_name[:-1]} is located inside the gRNA! ' \
                               'Make sure to examine the results carefully'

                sites_done.append(temp_general_site_name)

    return None

