from Bio import SeqIO

import Util
import LogicPrep

class Logics:
    def __init__(self):
        pass

    """
    checkSeqByChar : match sequences by char with rules
    :param
        dna_char :
        rule_char : rules with "A", "C", "G", "T", "U", "N", "R",...
    :return
        boolean
    """
    def checkSeqByChar(self, dna_char, rule_char):
        flag = False
        if rule_char == 'N':
            return True
        elif rule_char in 'ACGTU':
            if dna_char == rule_char:
                return True
        elif rule_char == 'R':
            if dna_char in 'AG':
                return True
        # elif rule_char == 'r':
        #     if dna_char in 'CT':
        #         return True
        """
        add more rules of "ACTGU"
        """
        return flag

    """
    match : match sequence with same length strings
    :param
        i : index of seq
        dna_seq : targeted DNA/RNA sequence 
        rule_str : rules with "ACGTU", "N", "R",...
    :return
        boolean
    """

    def match(self, i, dna_seq, rule_str):
        if len(dna_seq) == i:
            return True
        if self.checkSeqByChar(dna_seq[i], rule_str[i]):
            return self.match(i + 1, dna_seq, rule_str)
        else:
            return False

    def sort_n_merge_by_chr(self, init_merge, init_be):
        ref_path = init_merge[0]
        cdf_file = init_merge[1]
        a_or_c_idx = init_merge[2]
        a_c_rule = init_merge[3]
        work_dir = init_merge[4]
        top_n = init_merge[5]

        logic_prep = LogicPrep.LogicPreps()
        util = Util.Utils()

        trgt_seq_dict = logic_prep.get_target_seq_with_clvg_site(ref_path + cdf_file, init_be)
        chr_dict = logic_prep.target_seq_with_clvg_site_group_by_chromosome(trgt_seq_dict)

        cs9_score_dict = {}
        cs9_score_dict.update(logic_prep.get_deep_cas9_tupl(work_dir + "deep_cas_9/", "RANK_final_DeepCas9_0.txt",
                                                            "sample_0.txt"))
        cs9_score_dict.update(logic_prep.get_deep_cas9_tupl(work_dir + "deep_cas_9/", "RANK_final_DeepCas9_1.txt",
                                                            "sample_1.txt"))

        for chr_key, trnscrpt_list in chr_dict.items():
            result_list = []
            result_list = logic_prep.merge_cas9_abe_cbe_to_list(chr_key, [trnscrpt_list, {}, {},
                                                                          cs9_score_dict], result_list)

            sort_by_cas9_list = logic_prep.sort_by_idx_element(result_list, -3, [])

            """
            # extend TOP N lists to (top_n_abe_list, top_n_cbe_list)
            it needs filter out same context seq in different trnscrpt
            """
            util.make_tsv_after_sorting(work_dir + "output/mouse_seq_sorted_by_CAS9_" + chr_key, sort_by_cas9_list, init_be)

    def sort_n_merge_by_chr_one_file(self, init_merge, init_be):
        ref_path = init_merge[0]
        cdf_file = init_merge[1]
        a_or_c_idx = init_merge[2]
        a_c_rule = init_merge[3]
        work_dir = init_merge[4]
        top_n = init_merge[5]

        logic_prep = LogicPrep.LogicPreps()
        util = Util.Utils()

        trgt_seq_dict = logic_prep.get_target_seq_with_clvg_site(ref_path + cdf_file, init_be)
        chr_dict = logic_prep.target_seq_with_clvg_site_group_by_chromosome(trgt_seq_dict)

        cs9_score_dict = {}
        cs9_score_dict.update(logic_prep.get_deep_cas9_tupl(work_dir + "deep_cas_9/", "RANK_final_DeepCas9_0.txt",
                                                            "sample_0.txt"))
        cs9_score_dict.update(logic_prep.get_deep_cas9_tupl(work_dir + "deep_cas_9/", "RANK_final_DeepCas9_1.txt",
                                                            "sample_1.txt"))

        top_n_list = []
        for chr_key, trnscrpt_list in chr_dict.items():
            result_list = []
            result_list = logic_prep.merge_cas9_abe_cbe_to_list(chr_key, [trnscrpt_list, {}, {},
                                                                          cs9_score_dict], result_list)

            sort_by_cas9_list = logic_prep.sort_by_idx_element(result_list, -3, [])

            top_n_list.extend(sort_by_cas9_list[:top_n + 1])

        # make tsv file result
        util.make_tsv_after_sorting(work_dir + "output/mouse_seq_sorted_by_CAS9_top_" + str(top_n), top_n_list, init_be)
        # make excel result
        util.make_excel_after_sorting(work_dir + "output/mouse_seq_sorted_by_CAS9_top_" + str(top_n), top_n_list, init_be)


