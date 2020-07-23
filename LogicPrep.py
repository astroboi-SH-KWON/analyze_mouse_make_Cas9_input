from Bio import SeqIO
import re

import Logic
import Util

class LogicPreps:

    def __init__(self):
        self.ext_fa = ".fa"
        self.ext_dat = ".dat"
        self.ext_gtf = ".gtf"

    def check_strnd(self, loc_str):
        if "-" in loc_str:
            return "-"
        return "+"

    def get_target_seq_with_clvg_site_fr_fasta(self, path, init_arr):
        logic = Logic.Logics()
        tmp_dict = {}

        pam_seq = init_arr[0]
        add_seq1_len = init_arr[1]
        spacer_len = init_arr[2]
        add_seq2_len = init_arr[3]
        clvg_site = init_arr[4]
        pam_len = len(pam_seq)

        std_tot_len = add_seq1_len + spacer_len + pam_len + add_seq2_len

        for seq_record in SeqIO.parse(path, "fasta"):
            tot_cds_len = len(seq_record.seq)

            trncrpt_id = seq_record.id
            if trncrpt_id not in tmp_dict:
                tmp_dict.update({trncrpt_id: [seq_record.description]})

            tmp_p_str = ""
            idx = 0
            for c in seq_record.seq:
                idx += 1
                tmp_p_str = tmp_p_str + c.upper()

                if len(tmp_p_str) > std_tot_len:
                    tmp_p_str = tmp_p_str[-std_tot_len:]

                if len(tmp_p_str) == std_tot_len:
                    if 'N' not in tmp_p_str:
                        if logic.match(0, tmp_p_str[-(add_seq2_len + pam_len):-add_seq2_len], pam_seq):
                            tmp_dict[trncrpt_id].append(
                                [tmp_p_str, ((idx - clvg_site - pam_len - add_seq2_len) / tot_cds_len) * 100])

        return tmp_dict

    def target_seq_with_clvg_site_group_by_chromosome(self, trgt_seq_dict, deli_str=":"):
        result_dict = {}
        for trnscrpt_id, vals_arr in trgt_seq_dict.items():
            dscript = vals_arr[0]
            chrsm = "chr"
            try:
                chrsm += dscript.split(deli_str)[2]
            except ValueError as err:
                print(dscript + " : ", err)

            if chrsm in result_dict:
                if trnscrpt_id not in result_dict[chrsm]:
                    result_dict[chrsm].update({trnscrpt_id: vals_arr})
            else:
                result_dict.update({chrsm: {trnscrpt_id: vals_arr}})

        return result_dict

    def get_target_seq_with_clvg_site(self, path, init_arr):
        logic = Logic.Logics()
        tmp_dict = {}

        pam_seq = init_arr[0]
        add_seq1_len = init_arr[1]
        spacer_len = init_arr[2]
        add_seq2_len = init_arr[3]
        clvg_site = init_arr[4]
        pam_len = len(pam_seq)

        std_tot_len = add_seq1_len + spacer_len + pam_len + add_seq2_len

        for seq_record in SeqIO.parse(path, "fasta"):
            tot_cds_len = len(seq_record.seq)

            trncrpt_id = seq_record.id
            if trncrpt_id not in tmp_dict:
                tmp_dict.update({trncrpt_id: [seq_record.description]})

            tmp_p_str = ""
            idx = 0
            for c in seq_record.seq:
                idx += 1
                tmp_p_str = tmp_p_str + c.upper()

                if len(tmp_p_str) > std_tot_len:
                    tmp_p_str = tmp_p_str[-std_tot_len:]

                if len(tmp_p_str) == std_tot_len:
                    if 'N' not in tmp_p_str:
                        if logic.match(0, tmp_p_str[-(add_seq2_len + pam_len):-add_seq2_len], pam_seq):
                            tmp_dict[trncrpt_id].append([tmp_p_str, ((idx - clvg_site - pam_len - add_seq2_len) / tot_cds_len) * 100])

        return tmp_dict

    def get_cs9_scre(self, scre_txt_path):
        tmp_tpl = ()
        with open(scre_txt_path, "r") as f:
            f.readline().replace("\n", "")
            f.readline().replace("\n", "")
            f.readline().replace("\n", "")
            f.readline().replace("\n", "")
            # make result as tuple
            tmp_tpl += eval(f.readline().replace("\n", ""))
        return tmp_tpl

    def get_deep_cas9_tupl(self, path, scre_txt_path, seq_txt_path):
        tmp_dict = {}
        tmp_tpl = self.get_cs9_scre(path + scre_txt_path)

        with open(path + seq_txt_path) as f:
            f.readline().replace("\n", "")
            idx = 0
            while True:
                tmp_line = f.readline().replace("\n", "")
                if tmp_line == "":
                    break

                tmp_arr = tmp_line.split("\t")
                tmp_dict.update({tmp_arr[0]: tmp_tpl[idx]})
                idx += 1

        return tmp_dict

    def merge_cas9_abe_cbe_to_list(self, chr_key, init_dict_arr, result_list):
        trnscrpt_val = init_dict_arr[0]
        abe_score_dict = init_dict_arr[1]
        cbe_score_dict = init_dict_arr[2]
        cs9_score_dict = init_dict_arr[3]

        check_duple = set()
        for trnscrpt_id, vals_arr in trnscrpt_val.items():
            full_description = vals_arr[0]
            full_description_arr = full_description.split(" ")
            gene_id = full_description_arr[3].replace("gene:", "")
            strand = "+"
            if ":-1" in full_description_arr[2]:
                strand = "-"

            gene_nm = ""
            description = ""
            if "gene_symbol:" in full_description:
                if "description:" in full_description:
                    gene_nm = full_description[
                              full_description.index("gene_symbol:") + len("gene_symbol:"):full_description.index(
                                  "description:")]
                    description = full_description[full_description.index("description:") + len("description:"):]
                else:
                    gene_nm = full_description[
                              full_description.index("gene_symbol:") + len("gene_symbol:"):]

            for trgt_idx in range(1, len(vals_arr)):
                cntxt_seq = vals_arr[trgt_idx][0]

                # check same seq in check_duple
                if cntxt_seq in check_duple:
                    continue
                else:
                    check_duple.add(cntxt_seq)

                clvg_site = vals_arr[trgt_idx][1]
                cas9_score = 0
                abe_score = 0
                cbe_score = 0

                if cntxt_seq in cs9_score_dict:
                    cas9_score = cs9_score_dict[cntxt_seq]
                else:
                    print(cntxt_seq + " doesn't have CAS9 score")
                if cntxt_seq in abe_score_dict:
                    abe_score = abe_score_dict[cntxt_seq]
                else:
                    pass
                    # print(cntxt_seq + " doesn't have ABE score")
                if cntxt_seq in cbe_score_dict:
                    cbe_score = cbe_score_dict[cntxt_seq]
                else:
                    pass
                    # print(cntxt_seq + " doesn't have CBE score")

                result_list.append(
                    [chr_key, gene_nm, description, trnscrpt_id, gene_id, strand, cntxt_seq, clvg_site, float(cas9_score),
                     float(abe_score), float(cbe_score)])

        return result_list

    def sort_by_idx_element(self, trgt_list, idx, result_list):
        for tmp_list in sorted(trgt_list, key=lambda tmp_list: tmp_list[idx], reverse=True):
            result_list.append(tmp_list)
        return result_list


