import time
import os

import Util
import Logic
import LogicPrep
############### start to set env ################
# WORK_DIR = "D:/000_WORK/JangHyeWon_ShinJeongHong/20200527/WORK_DIR/"
WORK_DIR = os.getcwd() + "/"
PROJECT_NAME = WORK_DIR.split("/")[-2]

REF_PATH = "D:/000_WORK/000_reference_path/mouse/"
CDS_FILE = "cds/Mus_musculus.GRCm38.cds.all.fa"
DNA_FILE = "dna/"
ANNO_FILE = "genbank_anno/"


FRONT_WIN_LEN = 4
gRNA_LEN = 20
PAM_SEQ = "NGG"
BACK_WIN_LEN = 20
BATCH_SIZE = 2000000

# INIT_DEEP_PE = [PAM_SEQ, FRONT_WIN_LEN, gRNA_LEN, BACK_WIN_LEN]
A_or_C_IDX = [4, 10]
ACTG_RULE = ['A', 'C']
############## make_deep_pe_input ##############
BE_BACK_WIN_LEN = 20
CLEAVAGE_SITE = 3
MAX_MISMATCH = 3
REF_SRV_PATH = "FASTA/marmoset"
INIT_BE = [PAM_SEQ, FRONT_WIN_LEN, gRNA_LEN, BE_BACK_WIN_LEN, CLEAVAGE_SITE]
#################### top N #####################
TOP_N = 100
# TOP_N = 2000000
INIT_MERGE_BY_CHAR = [REF_PATH, CDS_FILE, A_or_C_IDX, ACTG_RULE, WORK_DIR, TOP_N]
############### end setting env ################

def sort_n_merge_by_chr():
    logic = Logic.Logics()
    logic.sort_n_merge_by_chr(INIT_MERGE_BY_CHAR, INIT_BE)

def make_deep_cas9_input():
    logic = Logic.Logics()
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    trgt_seq_dict = logic_prep.get_target_seq_with_clvg_site_fr_fasta(REF_PATH + CDS_FILE, INIT_BE)
    chr_dict = logic_prep.target_seq_with_clvg_site_group_by_chromosome(trgt_seq_dict)

    util.make_deep_cas9_input(WORK_DIR + "deep_cas_9/sample", [chr_dict], INIT_BE, BATCH_SIZE)



if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [" + PROJECT_NAME + "]>>>>>>>>>>>>>>>>>>")
    sort_n_merge_by_chr()  # 2
    # make_deep_cas9_input()  # 1
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))