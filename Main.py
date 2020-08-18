import time
import os

import Util
import Logic
import LogicPrep
############### start to set env ################
# WORK_DIR = "D:/000_WORK/JangHyeWon_ShinJeongHong/20200527/WORK_DIR/"
WORK_DIR = os.getcwd() + "/"
PROJECT_NAME = WORK_DIR.split("/")[-2]

# REF_PATH = "D:/000_WORK/000_reference_path/mouse/"
REF_PATH = "D:/000_WORK/000_reference_path/monkey/marmoset/"
CDS_FILE = "cds/Callithrix_jacchus.ASM275486v1.cds.all.fa"
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
INIT_BE = [PAM_SEQ, FRONT_WIN_LEN, gRNA_LEN, BE_BACK_WIN_LEN, CLEAVAGE_SITE, REF_PATH.split("/")[-2]]
#################### top N #####################
TOP_N = 1000
# TOP_N = 2000000
INIT_MERGE_BY_CHAR = [REF_PATH, CDS_FILE, A_or_C_IDX, ACTG_RULE, WORK_DIR, TOP_N]
############### end setting env ################

def sort_n_merge_by_chr():
    logic = Logic.Logics()
    # logic.sort_n_merge_by_chr(INIT_MERGE_BY_CHAR, INIT_BE)  # plan A
    logic.sort_n_merge_by_chr_one_file(INIT_MERGE_BY_CHAR, INIT_BE)  # plan B # TOP_N <= 100

def make_deep_cas9_input():
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    trgt_seq_dict = logic_prep.get_target_seq_with_clvg_site_fr_fasta(REF_PATH + CDS_FILE, INIT_BE)
    chr_dict = logic_prep.target_seq_with_clvg_site_group_by_chromosome(trgt_seq_dict)

    util.make_deep_cas9_input(WORK_DIR + "deep_cas_9/sample", [chr_dict], INIT_BE, BATCH_SIZE)



if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [" + PROJECT_NAME + "]>>>>>>>>>>>>>>>>>>")
    # make_deep_cas9_input()  # 1
    """
    make_deep_cas9_input() 돌리기전 deep_cas_9 폴더에 sample_0.txt, sample_1.txt 삭제
    # 1의 결과물 sample_0.txt, sample_1.txt을 각각 순차적으로
    D:\Github\monkey_PE2_efficiency_prediction_merge\DeepCas9\에
    DeepCas9_example_input.txt 으로 변경(header는 그대로, seq만 사용)

    Test.py 돌려서 RANK_final_DeepCas9.txt 결과 파일 생성
    """
    sort_n_merge_by_chr()  # 2
    """
    RANK_final_DeepCas9.txt 결과 파일들과 sample_0.txt, sample_1.txt를 번호 매칭
    ex) RANK_final_DeepCas9_0.txt 와 RANK_final_DeepCas9_1.txt
    chromosome별로 TOP_N개씩 모아서 한 파일로 만들기 : sort_n_merge_by_chr_one_file
    """

    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))