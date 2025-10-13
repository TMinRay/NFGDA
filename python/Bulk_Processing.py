import sys
clist=[
'KABX20200704_02',
'KABX20200705_21',
]


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python Bulk_Processing.py [convert|detect|plot]")
        sys.exit(1)

    mode = sys.argv[1].lower()

    if mode == "convert":
        from NF01_convert_V06_to_mat import *
        for case_name in clist:
            convert_v06_to_mat(v06_folder="../V06", case_id=case_name, mat_folder="../mat",
                               i_start=0, i_end=99)
    
    elif mode == "detect":
        from NFGDA import *
        for case_name in clist:
            nfgda_proc(case_name)
    
    elif mode == "plot":
        from NFFig import *
        for case_name in clist:
            nffig_proc(case_name)
    else:
        print(f"Unknown mode: {mode}")
        print("Valid options are: convert|detect|plot")
        sys.exit(1)