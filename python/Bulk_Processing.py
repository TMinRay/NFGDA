import sys
clist=[
# 'KABX20200704_02',
# 'KABX20200705_21',
# 'KABX20200721_19',
'KABX20230703_01',
# 'KABX20230704_03',
'KAMA20240227_21',

# 'KABR20140621_21',
# 'KABX20200704_02',
# 'KABX20200705_21',
# 'KABX20200707_01',
# 'KABX20200712_21',
# 'KABX20200715_23',
# 'KABX20200721_03',
# 'KABX20200721_19',
# 'KABX20200724_21',
# 'KABX20200726_19',
# 'KABX20200726_20',
# 'KABX20210702_21',
# 'KABX20210704_00',
# 'KABX20210705_05',
# 'KABX20210706_00',
# 'KABX20210706_23',
# 'KABX20210707_00',
# 'KABX20210707_01',
# 'KABX20210708_00',
# 'KABX20210708_23',
# 'KABX20210709_22',
# 'KABX20210711_22',
# 'KABX20230703_01',
# 'KABX20230704_03',
# 'KABX20230708_21',
# 'KABX20230711_20',
# 'KABX20230716_01',
# 'KABX20230720_23',
# 'KABX20230725_22',
# 'KABX20230727_01',
# 'KABX20230728_21',
# 'KABX20230729_01',
# 'KABX20250712_20',
# 'KAMA20240227_21',
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

    elif mode == "forecast":
        # from NF_step_forecast import *
        from NF_stochastic_forecast import *
        for case_name in clist:
            nfgda_forecast(case_name)
    elif mode == "operation":
        from NF_forecast_operation import *
        for case_name in clist:
            nfgda_forecast(case_name)
    else:
        print(f"Unknown mode: {mode}")
        print("Valid options are: convert|detect|plot|forecast")
        sys.exit(1)