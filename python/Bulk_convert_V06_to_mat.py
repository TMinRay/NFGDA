from NF01_convert_V06_to_mat import *
clist=[
# 'KABX20200704_02',
# 'KABX20200705_21',
'KABX20200707_01',
'KABX20200712_21',
'KABX20200715_23',
'KABX20200721_03',
'KABX20200721_19',
'KABX20200724_21',
'KABX20200726_19',
'KABX20210702_21',
'KABX20210704_00',
'KABX20210705_05',
'KABX20210706_00',
'KABX20210706_23',
'KABX20210707_01',
'KABX20210709_22']
for case_name in clist:
	convert_v06_to_mat(v06_folder="../V06", case_id=case_name, mat_folder="../mat",
	                   i_start=0, i_end=99)