% % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % NFGDA MATLAB CODES 
% % % % % % Yunsung Hwang
% % % % % % contact: hisnameys@gmail.com
% % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % see the subfiles below
% % % % % % Tested with KABR 2014.06.21 data
% % % % % % Take a look at ../V06/ folder
% % % % % % Take a look at ../MIGFAdat/ folder
% % % % % % Take a look at ../IMG/ --> FTC
% % % % % % Take a look at ../netcdf/ --> ncfiles
% % % % % % Take a look at ../mat/ --> matfiles

% % % NF01_convert_nc_to_mat

% NF01_convert_V06_to_mat

NF01_convert_to_cartesian
NF02_calc_5by5_SD
NF04_calc_DeltaZ
NF05_calc_BETA_LINEFEATURE
NF06_calc_6variables_preprocessing

NF07_handpick_region
NF07_obtaining_evaluation_box

% % -------------pre-processing
% % + obtaining evaluation box
NF08_making_training_dataset
% % + obtaining trainingdataset
% % -- defined in NF00_header
% % use command "anfisedit" to train NFsystem
tic
NF10_evalfuzzy_and_skel
toc
% % 1. obtain fuzzy output
% % 2. 3 by 3 window filtering
% % 3. skeleton thining
NF11_postprocessing_movingavg
% % postprocess --> obtain probable area of GF
NF12_final_output_skel
% % final output-thinned line (detected GF)
NF13_making_stats_evaluation
NF14_making_stats_scores_final

% --->NFGDA scores

% % % % % MIGFA
% run ./MIGFA/YS_ams_00_run_all.m;
% % % % % MIGFA
% NF16_MIGFA_stats_evaluation
% NF17_MIGFA_stats_scores_migfa_final
% % --->MIGFA scores


% % % % % Stats
% % % % % Stats
% % % % % Stats
NF99_STAT_3diff_regions_cases_histogram
NF99_STAT_GF_nonGF_histogram

% % % % FIGS --> ./NFFIG/
NFFIG_00_STAT_3region
NFFIG_01evlabox
NFFIG_02NFinputvari_plot_from_inputNF
NFFIG_03FINAL_output
NFFIG_04Membership_functions







