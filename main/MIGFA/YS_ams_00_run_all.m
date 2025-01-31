% need to modify MIGFA output into NFGDA
% for example in this case
% KABR20140621_215954_V06
% KABR20140621_220428_V06
% KABR20140621_220901_V06
% KABR20140621_221335_V06
% this will be converted to the folder ../MIGFA as
% KABR14062101
% KABR14062102
% KABR14062103
% KABR14062104
% because that is how NFGDA reads MIGFA info
% take a look at the folder ../MIGFA/20140621-KABR
% you need to do this by hand but you can make codes 
% to do that, I just did not have much time, though
% % % % % % % % % % % % % % % % % 
% take a look at ../MIGFA/folder
% % % % % % % % % % % % % % % % % 
% in case of using MIGFA as comparison
% Don't forget Mlat and Mlon in NF00_header
% which determines the location of radar
% 

YS_ams_00_cord_read_data
YS_ams_01_cord_convert
YS_ams_02_skel_fit