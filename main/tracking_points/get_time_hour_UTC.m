function time_hr = get_time_hour_UTC(ppi_name)
hh = ppi_name(14:15);
mm = ppi_name(16:17);
ss = ppi_name(18:19);

% disp(ppi_num);
% disp(ppi_name);
% fprintf("%2s:%2s:%2s\n",hh,mm,ss);
time_hr = str2num(hh) + str2num(mm)/60 + str2num(ss)/3600;
end