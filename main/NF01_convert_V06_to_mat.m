NF00_header;

V06PATH = "../V06";
matPATH = "../mat";
python_exe = "C:\Users\preci\Anaconda3\envs\pyart_matlab\python";

for cindex=1:numel(ttable(:,1))

    PUTDAT=ttable(cindex,:);

    startm=startt(cindex);
    endm=endt(cindex);

    curr_command = join([python_exe, '..\python_accessories\main\NF01_convert_V06_to_mat.py', ...
        '-v06_folder', V06PATH, ...
        '-case_id', PUTDAT, ...
        '-mat_folder', matPATH, ...
        '-i_start', (startm-1), ...
        '-i_end',  (endm-1)]);
    
    disp(curr_command);
    system(curr_command);

end


