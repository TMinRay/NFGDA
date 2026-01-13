function [scoremt] = gen_tot_score(a2,c_para,s_para,thrREF,numINT,scorediv,step)
    % cys_num = [cnum1,cnum2,ynum1,snum1,snum2,synum1];
    % cnum1 = cys_num(1);
    % cnum2 = cys_num(2);
    % ynum1 = cys_num(3);
    % snum1 = cys_num(4);
    % snum2 = cys_num(5);
    % synum1 = cys_num(6);
    cnum1 = c_para(1);
    cnum2 = c_para(2);
    csig1 = c_para(3);
    cfactor1 = c_para(4);
    cintersec1 = c_para(5);
    csig2 = c_para(6);
    cfactor2 = c_para(7);
    cintersec2 = c_para(8);
    cyfill = c_para(9);

    datacy = [-8:8].';
    datacx = zeros(17,1);
    datasy = [-7:2:-1,0,1:2:7,-7:2:-1,0,1:2:7].';
    datasx = [-4*ones(9,1);4*ones(9,1)];
    indexcnum = numel(datacx);
    [row_indices, col_indices] = find(a2>thrREF);
    inINT = (row_indices>numINT) & (row_indices<=401-numINT)...
        &(col_indices>numINT) & (col_indices<=401-numINT);
    row_indices = row_indices(inINT);
    col_indices = col_indices(inINT);
    cridx = row_indices.' + datacy;
    ccidx = col_indices.' + datacx;
    c_indices = sub2ind(size(a2), cridx, ccidx);
    cbox=a2(c_indices);
    indcb = cbox>thrREF;
    numtcb=sum(indcb,1);
    cbr = numtcb/indexcnum;
    row_indices = row_indices(cbr>0.5);
    col_indices = col_indices(cbr>0.5);
    cbox = cbox(:,cbr>0.5);

    sridx = row_indices.' + datasy;
    scidx = col_indices.' + datasx;
    s_indices = sub2ind(size(a2), sridx, scidx);
    sbox=a2(s_indices);

    llscore=zeros(size(cbox));
    if max(cbox<=cnum1,[],'all')
        llscore(cbox<=cnum1) = gaussmf(cbox(cbox<=cnum1),[csig1 cnum1])...
        *cfactor1+cintersec1;
    end
    llscore(cbox>cnum1 & cbox<=cnum2) = cyfill;
    if max(cbox>cnum2,[],'all')
        llscore(cbox>cnum2) = gaussmf(cbox(cbox>cnum2),[csig2 cnum2])...
         *cfactor2+cintersec2;
    end
    clscore=sum(llscore,1,"omitmissing");

    snum1 = s_para(1);
    snum2 = s_para(2);
    ssig1 = s_para(3);
    sfactor1 = s_para(4);
    sintersec1 = s_para(5);
    ssig2 = s_para(6);
    sfactor2 = s_para(7);
    sintersec2 = s_para(8);
    syfill = s_para(9);

    if step==1
        con1 = sbox>=snum1 & sbox<=snum2;
        con2 = sbox>snum2;
    elseif step==2
        con1 = sbox>=snum1 & sbox<snum2;
        con2 = sbox>=snum2;
    end

    ssscore=zeros(size(sbox));
    ssscore(sbox<snum1) = syfill;
    if max(con1,[],'all')
            ssscore(con1) = ...
                gaussmf(sbox(con1),[ssig1 snum1])...
                 *sfactor1 + sintersec1;
    end
    if max(con2,[],'all')
        ssscore(con2)=gaussmf(sbox(con2),[ssig2 snum2])...
             *sfactor2 + sintersec2;
    end

    sdscore=sum(ssscore,1,"omitmissing");

    pretotscore=(sdscore+clscore);
    pretotscore=pretotscore./scorediv;
    pretotscore(pretotscore<0)=0;
    scoremt = zeros(401,401);
    lin_indices = sub2ind(size(a2), row_indices, col_indices);
    scoremt(lin_indices)=pretotscore;
end