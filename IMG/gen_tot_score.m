function [scoremt] = gen_tot_score(a2,cys_num,thrREF)
    % cys_num = [cnum1,cnum2,ynum1,snum1,snum2,synum1];
    cnum1 = cys_num(1);
    cnum2 = cys_num(2);
    ynum1 = cys_num(3);
    snum1 = cys_num(4);
    snum2 = cys_num(5);
    synum1 = cys_num(6);
    numINT = 10;
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
        llscore(cbox<=cnum1) = gaussmf(cbox(cbox<=cnum1),[3 cnum1])*(2*ynum1-1)-1;
    end
    llscore(cbox>cnum1 & cbox<=cnum2) = ynum1+1;
    if max(cbox>cnum2,[],'all')
        llscore(cbox>cnum2) = gaussmf(cbox(cbox>cnum2),[12 cnum2])...
         *(2*ynum1)-(ynum1);
    end
    clscore=sum(llscore,1,"omitmissing");


    ssscore=zeros(size(sbox));
    ssscore(sbox<snum1) = synum1;
    if max(sbox>=snum1 & sbox<=snum2,[],'all')
        ssscore(sbox>=snum1 & sbox<=snum2) = ...
            gaussmf(sbox(sbox>=snum1 & sbox<=snum2),[5 snum1])...
             *(synum1+1)-1;
    end
    if max(sbox>snum2,[],'all')
        ssscore(sbox>snum2)=gaussmf(sbox(sbox>snum2),[5 snum2])...
             *(synum1+2)-3;
    end

    sdscore=sum(ssscore,1,"omitmissing");

    pretotscore=(sdscore+clscore);
    pretotscore=pretotscore./(3*17+1*18);
    pretotscore(pretotscore<0)=0;
    scoremt = zeros(401,401);
    lin_indices = sub2ind(size(a2), row_indices, col_indices);
    scoremt(lin_indices)=pretotscore;
end