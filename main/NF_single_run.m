tic
NF00_header;

debugmat = true;

export_preds_dir = ['./tracking_points/nf_preds'];
xg = x(2:end,:);
yg = y(2:end,:);
% stdCELLcxout=['./NF03_STDdatax.mat'];
% stdCELLcyout=['./NF03_STDdatay.mat'];
% load(stdCELLcxout,'SD_displace_x');
% load(stdCELLcyout,'SD_displace_y');
SD_displace_x = [-1,-1,-1,0,0,0,1,1,1].';
SD_displace_y = [-1,0,1,-1,0,1,-1,0,1].';
indexnum=numel(SD_displace_x);
stdINT=2;

[gcol_indices,grow_indices]=meshgrid(1:401);
grow_indices = grow_indices(:);
gcol_indices = gcol_indices(:);

load('../IMG/tmpCELLdatax.mat','datacx');
load('../IMG/tmpCELLdatay.mat','datacy');
load('../IMG/tmpCELLdatax2.mat','datacx2');
load('../IMG/tmpCELLdatay2.mat','datacy2');
ndatacx = numel(datacx);
crsize = 5;
cellINT = crsize + 2;
widecellINT =crsize+4;

s2xnum = [10 15];
s2ynum = [-3 1];

s2xdel = s2xnum(2)-s2xnum(1);
s2ydel = s2ynum(2)-s2ynum(1);

s2g = s2ydel/s2xdel;
s2gc = s2ynum(2)-s2g*s2xnum(2);


h = fspecial('gaussian', 11,1);
debug_plot = true;
warning('off', 'all');
fuzzGST=readfis('./NF00ref_YHWANG.fis');
se = strel('ball',5,1);

avgINT =8;
disccx = zeros(17,17);
disccy = zeros(17,17);
for i=1:17
    cxout=['./movingavg/ccx' num2str(i,'%02i') '.mat'];
    cyout=['./movingavg/ccy' num2str(i,'%02i') '.mat'];
    load(cxout,'ccx');
    load(cyout,'ccy');
    disccx(i,:)=ccx;
    disccy(i,:)=ccy;
end
nccx=numel(ccx);
disccx = reshape(disccx,1,17,17);
disccy = reshape(disccy,1,17,17);

disccx = zeros(1,17,17);
disccy = zeros(1,17,17);
for i =1:17
disccx(:,i,:)=ceil([-8:8]*cosd(90/8*(i-1)))
disccy(:,i,:)=ceil([-8:8]*sind(90/8*(i-1)))
end
v6m_path = fullfile(matPATH,'POLAR',case_name);
v6m_list = {dir(fullfile(v6m_path,'polar*mat')).name};
% v6m_list([1,2]) = [];
PUTDAT = case_name;
% for cindex=1:numel(ttable(:,1));
%     PUTDAT=ttable(cindex,:);
cindex = 1;
startm=startt(cindex)+1;
endm=endt(cindex);

exp_preds_event = fullfile(export_preds_dir,PUTDAT);
if not(isfolder(exp_preds_event))
    mkdir(exp_preds_event);
end
% mrout1=[matPATH '/POLAR/polar' PUTDAT num2str(startm-1,'%02i') '.mat'];
% t0=load(mrout1, 'PARROT');

t0 = load(fullfile(v6m_path,v6m_list{startm-1}), 'PARROT');
for m=startm:endm
%%%%%%%%%%%%%%      NF01_convert_to_cartesian
    % mpolarout=[matPATH '/POLAR/polar' PUTDAT num2str(m,'%02i') '.mat'];
    % load(mpolarout, 'PARROT');
    load(fullfile(v6m_path,v6m_list{m}), 'PARROT');
    PARITP=zeros(401,401,6);
    varnum=[1 2 3 4 5 6];
    sdphi=zeros(400,720);
    phi=double(PARROT(:,:,4));
    phi(phi<0)=nan;
    phi(phi>360)=nan;
    buf =zeros(400,720,5);
    % for displaceR=-2:2
    %     buf(3:end-2,:,displaceR+3)=phi(3+displaceR:end-2+displaceR,:);
    % end
    % sdphi(3:end-2,:) = std(buf(3:end-2,:,:),0,3,"omitmissing");
    for displaceR=-2:2
        buf(5:Gate2-2,:,displaceR+3)=phi(5+displaceR:Gate2-2+displaceR,:);
    end
    sdphi(5:Gate2-2,:) = std(buf(5:Gate2-2,:,:),0,3,"omitmissing");
    PARROT(:,:,4)=sdphi;

    val = PARROT(2:end,:,1);
    F = scatteredInterpolant(xg(:),yg(:),val(:));
    for i=1:6
        val = PARROT(2:end,:,i);
        F.Values = double(val(:));
        PARITP(:,:,i) = F(xi2,yi2);
    end
    if debugmat
        mcartout=[matPATH '/DEBUG/cart' PUTDAT num2str(m,'%02i') '.mat'];
        save(mcartout, 'PARITP');
    end
%%%%%%%%%%%%%%      NF01_convert_to_cartesian
%%%%%%%%%%%%%%      NF02_calc_5by5_SD
    % mcartout=[matPATH '/CART/cart' PUTDAT num2str(m,'%02i') '.mat'];
    % load(mcartout, 'PARITP');
    stda = zeros(401,401,3);
    a2 = double(PARITP(:,:,2));
    [col_indices,row_indices]=meshgrid(stdINT+1:401-stdINT);
    cridx = row_indices(:).' + SD_displace_y;
    ccidx = col_indices(:).' + SD_displace_x;
    c_indices = sub2ind(size(a2), cridx, ccidx);
    cbox=a2(c_indices);
    indcb = ~isnan(cbox);
    numtcb=sum(indcb,1);
    cbr = numtcb/indexnum;
    row_indices = row_indices(cbr>=0.3);
    col_indices = col_indices(cbr>=0.3);
    cbox = cbox(:,cbr>=0.3);
    lin_indices = sub2ind(size(a2), row_indices, col_indices);
    cboxstd = std(cbox,0,1,"omitmissing");
    buf = zeros(401,401);
    buf(lin_indices)=cboxstd;
    stda(:,:,2) = buf;
    if debugmat
        mstdoutv=[matPATH '/DEBUG/std' PUTDAT num2str(m,'%02i') '.mat'];
        save(mstdoutv,'stda');
    end
%%%%%%%%%%%%%%      NF02_calc_5by5_SD

%%%%%%%%%%%%%%      NF04_calc_DeltaZ
    % mrout1=[matPATH '/POLAR/polar' PUTDAT num2str(m-1,'%02i') '.mat'];
    % load(mrout1, 'PARROT');
    zt0 = t0.PARROT(:,:,1);
    zt0(isnan(zt0)) = 0;

    % mrout2=[matPATH '/POLAR/polar' PUTDAT num2str(m,'%02i') '.mat']; 
    % load(mrout2, 'PARROT');
    zt1(:,:) = PARROT(:,:,1);
    zt1(isnan(zt1)) = 0;
    dif2 = zt1 - zt0;
    if debugmat
        dout=[matPATH '/DEBUG/delz' PUTDAT num2str(m,'%02i') '.mat'];
        save(dout,'dif2');
    end
%%%%%%%%%%%%%%      NF04_calc_DeltaZ

%%%%%%%%%%%%%%      NF05_calc_BETA_LINEFEATURE
%%%%%%%%%%%%%%      ../IMG/exe_3_img_exe.m
    oriz=double(PARROT(:,:,1));
    orirot=double(dif2);
    % degshit=rotdegree/angint;
    rotgz = [];
    rotz = zeros(size(oriz));
    ztotscore=zeros(401,401,rotnum);
    zoriginscore=zeros(401,401,rotnum);
    delztotscore=zeros(401,401,rotnum);
    originscore=zeros(401,401,rotnum);
    for i=1:rotnum
        indi=rotdegree*(i-1)/angint;
        rotz(:,1:720-indi)=oriz(:,indi+1:720);
        rotz(:,720-indi+1:720)=oriz(:,1:indi);
        rotgz(:,:,i) = griddata(x,y,rotz,xi2,yi2);

        roted(:,1:720-indi)=orirot(:,indi+1:720);
        roted(:,720-indi+1:720)=orirot(:,1:indi);
        rotitp(:,:,i) = griddata(x,y,roted,xi2,yi2);

        ztotscore(:,:,i)=gen_tot_score(rotgz(:,:,i), ...
            [15, 20, 3, 3, -1, 12, 4, -2, 3], ...
            [0, 5, 5, 2,-1, 5, 3,-3, 1], ...
            thrREF,10,(3*17+1*18),1);

        delztotscore(:,:,i)=gen_tot_score(rotitp(:,:,i), ...
            [5,10,4,3,-2,9,4,-3,2], ...
            [-10,5,5,2,-1,8,2,-3,1], ...
            thrdREF, 8, (2*17+1*18),2);

        origindeg = rotbackrad*(i-1);
        zoriginscore(:,:,i) = rot_score_back(ztotscore(:,:,i),origindeg,grow_indices,gcol_indices);
        originscore(:,:,i) = rot_score_back(delztotscore(:,:,i),origindeg,grow_indices,gcol_indices);

    end
    linez=max(zoriginscore,[],3);
    linedelz=max(originscore,[],3);

    a2=rotgz(:,:,1);
    clscore=zeros(401,401);
    [row_indices, col_indices] = find(a2>cellthresh);
    inINT = (row_indices>cellINT) & (row_indices<=401-cellINT)...
        &(col_indices>cellINT) & (col_indices<=401-cellINT);
    row_indices = row_indices(inINT);
    col_indices = col_indices(inINT);
    cridx = row_indices.' + datacy;
    ccidx = col_indices.' + datacx;
    c_indices = sub2ind(size(a2), cridx, ccidx);
    cbox=a2(c_indices);
    indcb = cbox>cellthresh;
    numtcb=sum(indcb,1);
    cbr = numtcb/ndatacx;
    row_indices = row_indices(cbr>cbcellthrsh);
    col_indices = col_indices(cbr>cbcellthrsh);
    cbox = cbox(:,cbr>cbcellthrsh);

    llscore=zeros(size(cbox));
    llscore(cbox<=s2xnum(1))=s2ynum(1);
    pp = cbox>=s2xnum(1) & cbox<s2xnum(2);
    if max(pp,[],'all')
        llscore(pp)=s2g*cbox(pp)+s2gc;
    end
    llscore(cbox>=s2xnum(2))=s2ynum(2);
    clscore = sum(llscore,1,"omitmissing");

    clscore = clscore/113;
    scoremt = zeros(401,401);
    lin_indices = sub2ind(size(a2), row_indices, col_indices);
    scoremt(lin_indices) = clscore;
    totscore = zeros(401,401);
    totscore(9:392,9:392) = scoremt(9:392,9:392);
    CELLline=medfilt2(totscore,[11 11]);

    a2 = CELLline;
    [row_indices, col_indices] = find(a2>cellcsrthresh);
    inINT = (row_indices>widecellINT) & (row_indices<=401-widecellINT)...
        &(col_indices>widecellINT) & (col_indices<=401-widecellINT);
    row_indices = row_indices(inINT);
    col_indices = col_indices(inINT);
    cridx = row_indices.' + datacy;
    ccidx = col_indices.' + datacx;
    c_indices = sub2ind(size(a2), cridx, ccidx);
    cbox=double(a2(c_indices)>cellcsrthresh);
    numtcb=sum(cbox,1);
    cbr = numtcb/ndatacx;
    row_indices = row_indices(cbr<1);
    col_indices = col_indices(cbr<1);
    cridx = row_indices.' + datacy2;
    ccidx = col_indices.' + datacx2;
    c_indices = sub2ind(size(a2), cridx, ccidx);
    a2(c_indices)=1;
    widecellz=double(a2>0.5);
%%%%%%%%%%%%%%      ../IMG/exe_3_img_exe.m
    pbeta = (linez+linedelz)./2;
    pbeta(isnan(PARITP(:,:,1)))=nan;
    beta = pbeta-widecellz;
    beta(beta<0)=0;
%%%%%%%%%%%%%%      NF05_calc_BETA_LINEFEATURE

%%%%%%%%%%%%%%      NF06_calc_6variables_preprocessing
    inputNF=zeros(401,401,12);
    % inputNF(:,:,1)=PARITP(:,:,1);
    % inputNF(:,:,2)=beta;
    % inputNF(:,:,3)=PARITP(:,:,6);
    % inputNF(:,:,4)=PARITP(:,:,5);
    % inputNF(:,:,5)=stda(:,:,2);
    % inputNF(:,:,6)=PARITP(:,:,4);

    %%%%%%%%%%%%%%      direct layout for fis input
    inputNF(:,:,1)=beta;
    inputNF(:,:,2)=PARITP(:,:,1);
    inputNF(:,:,3)=PARITP(:,:,5);
    inputNF(:,:,4)=PARITP(:,:,6);
    inputNF(:,:,5)=stda(:,:,2);
    inputNF(:,:,6)=PARITP(:,:,4);
    %%%%%%%%%%%%%%      direct layout for fis input
    pnan=(isnan(inputNF));
    pnansum=max(pnan,[],3);
    inputNF(repmat(pnansum, 1, 1, size(inputNF, 3))) = nan;
%%%%%%%%%%%%%%      NF06_calc_6variables_preprocessing

%%%%%%%%%%%%%%      NF10_evalfuzzy_and_skel
    dINPUT=reshape(inputNF(:,:,1:6),401*401,6);
    % [incoef,outcoef,rulelogic] = fis2python(fuzzGST,'NF00ref_YHWANG_fis4python.mat');
    PREGST=evalfis(dINPUT,fuzzGST);
    hGST=reshape(PREGST,401,401);
    numthr=0.24;
    hGST = double(hGST>=numthr);
    hh=hGST;
    hGST=medfilt2(hGST,[3 3]);
    smoothedhGST = imfilter(hGST,h,'replicate');
    smoothedhGST(smoothedhGST<0.1)=nan;
    skel_nfout = bwmorph(double(smoothedhGST>0.3), 'skel', inf);

    if debugmat
        mGSTout=[matPATH '/DEBUG/nf' PUTDAT num2str(m,'%02i') '.mat'];
        save(mGSTout,'hh','hGST','smoothedhGST','skel_nfout');
    end
%%%%%%%%%%%%%%      NF10_evalfuzzy_and_skel

%%%%%%%%%%%%%%      NF07_obtaining_evaluation_box

    if no_eval
        evalbox = zeros(401,401);
    else
        mhandpick=[ matPATH '/HANDPICK/handpick' PUTDAT num2str(m,'%02i') '.mat'];
        load(mhandpick,'evalbox');
        % load(mhandpick,'handpick');
        % pevalbox = bwmorph(double(handpick), 'skel', inf);
        % evalbox = double(imdilate(double(pevalbox),se)>1);
    end

    if debugmat
        mevalbox=[ matPATH '/DEBUG/newevalbox' PUTDAT num2str(m,'%02i') '.mat'];
        save(mevalbox,'evalbox');
    end
%%%%%%%%%%%%%%      NF07_obtaining_evaluation_box

%%%%%%%%%%%%%%      NF11_postprocessing_movingavg
    ppresult=zeros(401,401);
    a2=hGST;

    [row_indices, col_indices] = find(a2>0);
    inINT = (row_indices>avgINT) & (row_indices<=401-avgINT)...
        &(col_indices>avgINT) & (col_indices<=401-avgINT);
    row_indices = row_indices(inINT);
    col_indices = col_indices(inINT);

    cridx = reshape(row_indices,numel(row_indices),1,1) + disccy;
    ccidx = reshape(col_indices,numel(col_indices),1,1) + disccx;
    c_indices = sub2ind(size(a2), cridx, ccidx);
    cbox=a2(c_indices);
    indcb = cbox>0;
    numtcb=sum(indcb,3);
    cbr = numtcb/nccx;
    validcenter = max(cbr>0.1,[],2);
    row_indices = row_indices(validcenter);
    col_indices = col_indices(validcenter);
    cbox = cbox(validcenter,:,:);
    mc = mean(cbox,3,"omitmissing");
    mc(~(cbr>0.1))=0;
    buf = zeros(401,401);
    lin_indices = sub2ind(size(a2), row_indices, col_indices);
    buf(lin_indices) = max(mc,[],2);
    ppresult(avgINT+1:401-avgINT,avgINT+1:401-avgINT) = buf(avgINT+1:401-avgINT,avgINT+1:401-avgINT);
%%%%%%%%%%%%%%      NF11_postprocessing_movingavg

%%%%%%%%%%%%%%      NF12_final_output_skel
    pskel_nfout = double(imdilate(double(ppresult>=0.6),se)>1);
    pskel = pskel_nfout.*hh;
    skel_nfout = bwmorph(double(pskel), 'skel', inf);
    skel_nfout2 = bwareaopen(skel_nfout,10);
%%%%%%%%%%%%%%      NF12_final_output_skel

%%%%%%%%%%%%%%      NFFIG_03FINAL_output
    % REF=inputNF(:,:,1);
    REF=inputNF(:,:,2);
    nfout=skel_nfout2;
    matout = fullfile(exp_preds_event,['nf_pred' v6m_list{m}(6:end)]);
    % matout = [exp_preds_event '\' num2str(m,'%02i') '.mat'];
    save(matout,"xi2","yi2","REF","nfout","evalbox");
%%%%%%%%%%%%%%      NFFIG_03FINAL_output

    if debugmat
        mROTout=[matPATH '/DEBUG/' PUTDAT num2str(m,'%02i') '.mat'];
        save(mROTout, 'rotgz', 'rotitp', 'ztotscore', 'delztotscore', ...
            'zoriginscore', 'linez', 'linedelz', 'CELLline', 'widecellz', ...
            'beta', 'inputNF', 'ppresult', 'skel_nfout', 'skel_nfout2');
    end
    t0.PARROT = PARROT;
end
% end
toc