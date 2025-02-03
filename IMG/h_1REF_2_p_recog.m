% run ../main/NF00_header
% 
% 
% load('tmpREFcdatax.mat','datacx');
% load('tmpREFcdatay.mat','datacy');
% load('tmpREFsdatax.mat','datasx');
% load('tmpREFsdatay.mat','datasy');
% indexcnum=numel(datacx);
% indexsnum=numel(datasx);
% refcentersize=8;
% refapart=5;
% 
% for cindex=1:numel(ttable(:,1));
% 
%     PUTDAT=[ ttable(cindex,:) ];
%     startm=startt(cindex)+1;
%     endm=endt(cindex);
%     for m=startm:endm
%         mROTout=[matPATH '/IMG/Z/rotz' PUTDAT num2str(m,'%02i') '.mat'];
%         load(mROTout, 'rotitp');
%         mSCRout=[matPATH '/IMG/Z/scrz' PUTDAT num2str(m,'%02i') '.mat'];
%         llscore=[];
%         ssscore=[];
%         totscore=zeros(401,401,rotnum); 
%         for i=1:rotnum
%             a2=rotitp(:,:,i);
%             numINT = refcentersize+2;
%             [row_indices, col_indices] = find(a2>thrREF);
%             inINT = (row_indices>numINT) & (row_indices<=401-numINT)...
%                 &(col_indices>numINT) & (col_indices<=401-numINT);
%             row_indices = row_indices(inINT);
%             col_indices = col_indices(inINT);
%             cridx = row_indices.' + datacy;
%             ccidx = col_indices.' + datacx;
%             c_indices = sub2ind(size(a2), cridx, ccidx);
%             cbox=a2(c_indices);
% 
%             indcb = cbox>thrREF;   
%             numtcb=sum(indcb,1);
%             cbr = numtcb/indexcnum;
%             row_indices = row_indices(cbr>0.5);
%             col_indices = col_indices(cbr>0.5);
%             cbox = cbox(:,cbr>0.5);
% 
%             sridx = row_indices.' + datasy;
%             scidx = col_indices.' + datasx;
%             s_indices = sub2ind(size(a2), sridx, scidx);
%             sbox=a2(s_indices);
%             kkk=0;   
%             h_1REF_2_SCRfunc;
% 
%             lin_indices = sub2ind(size(a2), row_indices, col_indices);
%             pretotscore=(sdscore+clscore);
%             pretotscore=pretotscore./(3*17+1*18);
%             pretotscore(pretotscore<0)=0;
%             totscore2 = zeros(401,401);
%             totscore2(lin_indices)=pretotscore;
%             totscore(:,:,i)=totscore2;
%         end
%         save(mSCRout, 'totscore');
%     end
% end


run ../main/NF00_header

refcentersize=8;
refapart=5;
numINT = refcentersize+2;


for cindex=1:numel(ttable(:,1));
    PUTDAT=[ ttable(cindex,:) ];
    startm=startt(cindex)+1;
    endm=endt(cindex);
    for m=startm:endm
        mROTout=[matPATH '/IMG/Z/rotz' PUTDAT num2str(m,'%02i') '.mat'];
        load(mROTout, 'rotitp');
        mSCRout=[matPATH '/IMG/Z/scrz' PUTDAT num2str(m,'%02i') '.mat'];
        totscore=zeros(401,401,rotnum); 
        for i=1:rotnum
            a2=rotitp(:,:,i);
            totscore(:,:,i)=gen_tot_score(a2, ...
                [15, 20, 3, 3, -1, 12, 4, -2, 3], ...
                [0, 5, 5, 2,-1, 5, 3,-3, 1], ...
                thrREF,10,(3*17+1*18),1);
        end
        save(mSCRout, 'totscore');
    end
end


