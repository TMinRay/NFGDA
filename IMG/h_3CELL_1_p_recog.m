run ../main/NF00_header

load('tmpCELLdatax.mat','datacx');
load('tmpCELLdatay.mat','datacy');
indexcnum=numel(datacx);
crsize=5;
numINT =crsize+2;

s2xnum=[10 15];
s2ynum=[-3 1];

s2xdel=s2xnum(2)-s2xnum(1);
s2ydel=s2ynum(2)-s2ynum(1);

s2g=s2ydel/s2xdel;
s2gc=s2ynum(2)-s2g*s2xnum(2);

for cindex=1:numel(ttable(:,1));
    
    PUTDAT=[ ttable(cindex,:) ];
    startm=startt(cindex)+1;
    endm=endt(cindex);
    for m=startm:endm
    mROTout=[matPATH '/IMG/Z/rotz' PUTDAT num2str(m,'%02i') '.mat'];
    load(mROTout, 'rotitp');
    mLINEout=[matPATH '/LINE/CELLZ/cellz' PUTDAT num2str(m,'%02i') '.mat'];
    for i=1:1
        indi=rotdegree*(i-1);
        a2=rotitp(:,:,i);
        clscore=zeros(401,401);
        % sdscore=zeros(401,401);
        inda2=(a2>cellthresh);

        [row_indices, col_indices] = find(a2>cellthresh);
        inINT = (row_indices>numINT) & (row_indices<=401-numINT)...
            &(col_indices>numINT) & (col_indices<=401-numINT);
        row_indices = row_indices(inINT);
        col_indices = col_indices(inINT);
        cridx = row_indices.' + datacy;
        ccidx = col_indices.' + datacx;
        c_indices = sub2ind(size(a2), cridx, ccidx);
        cbox=a2(c_indices);
        indcb = cbox>cellthresh;
        numtcb=sum(indcb,1);
        cbr = numtcb/indexcnum;
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


        % totscore(:,1:numINT)=nan;
        % totscore(:,401-numINT+1:401)=nan;
        % totscore(1:numINT,:)=nan;
        % totscore(401-numINT+1:401,:)=nan;
        % for ii=1+numINT:401-numINT
        %     for jj=1+numINT:401-numINT
        %         if (inda2(ii,jj)==1)
        %         cbox=[];
        %             for kind=1:indexcnum
        %                 cadd=[a2(ii+datacy(kind),jj+datacx(kind))];
        %                 cbox=[cbox; cadd;];
        %                 clear cadd;
        %             end
        %             indcb=((cbox)>cellthresh);
        %             numtcb=sum(indcb(:));
        %             cbr=numtcb/indexcnum;
        %             if (cbr>cbcellthrsh)
        %             numtcb=indexcnum;
        %             % kkk=0;

        %             % h_3CELL_1_SCRfunc;


        %             end
        %         else
        %         clscore(ii,jj)=0;
        %         end
        %     end %end Azimuth
        % end
        % pretotscore=clscore;
        % totscore(:,:,i)=pretotscore./113;
        % clear clscore;
        

        clscore = clscore/113;
        scoremt = zeros(401,401);
        lin_indices = sub2ind(size(a2), row_indices, col_indices);
        scoremt(lin_indices) = clscore;
        totscore(:,:,i) = scoremt;
    end
    totscore(1:8,:,:)=0;
    totscore(393:401,:,:)=0;
    totscore(:,393:401,:)=0;
    totscore(:,1:8,:)=0;
    CELLline=medfilt2(totscore,[11 11]);
    save(mLINEout, 'CELLline');
    end
end
