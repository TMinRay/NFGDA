run ../main/NF00_header

load('tmpCELLdatax.mat','datacx');
load('tmpCELLdatay.mat','datacy');
load('tmpCELLdatax2.mat','datacx2');
load('tmpCELLdatay2.mat','datacy2');
indexnum=numel(datacx);
indexnum2=numel(datacx2);
crsize=5;
crsize2=crsize+2;
rsize=(crsize2);

for cindex=1:numel(ttable(:,1));
    
    PUTDAT=[ ttable(cindex,:) ];
    startm=startt(cindex)+1;
    endm=endt(cindex);
    for m=startm:endm
        mLINEout=[matPATH '/LINE/CELLZ/cellz' PUTDAT num2str(m,'%02i') '.mat'];
        load(mLINEout, 'CELLline');
        cellscr2=CELLline;
        tnum=1;
        indcell=(cellscr2>cellcsrthresh).*tnum;
        mLINEout2=[matPATH '/LINE/CELLZ/widecellz' PUTDAT num2str(m,'%02i') '.mat'];
        numINT =rsize+2;
        totscore(:,1:numINT)=nan;
        totscore(:,401-numINT+1:401)=nan;
        totscore(1:numINT,:)=nan;
        totscore(401-numINT+1:401,:)=nan;
        for ii=1+numINT:401-numINT
            for jj=1+numINT:401-numINT
                if (indcell(ii,jj)==tnum)
                cbox=[];
                for kind=1:indexnum
                    cadd=[indcell(ii+datacy(kind),jj+datacx(kind))];
                    cbox=[cbox; cadd;];
                    clear cadd;
                end
                indcb=((cbox)==tnum);
                numtcb=sum(indcb(:));
                cbr=numtcb/indexnum;
                if (cbr<1)
                    for kind=1:indexnum2
                        [cellscr2(ii+datacy2(kind),jj+datacx2(kind))]=1;
                    end

                end
                end
            end
        end
        widecellz=(cellscr2>0.5).*1;
        save(mLINEout2, 'widecellz');
    end
end
