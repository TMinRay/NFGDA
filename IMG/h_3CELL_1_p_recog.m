run ../main/NF00_header

load('tmpCELLdatax.mat','datacx');
load('tmpCELLdatay.mat','datacy');
indexcnum=numel(datacx);
crsize=5;
numINT =crsize+2;

   
for cindex=1:numel(ttable(:,1));
    
    PUTDAT=[ ttable(cindex,:) ];
    startm=startt(cindex)+1;
    endm=endt(cindex);
    for m=startm:endm
   
     
    mROTout=[matPATH '/IMG/Z/rotz' PUTDAT num2str(m,'%02i') '.mat'];
    load(mROTout, 'rotitp');
    mLINEout=[matPATH '/LINE/CELLZ/cellz' PUTDAT num2str(m,'%02i') '.mat'];
    llscore=[];
    ssscore=[];
    totscore=[];
    for i=1:1
        indi=rotdegree*(i-1);
        a2=rotitp(:,:,i);
        clscore=zeros(401,401);
        sdscore=zeros(401,401);
        inda2=(a2>cellthresh);
        
        totscore(:,1:numINT)=nan;
        totscore(:,401-numINT+1:401)=nan;
        totscore(1:numINT,:)=nan;
        totscore(401-numINT+1:401,:)=nan;
        for ii=1+numINT:401-numINT
            for jj=1+numINT:401-numINT
                if (inda2(ii,jj)==1)
                cbox=[];
                    for kind=1:indexcnum
                        cadd=[a2(ii+datacy(kind),jj+datacx(kind))];
                        cbox=[cbox; cadd;];
                        clear cadd;
                    end
                    indcb=((cbox)>cellthresh);
                    numtcb=sum(indcb(:));
                    cbr=numtcb/indexcnum;
                    if (cbr>cbcellthrsh)
                    numtcb=indexcnum;
                    kkk=0;
                    h_3CELL_1_SCRfunc;
                    end
                else
                sdscore(ii,jj)=0; clscore(ii,jj)=0;
                end
            end %end Azimuth
        end
        pretotscore=clscore;
        totscore(:,:,i)=pretotscore./113;
        clear sdscore clscore;
        totscore(1:8,:,i)=0;
        totscore(393:401,:,i)=0;
        totscore(:,393:401,i)=0;
        totscore(:,1:8,i)=0;
    end
    CELLline=medfilt2(totscore,[11 11]);
    save(mLINEout, 'CELLline');
    end
end
