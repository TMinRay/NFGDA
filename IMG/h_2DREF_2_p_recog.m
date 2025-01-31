run ../main/NF00_header

load('tmpdREFcdatax.mat','datacx');
load('tmpdREFcdatay.mat','datacy');
load('tmpdREFsdatax.mat','datasx');
load('tmpdREFsdatay.mat','datasy');
indexcnum=numel(datacx);
indexsnum=numel(datasx);
refcentersize=8;
refapart=5;
numINT =refcentersize;
nunapart=refapart;
    

for cindex=1:numel(ttable(:,1));
    
     PUTDAT=[ ttable(cindex,:) ];
     startm=startt(cindex)+1;
     endm=endt(cindex);
  
for m=startm:endm
      
    mDROTout=[matPATH '/IMG/DELZ/rotdelz' PUTDAT num2str(m,'%02i') '.mat'];
    load(mDROTout, 'rotitp');
    mSCRout=[matPATH '/IMG/DELZ/scrdelz' PUTDAT num2str(m,'%02i') '.mat'];
    llscore=[];
    ssscore=[];
    totscore=zeros(401,401,rotnum);
    totscore2=zeros(401,401);
    for i=1:rotnum 
    indi=rotdegree*(i-1);
    a2=rotitp(:,:,i);
    clscore=zeros(401,401);
    sdscore=zeros(401,401);
    inda2=(rotitp(:,:,i)>thrdREF);
    totscore(:,1:numINT,i)=nan;
    totscore(:,401-numINT+1:401,i)=nan;
    totscore(1:numINT,:,i)=nan;
    totscore(401-numINT+1:401,:,i)=nan;
    for ii=1+numINT:401-numINT   
        for jj=1+numINT:401-numINT
            cbox=[];
            sbox=[];
            if (inda2(ii,jj)==1)
            cbox=a2(ii-numINT:ii+numINT,jj);
            for kind=1:indexsnum
                sadd=[a2(ii+datasy(kind),jj+datasx(kind))];
                sbox=[sbox; sadd;];
                clear sadd;
            end
            indcb=((cbox)>=thrdREF);
            numtcb=sum(indcb(:));
            cbr=numtcb/indexcnum;
            if (cbr>=0.5)
            kkk=0;      
            h_2DREF_2_SCRfunc;
            else
            sdscore(ii,jj)=0; clscore(ii,jj)=0;  
            end
            else
            sdscore(ii,jj)=0; clscore(ii,jj)=0;  
            end
            clear cbox sbox;
        end %end Azimuth
    end
    pretotscore=(sdscore+clscore);
    totscore2=pretotscore./(2*17+1*18);
    for ki=1:401
        for kj=1:401
            if totscore2(ki,kj)<0
            totscore2(ki,kj)=0;
            end
        end
    end
    totscore(:,:,i)=totscore2;
    end
    clear totscore2;
    save(mSCRout, 'totscore');
    clear totscore;
    end
end

