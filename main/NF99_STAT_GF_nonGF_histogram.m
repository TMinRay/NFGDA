clear all;
close all;
NF00_header;
totnumag=0;
totnumao=0;
totg=cell(1,6);
toto=totg;

for cindex=1:numel(ttable(:,1));

    PUTDAT=ttable(cindex,:);
    startm=startt(cindex)+1;
    endm=endt(cindex);
  
    for m=startm:endm
        tot=zeros(401,401,6);
        minputNF=[matPATH '/CART/inputNF' PUTDAT num2str(m,'%02i') '.mat']; 
%         exist(minputNF)
        load(minputNF,'inputNF')
        nokink=inputNF(:,:,3);
        nokink(nokink>=7.98)=nan;
        inputNF(:,:,3)=nokink;
        clear nokink;
        nokink=inputNF(:,:,5);
        nokink(nokink==0)=nan;
        inputNF(:,:,5)=nokink;
        clear nokink;
        
        tot=inputNF;
        clear inputNF;
        mevalbox=[ matPATH '/EVAL/newevalbox' PUTDAT num2str(m,'%02i') '.mat']; 
        load(mevalbox,'evalbox');  
%         figure(m)
%         pcolor(evalbox)
%         shading flat
    
        gustin=evalbox;
        oth=logical(gustin.*(-1)+1);
        gst=logical(gustin);
        clear gustin;
        for k=1:6
            aa=tot(:,:,k);
            agst=aa(gst);
            aoth=aa(oth);
            agnidx=~isnan(agst);
            if k==2
            agnidx2=(agst~=0);  
            for kkk=1:numel(agnidx2)
                agnidx(kkk)=agnidx(kkk)*agnidx2(kkk);          
            end     
            end
            aonidx=~isnan(aoth);
            totg{k}=[totg{k}; agst(agnidx);];
            toto{k}=[toto{k}; aoth(aonidx);];
            clear aa agst aoth;
        end
        clear tot;
    end
end
mANALGSTout=[matPATH '/HIST/stat2regions.mat'];
save(mANALGSTout,'totg','toto');
