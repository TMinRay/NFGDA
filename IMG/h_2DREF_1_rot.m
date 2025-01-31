run ../main/NF00_header

for cindex=1:numel(ttable(:,1));
    
     PUTDAT=[ ttable(cindex,:) ];
     startm=startt(cindex)+1;
     endm=endt(cindex);
  
for m=startm:endm
    t=m;
    
    dout=[matPATH '/DELZ/delz' PUTDAT num2str(m,'%02i') '.mat'];
    load(dout,'dif2');
    orirot=double(dif2);    
    degshit=rotdegree/angint;
    mDROTout=[matPATH '/IMG/DELZ/rotdelZ' PUTDAT num2str(m,'%02i') '.mat'];
    roted=[]; 
    rotitp=[]; 
    
    for i=1:rotnum
        indi=rotdegree*(i-1);
        if i==1
        roted(:,1:720,i)=orirot(:,1:720);
        else
        roted(:,1:720-degshit*(i-1),i)=orirot(:,degshit*(i-1)+1:720);
        roted(:,720-degshit*(i-1)+1:720,i)=orirot(:,1:degshit*(i-1));
        end
        rotitp(:,:,i) = griddata(x,y,roted(:,:,i),xi2,yi2);
        clear roted;  
     end
     save(mDROTout, 'rotitp');
     clear rotitp;
end
end
     