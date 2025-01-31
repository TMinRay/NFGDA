NF00_header;
[x2,y2]=meshgrid(-100:0.5:100,-100:0.5:100);
indcirc=zeros(401,401);
for i=1:401
    for j=1:401
        if sqrt(x2(i,j)^2+y2(i,j)^2)<=70
        indcirc(i,j)=1;
        else
        indcirc(i,j)=0;
        end
    end
end


    lnum=0;
    for pp=1:numel(ttable(:,1));
      lnum=lnum+endt(pp)-startt(pp);
    end
    lscr=zeros(lnum,numel(ttable(:,1)));
    oo=1;
    for cindex=1:numel(ttable(:,1));
    PUTDAT=ttable(cindex,:);
    startm=startt(cindex)+1;
    endm=endt(cindex);
    for m=startm:endm
        lscr(oo,1)=(cindex);
        lscr(oo,2)=m;
        oo=oo+1;
    end
    end
    oo=1;
    t=1;
    for cindex=1:numel(ttable(:,1));

    PUTDAT=ttable(cindex,:);
    startm=startt(cindex)+1;
    endm=endt(cindex);
    for m=startm:endm
        
         mevalbox=[ matPATH '/EVAL/newevalbox' PUTDAT num2str(m,'%02i') '.mat']; 
        if exist(mevalbox)>0
        load(mevalbox,'evalbox');   
        else
           evalbox=zeros(401,401); 
        end
%         lscr(oo,numel(ttable(:,1)))=sum(evalbox(:))/10;
        mindxy=['../mat/OUT/MIGFAfinal' PUTDAT num2str(m,'%02i') '.mat'];
        if exist(mindxy)>0
        load(mindxy,'skel_nfout');    
        else
            skel_nfout=zeros(401,401);
        end


        evalbox2=evalbox.*indcirc;            
        farray2=skel_nfout.*indcirc;
    
        truescr=evalbox2;
        falsescr=-1.*(truescr-1);
        tp=farray2.*truescr;  
        fp=farray2.*falsescr;  
    
        if sum(truescr(:))>0 & sum(tp(:))>0
            lscr(oo,3)=1;
            lscr(oo,7)=sum(tp(:));
            lscr(oo,8)=sum(fp(:));
            lscr(oo,9)=(sum(truescr(:)))/10;
            
            elseif sum(truescr(:))>0 & sum(tp(:))==0   
            lscr(oo,4)=1;
            elseif sum(truescr(:))==0 & sum(tp(:))==0 & sum(fp(:))==0   
            lscr(oo,6)=1;
            elseif sum(truescr(:))==0 & sum(fp(:))>0   
            lscr(oo,5)=1;
%             sum(fp(:));
%             if sum(fp(:))~=(401*401)
            lscr(oo,8)=sum(fp(:));
%             end
            
        end
        clear trueregion truescr falsescr fp tp ;
        oo=oo+1;   
    end
    end
    lscrfile2=[matPATH '/EVAL/neweval_exmigfa.mat'];
    save(lscrfile2,'lscr')   
% end
