run '../NF00_header';
[xi3, yi3]=meshgrid(-100:0.1:100,-100:0.1:100);

origindeg=deg2rad(90);    
mat11=cos(origindeg);
mat12=sin(origindeg);
mat21=-1*mat12;
mat22=mat11;
backprocess=[mat11 mat12; mat21 mat22;];
se=strel('ball',6,0);
    
for cindex=1:numel(ttable(:,1));

    PUTDAT=ttable(cindex,:);
    startm=startt(cindex)+1;
    endm=endt(cindex);
for m=startm:endm
% for m=startm:startm

    t=m;
    
Mresultf=[ '../../mat/MIGFAOUT/new' PUTDAT  ...
    num2str(t,'%02i')  '.mat'];

mfile=[ '../../MIGFAdat/' PUTDAT ...
                  num2str(m,'%02i') '.dat'];

 fid   = fopen(mfile,'r');
 M=textscan(fid,'%d%f%f','delimiter','\t');
 checkval=M{1};
 if sum(checkval(:))>1
     
 Mlat=(M{3}-Mlath(cindex)).*100;
 Mlon=(M{2}-Mlonh(cindex)).*100;
 fclose(fid);

      txnew=ceil(Mlat)-0.1;
      tynew=ceil(Mlon)-0.1;
     
      indx=round(1+((txnew+100)./0.1));
      indy=round(1+((tynew+100)./0.1));
      
      numindx=numel(indx);

       ptrueregion=zeros(2001,2001).*1;

       for ki=1:numindx
%           ptrueregion(indy(ki)-1:indy(ki)+1,indx(ki)-1:indx(ki)+1)=1;
       ptrueregion(indx(ki),indy(ki))=1;
       end
       trueregion=imdilate(ptrueregion,se);
       clear indx indy txnew tynew;
      
      Mresult=trueregion;
      clear trueregion;
      xx=xi3(logical(Mresult));
      yy=yi3(logical(Mresult));
      
%        figure(m)
%       for i=1:numel(M{1})
%       plot(xx,yy,'k.');
%       shading flat;
%        hold on;
%       end
%       close(m)
      
      xlim([-70 70]);
      ylim([-70 70]);
 else
     Mresult=zeros(2001,2001);
 end
%  clf
 save(Mresultf,'Mresult');
 
%  clear M Mlat Mlon Mresult;
     
end

end

  


