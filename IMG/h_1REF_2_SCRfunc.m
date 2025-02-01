cnum1=15;
cnum2=20;
ynum1=2;
snum1=0;
snum2=5;
synum1=1;
llscore=zeros(size(cbox));
if max(cbox<=cnum1)
    llscore(cbox<=cnum1) = gaussmf(cbox(cbox<=cnum1),[3 cnum1])*(2*ynum1-1)-1;
end
llscore(cbox>cnum1 & cbox<=cnum2) = ynum1+1;
if max(cbox>cnum2)
    llscore(cbox>cnum2) = gaussmf(cbox(cbox>cnum2),[12 cnum2])...
     *(2*ynum1)-(ynum1);
end
% llscore=zeros(size(cbox));
% for pp=1:indexcnum
% 
%      if cbox(pp)<=cnum1
%      llscore(pp)=gaussmf(cbox(pp),[3 cnum1])*(2*ynum1-1)-1;
% %      llscore(pp)=gaussmf(cbox(pp),[3 cnum1])...
% 
%      elseif cbox(pp)>cnum1 & cbox(pp)<=cnum2
%      llscore(pp)=ynum1+1;
%      elseif cbox(pp)>cnum2
%      llscore(pp)=gaussmf(cbox(pp),[12 cnum2])...
%          *(2*ynum1)-(ynum1);
%      end
% end

if kkk==0       
    clscore(ii,jj)=sum(llscore,"omitmissing");
end

ssscore=zeros(size(sbox));
ssscore(sbox<snum1) = synum1;
if max(sbox>=snum1 & sbox<=snum2)
    ssscore(sbox>=snum1 & sbox<=snum2) = ...
        gaussmf(sbox(sbox>=snum1 & sbox<=snum2),[5 snum1])...
         *(synum1+1)-1;
end
if max(sbox>snum2)
    ssscore(sbox>snum2)=gaussmf(sbox(sbox>snum2),[5 snum2])...
         *(synum1+2)-3;
end
% ssscore=zeros(size(sbox));
%   for ppp=1:indexsnum
%      if sbox(ppp)<snum1
%          ssscore(ppp)=synum1;
%      elseif sbox(ppp)>=snum1 & sbox(ppp)<=snum2
%      ssscore(ppp)=gaussmf(sbox(ppp),[5 snum1])...
%          *(synum1+1)-1;
% %            elseif sbox(ppp)>=snum1 & sbox(ppp)<=snum2
% %      ssscore(ppp)=gaussmf(sbox(ppp),[5 snum1])...
% %          *(synum1+2)-1;   
%      elseif sbox(ppp)>=snum2
%      ssscore(ppp)=gaussmf(sbox(ppp),[5 snum2])...
%          *(synum1+2)-3;
% %      elseif sbox(ppp)>snum2
% %      ssscore(ppp)=gaussmf(sbox(ppp),[10 snum2])...
% %          *(synum1+2)-3-1;
% 
% 
%      end
%  end
 
if kkk==0   
    sdscore(ii,jj)=sum(ssscore,"omitmissing");
end