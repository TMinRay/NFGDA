


cnum1=5;
cnum2=10;

ynum1=2;
snum1=-10;
snum2=5;
synum1=2;

% kkk=0;


 for pp=1:indexcnum
     if cbox(pp)<=cnum1
     llscore(pp)=gaussmf(cbox(pp),[4 cnum1])...
         *(2*ynum1-1)-ynum1;
     elseif cbox(pp)>cnum1 & cbox(pp)<=cnum2
     llscore(pp)=ynum1;
     elseif cbox(pp)>cnum2
     llscore(pp)=gaussmf(cbox(pp),[9 cnum2])...
         *(2*ynum1+1-1)-(ynum1)-1;
     end
         
 end
  if kkk==0   
     clscore(ii,jj)=(sum(llscore(:),"omitmissing"));
  end
 for ppp=1:indexsnum
     if sbox(ppp)<snum1
         ssscore(ppp)=synum1-1;
     elseif sbox(ppp)>=snum1 & sbox(ppp)<snum2
     ssscore(ppp)=gaussmf(sbox(ppp),[5 snum1])...
         *(synum1)-(synum1)+1;
         
     elseif sbox(ppp)>=snum2
     ssscore(ppp)=gaussmf(sbox(ppp),[8 snum2])...
         *(synum1)-(synum1)-1;
     end
 end
    if kkk==0
     sdscore(ii,jj)=(sum(ssscore(:),"omitmissing"));     
    end
   