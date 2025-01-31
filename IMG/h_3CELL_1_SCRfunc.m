s2xnum=[10 15];
s2ynum=[-3 1];
[n2 size2xnum]=size(s2xnum);
s2xdel=zeros(size2xnum,1);
s2ydel=zeros(size2xnum,1);
s2g=zeros(size2xnum,1);
s2gc=zeros(size2xnum,1);

for r=2:size2xnum
    s2xdel(r)=s2xnum(r)-s2xnum(r-1);
    s2ydel(r)=s2ynum(r)-s2ynum(r-1);
end    

for r=2:size2xnum
    s2g(r)=s2ydel(r)/s2xdel(r);
    s2gc(r)=s2ynum(r)-s2g(r)*s2xnum(r);    
end    


%  for pp=1:indexcnum
%      if cbox(pp)<=cnum1
%      llscore(pp)=ynum1;
%      elseif cbox(pp)>=cnum1 
%      llscore(pp)=-1*gaussmf(cbox(pp),[3 cnum1])*4+1;
%         
%      end
%          
%  end
% %  
%   if kkk==0   
%      clscore(ii,jj)=(nansum(llscore(:)));
%   end
  
  for pp=1:indexcnum
     if cbox(pp)<=s2xnum(1)
     llscore(pp)=s2ynum(1);
     elseif cbox(pp)>=s2xnum(1) & cbox(pp)<s2xnum(2)
     llscore(pp)=s2g(2)*cbox(pp)+s2gc(2);
     elseif cbox(pp)>=s2xnum(2)
     llscore(pp)=s2ynum(2);
     end
 end
    
  if kkk==0   
     clscore(ii,jj)=(sum(llscore(:),"omitmissing"));     
  end