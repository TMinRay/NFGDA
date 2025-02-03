s2xnum=[10 15];
s2ynum=[-3 1];
[n2 size2xnum]=size(s2xnum);

s2xdel=s2xnum(2)-s2xnum(1);
s2ydel=s2ynum(2)-s2ynum(1);

s2g=s2ydel/s2xdel;
s2gc=s2ynum(2)-s2g*s2xnum(2);


llscore=zeros(size(cbox));
llscore(cbox<=s2xnum(1))=s2ynum(1);
pp = cbox>=s2xnum(1) & cbox<s2xnum(2);
llscore(pp)=s2g(2)*cbox(pp)+s2gc(2);
llscore(cbox>=s2xnum(2))=s2ynum(2);

 %  for pp=1:indexcnum
 %     if cbox(pp)<=s2xnum(1)
 %     llscore(pp)=s2ynum(1);
 %     elseif cbox(pp)>=s2xnum(1) & cbox(pp)<s2xnum(2)
 %     llscore(pp)=s2g(2)*cbox(pp)+s2gc(2);
 %     elseif cbox(pp)>=s2xnum(2)
 %     llscore(pp)=s2ynum(2);
 %     end
 % end
    
if kkk==0
    clscore(ii,jj)=sum(llscore,"omitmissing");
end