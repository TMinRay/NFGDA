clear all;
close all;
run ../../hisnameys/header
tmpCELLcxout2=['./tmpmat/tmpCELLdatax2.mat'];
tmpCELLcyout2=['./tmpmat/tmpCELLdatay2.mat'];
tmppout=['./tmppng/tmpCELL2.png'];

   
rsize=crsize2+2;


gsize=rsize*2+1;
orip=rsize+1;
ori=zeros(gsize+1,gsize+1);
ori(orip,orip)=2;

 se = strel('disk',rsize);
 dori=imdilate(ori,se);
%  dori(orip,orip)=2;
dori(rsize+1,1)=2;
dori(rsize+1,rsize*2+1)=2;
dori(1,rsize+1)=2;
dori(rsize*2+1,rsize+1)=2;

 for j=1:gsize+1
     for k=1:gsize+1
 x3(j,k)=k-(rsize+1);
 y3(j,k)=j-(rsize+1);
     end
 end
 
  for j=1:gsize+1
     for k=1:gsize+1
         if dori(j,k)==0
             dori(j,k)=nan;
         end
     end
  end
  figure(1)
  pcolor(x3,y3,dori);
  

 caxis([0 2])
 set(figure(1),'Position', [ 100 100 400 400 ] );
 set(figure(1), 'PaperPositionMode','auto') 
 title('CELL tmp.')
 saveas(gcf,tmppout,'png');


 
 [indy,indx]=find(dori==2);
 
 datacx2=indx-(rsize+1);
 datacy2=[y3(indy)];
 
 save(tmpCELLcxout2,'datacx2');
 save(tmpCELLcyout2,'datacy2');
 
 
%  indexnum=numel(indy);

 
 
 
