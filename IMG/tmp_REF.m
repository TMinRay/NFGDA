clear all;
close all;
run ../hisnameys/header


tmpREFcxout=['./tmpmat/tmpREFcdatax.mat'];
tmpREFcyout=['./tmpmat/tmpREFcdatay.mat'];
tmpREFsxout=['./tmpmat/tmpREFsdatax.mat'];
tmpREFsyout=['./tmpmat/tmpREFsdatay.mat'];
tmppout=['./tmppng/tmpREF.png'];


   
linesize=refcentersize;
apartsize=refapart-1;

gsize=linesize*2+1;
orip=apartsize*2+1;
ori=zeros(gsize+1,orip+1);

ori(1:gsize,apartsize+1)=2;

ori(2*(1:linesize),apartsize-(apartsize-1))=1;
ori(2*(1:linesize),apartsize*2+1)=1;
ori(linesize+1,apartsize-(apartsize-1))=1;
ori(linesize+1,apartsize*2+1)=1;

 for j=1:gsize+1
     for k=1:orip+1
 x3(j,k)=k-(apartsize+1);
 y3(j,k)=j-(linesize+1);
     end
 end
 
  for j=1:gsize+1
     for k=1:orip+1
         if ori(j,k)==0
             ori(j,k)=nan;
         end
     end
  end
 figure(1) 
 pcolor(x3,y3,(ori));
 set(figure(1),'Position', [ 100 100 400 800 ] );
 set(figure(1), 'PaperPositionMode','auto') 
 title('REF tmp.')
 saveas(gcf,tmppout,'png');
 
 
  [indcy,indcx]=find(ori==2);
  [indsy,indsx]=find(ori==1);
  
 
 datacx=indcx-(apartsize+1);
 datacy=[y3(indcy)];
 datasx=indsx-(apartsize+1);
 datasy=[y3(indsy)];

 
 save(tmpREFcxout,'datacx');
 save(tmpREFcyout,'datacy');
 save(tmpREFsxout,'datasx');
 save(tmpREFsyout,'datasy');
 
 
%  
%  indexnum=numel(indy);
%  
%  for i=1:indexnum
% %      i 
% %      datax(i)
% %      i 
% %      datay(i)
%  end
%  
%  
 
