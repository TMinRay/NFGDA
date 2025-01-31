clear all;
close all;
run ../hisnameys/header


tmpDTDcxout=['./tmpmat/tmpDTDcdatax.mat'];
tmpDTDcyout=['./tmpmat/tmpDTDcdatay.mat'];
tmpDTDsxout=['./tmpmat/tmpDTDsdatax.mat'];
tmpDTDsyout=['./tmpmat/tmpDTDsdatay.mat'];
tmppout=['./tmppng/tmpDTD.png'];


   
linesize=refcentersize;
apartsize=refapart;

gsize=linesize*2+1;
orip=apartsize*2+1;
ori=zeros(gsize+1,orip+1);

ori(1:gsize,apartsize+1)=2;
ori(3:6,apartsize+2)=2;
ori(3:6,apartsize)=2;

ori(12:15,apartsize+2)=2;
ori(12:15,apartsize)=2;

dists=5;
ori(4:14,apartsize+dists)=1;
ori(4:14,apartsize-(dists-2))=1;

ori(15:17,apartsize+(dists+1))=1;
ori(15:17,apartsize-(dists-1))=1;

ori(1:3,apartsize+(dists+1))=1;
ori(1:3,apartsize-(dists-1))=1;

% ori(2*(1:linesize),apartsize-(apartsize-1))=1;
% ori(2*(1:linesize),apartsize*2+1)=1;
% ori(linesize+1,apartsize-(apartsize-1))=1;
% ori(linesize+1,apartsize*2+1)=1;


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

 
 save(tmpDTDcxout,'datacx');
 save(tmpDTDcyout,'datacy');
 save(tmpDTDsxout,'datasx');
 save(tmpDTDsyout,'datasy');
 
 
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
 
