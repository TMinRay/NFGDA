clear all;
close all;
run ../hisnameys/header


tmpSMTHcxout=['./tmpmat/tmpSMTHcdatax.mat'];
tmpSMTHcyout=['./tmpmat/tmpSMTHcdatay.mat'];
tmpSMTHsxout=['./tmpmat/tmpSMTHsdatax.mat'];
tmpSMTHsyout=['./tmpmat/tmpSMTHsdatay.mat'];
tmppout=['./tmppng/tmpSMTH.png'];


   
linesize=refcentersize+1;
apartsize=refapart;

gsize=linesize*2+1;
orip=apartsize*2+1;
ori=zeros(gsize+1,orip+1);

ori(1:gsize,apartsize+1)=2;
ori(3:7,apartsize+2)=2;
ori(3:7,apartsize)=2;

ori(13:17,apartsize+2)=2;
ori(13:17,apartsize)=2;

dists=5;
ori(5:15,apartsize+dists)=1;
ori(5:15,apartsize-(dists-2))=1;

ori(16:19,apartsize+(dists+1))=1;
ori(16:19,apartsize-(dists-1))=1;

ori(1:4,apartsize+(dists+1))=1;
ori(1:4,apartsize-(dists-1))=1;

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

 
 save(tmpSMTHcxout,'datacx');
 save(tmpSMTHcyout,'datacy');
 save(tmpSMTHsxout,'datasx');
 save(tmpSMTHsyout,'datasy');
 
 
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
 
