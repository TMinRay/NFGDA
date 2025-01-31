% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % Yunsung Hwang 2015 
% % % % % % % step 11: post-processing; moving averaging
% % % % % % % 1. link fragmented lines + 2. further removing false alarms
% % the purpose of this process is additional removing remaing false alarms
% % simple moving averaging with 16 different templates (11.25 deg)
% % It is assumed that templates represent 11.25 rotation properly
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % should have subfolder "movingavg" since information of different
% templates are in that folder
% never delete the folder "movingavg" unless a new post-processing is
% proposed
clear all
close all
NF00_header
  
for cindex=1:numel(ttable(:,1));
    PUTDAT=ttable(cindex,:);
    startm=startt(cindex)+1;
    endm=endt(cindex);
   
   for m=startm:endm
%       for m=startm:startm    
        mGSTout=[matPATH '/NFRESULT/nf' PUTDAT num2str(m,'%02i') '.mat'];
        load(mGSTout,'hh','hGST','smoothedhGST','skel_nfout');
        mppout=[matPATH '/POSTPROCESS/pp' PUTDAT num2str(m,'%02i') '.mat'];
         pppout=['./check_png/smootheda2' PUTDAT num2str(m,'%02i') '.png'];

%         
%         ppresult=zeros(401,401);
% %         a2=skel_nfout;
        se = strel('ball',10,1);
        pctlth=prctile(reshape(smoothedhGST,401*401,1),80);
        pskel_nfout = double(imdilate(double(smoothedhGST>pctlth),se)-1);
        pctlth2=prctile(reshape(smoothedhGST,401*401,1),70);
        
%         skel_nfout2=medfilt2(skel_nfout,[3 3]);
%         pctlth=prctile(reshape(smoothedhGST,401*401,1),75)
%         pskel_nfout = double(imdilate(double(smoothedhGST>0.5),se));  
        
%         figure(m)
%         pcolor(xi2,yi2,double(pskel_nfout)-1)
%         shading flat;
%         xlim([-100 100])
%         ylim([-100 100])
%         colorbar
%         title('Postprocessing_2 moving average (this will be converted to 0 and 1s)')

%         h = fspecial('gaussian', 11,1);


% % % % % % % % % % % % % % % % % % % %           
% % % % % % % % % % % % % % % % % % %           
% % % % % % % % % % % % % % % % % % % %           
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %           
% % % % % % % % % % % % % % % % % % %           
% % % % % % % % % % % % % % % % % % % %           
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %           
% % % % % % % % % % % % % % % % % % %           
% % % % % % % % % % % % % % % % % % % %           
% % % % % % % % % % % % % % % % % % % % 
% 
         ppresult=zeros(401,401);
        a2=hGST;
        inda2=(a2>0);
          
% % % % % % % %         smootheda2 = imfilter(double(smoothedhGST>0.3),h,'replicate');
        numINT =8;
        ppresult(:,1:numINT)=0;
        ppresult(:,401-numINT+1:401)=0;
        ppresult(1:numINT,:)=0;
        ppresult(401-numINT+1:401,:)=0;
        longsize=21;        
        for ii=1+numINT:401-numINT
            for jj=1+numINT:401-numINT  
                if inda2(ii,jj)>0;
                meannum8=zeros(1,17);
                for i=1:17
                    cxout=['./movingavg/ccx' num2str(i,'%02i') '.mat'];
                    cyout=['./movingavg/ccy' num2str(i,'%02i') '.mat'];
                    load(cxout,'ccx');
                    load(cyout,'ccy');
                    indexcnum=numel(ccx);
                    cbox=[];         
                        for kind=1:indexcnum
                            cadd=[a2(ii+ccx(kind),jj+ccy(kind))];
                            cbox=[cbox; cadd;];
                            clear cadd;
                        end
                        indcb=((cbox)>0);
                        numtcb=sum(indcb(:));
                        cbr=numtcb/(indexcnum);
                        if (cbr>0.1)
                        meannum8(i)=mean(cbox(:),"omitmissing");
                        clear cbox sbox totbox;
                        else
                        meannum8(i)=0;
                        end
                end                
                ppresult(ii,jj)=max(meannum8(:));
                clear meannum8
                else
                ppresult(ii,jj)=0;    
                end
            end 
        end
        save(mppout, 'ppresult');
       
        
        figure(cindex)
        pcolor(xi2,yi2,double(ppresult))
%         pcolor(xi2,yi2,double(inda2))
        
        shading flat;
        xlim([-100 100])
        ylim([-100 100])
        colorbar
        caxis([0 1])
        title('Postprocessing_2 moving average (this will be converted to 0 and 1s)')
        saveas(gcf,pppout,'png')
        close(cindex)
        clear ppresult hGST2; 
   end
end
