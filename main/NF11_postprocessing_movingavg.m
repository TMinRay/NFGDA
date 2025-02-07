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
se = strel('ball',10,1);
numINT =8;

disccx = zeros(17,17);
disccy = zeros(17,17);
for i=1:17
    cxout=['./movingavg/ccx' num2str(i,'%02i') '.mat'];
    cyout=['./movingavg/ccy' num2str(i,'%02i') '.mat'];
    load(cxout,'ccx');
    load(cyout,'ccy');
    disccx(i,:)=ccx;
    disccy(i,:)=ccy;
end
indexcnum=numel(ccx);
disccx = reshape(disccx,1,17,17);
disccy = reshape(disccy,1,17,17);

for cindex=1:numel(ttable(:,1));
    PUTDAT=ttable(cindex,:);
    startm=startt(cindex)+1;
    endm=endt(cindex);
   
   for m=startm:endm
        mGSTout=[matPATH '/NFRESULT/nf' PUTDAT num2str(m,'%02i') '.mat'];
        load(mGSTout,'hh','hGST','smoothedhGST','skel_nfout');
        mppout=[matPATH '/POSTPROCESS/pp' PUTDAT num2str(m,'%02i') '.mat'];
        pppout=['./check_png/smootheda2' PUTDAT num2str(m,'%02i') '.png'];
        % pctlth=prctile(reshape(smoothedhGST,401*401,1),80);
        % pskel_nfout = double(imdilate(double(smoothedhGST>pctlth),se)-1);
        % pctlth2=prctile(reshape(smoothedhGST,401*401,1),70);

        ppresult=zeros(401,401);
        a2=hGST;
        % inda2=(a2>0);
        % ppresult(:,1:numINT)=0;
        % ppresult(:,401-numINT+1:401)=0;
        % ppresult(1:numINT,:)=0;
        % ppresult(401-numINT+1:401,:)=0;

        [row_indices, col_indices] = find(a2>0);
        inINT = (row_indices>numINT) & (row_indices<=401-numINT)...
            &(col_indices>numINT) & (col_indices<=401-numINT);
        row_indices = row_indices(inINT);
        col_indices = col_indices(inINT);

        cridx = reshape(row_indices,numel(row_indices),1,1) + disccy;
        ccidx = reshape(col_indices,numel(col_indices),1,1) + disccx;
        c_indices = sub2ind(size(a2), cridx, ccidx);
        cbox=a2(c_indices);
        indcb = cbox>0;
        numtcb=sum(indcb,3);
        cbr = numtcb/indexcnum;
        validcenter = max(cbr>0.1,[],2);
        row_indices = row_indices(validcenter);
        col_indices = col_indices(validcenter);
        cbox = cbox(validcenter,:,:);
        mc = mean(cbox,3,"omitmissing");
        mc(~(cbr>0.1))=0;
        buf = zeros(401,401);
        lin_indices = sub2ind(size(a2), row_indices, col_indices);
        buf(lin_indices) = max(mc,[],2);
        ppresult(numINT+1:401-numINT,numINT+1:401-numINT) = buf(numINT+1:401-numINT,numINT+1:401-numINT);
        % for ii=1+numINT:401-numINT
        %     for jj=1+numINT:401-numINT  
        %         if inda2(ii,jj)>0;
        %             meannum8=zeros(1,17);
        %             for i=1:17
        %                 cxout=['./movingavg/ccx' num2str(i,'%02i') '.mat'];
        %                 cyout=['./movingavg/ccy' num2str(i,'%02i') '.mat'];
        %                 load(cxout,'ccx');
        %                 load(cyout,'ccy');
        %                 indexcnum=numel(ccx);
        %                 cbox=[];         
        %                 for kind=1:indexcnum
        %                     cadd=[a2(ii+ccx(kind),jj+ccy(kind))];
        %                     cbox=[cbox; cadd;];
        %                     clear cadd;
        %                 end
        %                 indcb=((cbox)>0);
        %                 numtcb=sum(indcb(:));
        %                 cbr=numtcb/(indexcnum);
        %                 if (cbr>0.1)
        %                     meannum8(i)=mean(cbox(:),"omitmissing");
        %                     clear cbox sbox totbox;
        %                 else
        %                     meannum8(i)=0;
        %                 end
        %             end
        %             ppresult(ii,jj)=max(meannum8(:));
        %             clear meannum8
        %         else
        %             ppresult(ii,jj)=0;
        %         end
        %     end 
        % end
        save(mppout, 'ppresult');
       
        
%         figure(cindex)
%         pcolor(xi2,yi2,double(ppresult))
% %         pcolor(xi2,yi2,double(inda2))
        
%         shading flat;
%         xlim([-100 100])
%         ylim([-100 100])
%         colorbar
%         caxis([0 1])
%         title('Postprocessing_2 moving average (this will be converted to 0 and 1s)')
%         saveas(gcf,pppout,'png')
%         close(cindex)
%         clear ppresult hGST2; 
   end
end
