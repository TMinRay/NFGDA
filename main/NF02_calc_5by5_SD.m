% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % Yunsung Hwang 2015 
% % % % % % % step 3: calculating SD of v_r and phi_DP
% % % % % % % file names can be changed
% % % % % % % mat files will be written in ../mat/SD/
% % % % % % % ttable contains info about "RADAR,YY,MM,DD
% % % % % % % % % % % % % % % % % % % % % % % %
NF00_header;

stdCELLcxout=['./NF03_STDdatax.mat'];
stdCELLcyout=['./NF03_STDdatay.mat'];
load(stdCELLcxout,'datacx');
load(stdCELLcyout,'datacy');
indexnum=numel(datacx);
stdrsize=2;
numINT =stdrsize;
% % % % % % % stdCELLcxout and stdCELLcyout contain info of pixels of 5 by
% 5 window such as (i+1, j+1), (i+2, j+1), (i+1, j+2) ...


for cindex=1:numel(ttable(:,1));

    PUTDAT=ttable(cindex,:);
    startm=startt(cindex)+1;
    endm=endt(cindex);
    for m=startm:endm
        mcartout=[matPATH '/CART/cart' PUTDAT num2str(m,'%02i') '.mat']; 
        load(mcartout, 'PARITP'); 
        mstdoutv=[matPATH '/SD/std' PUTDAT num2str(m,'%02i') '.mat']; 
        % a=zeros(401,401,3);
        stda = zeros(401,401,3);
        % a(:,:,1)=double(PARITP(:,:,1));
        a2 = double(PARITP(:,:,2));
        % a(:,:,3)=double(PARITP(:,:,4));   
% % % % % %         1-->ref 
% % % % % %         2-->radialvelocity
% % % % % %         4-->phidp
% % % % % %  SD(phidp) is modified to be obtained while converting
% % % % % %  here only calculating SD(v_r) 
% % % % % %  to modify, you need to modify number of array a(:,:,5) or
% a(:,:,6) you can obtain texture of varialbes as many as you want
% here only obtain SD(vr) you need to also modify kk=2:2 to something else 
% to obtain textures of variables
        [col_indices,row_indices]=meshgrid(numINT+1:401-numINT);
        cridx = row_indices(:).' + datacy;
        ccidx = col_indices(:).' + datacx;
        c_indices = sub2ind(size(a2), cridx, ccidx);
        cbox=a2(c_indices);
        indcb = ~isnan(cbox);
        numtcb=sum(indcb,1);
        cbr = numtcb/indexnum;
        row_indices = row_indices(cbr>=0.3);
        col_indices = col_indices(cbr>=0.3);
        cbox = cbox(:,cbr>=0.3);
        lin_indices = sub2ind(size(a2), row_indices, col_indices);
        cboxstd = std(cbox,0,1,"omitmissing");
        buf = zeros(401,401);
        buf(lin_indices)=cboxstd;
        stda(:,:,2) = buf;
        % % totscore(:,1:numINT)=nan;
        % % totscore(:,401-numINT+1:401)=nan;
        % % totscore(1:numINT,:)=nan;
        % % totscore(401-numINT+1:401,:)=nan;
        % for kk=2:2 
        %     for ii=1+numINT:401-numINT
        %         for jj=1+numINT:401-numINT
        %         cboxv=[];
        %         for kind=1:indexnum
        %             caddv=[a(ii+datacy(kind),jj+datacx(kind),kk)];
        %             cboxv=[cboxv; caddv;];
        %         clear caddv caddp;
        %         end
        %         indcbv=(~isnan(cboxv));
        %         numtcbv=sum(indcbv(:));
        %         cbrv=numtcbv/indexnum;
        %         if (cbrv>=0.3)
        %         stda(ii,jj,kk)=std(cboxv,"omitmissing");
        %         else
        %         stda(ii,jj,kk)=0;
        %         end   
        %         end %end Azimuth
        %     end
        % end
        save(mstdoutv,'stda');
        clear stda  ; 
       end
end

