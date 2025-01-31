% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % Yunsung Hwang 2015 
% % % % % % % step5: combining line features estimated using ../IMG/ to obtain 
% beta as final results
% % % % % % % % % % % % % % % % % % % % % % % %

NF00_header;
% % % % % go to ../IMG/
% % % % % take a look at what exe_3_img_exe.m
% % % % % see detailed sub-routines in exe_3_img_exe.m
% % % % % FTC --> template, scoring functions are defined there
% % % % % it will take time
% % % % % take a look at Y.Hwang's Master's thesis
% % % % % This will take a huge amount of time 
% % % % % and this is where you can improve the algorithm
% % % % % rotation and conversion (polar-Cart) will take time
% % % % % hint: apply number of templates intead of rotating/converting
% % % % % Good luck!!

run ../IMG/exe_3_img_exe.m;

for cindex=1:numel(ttable(:,1));
    
    PUTDAT=[ ttable(cindex,:) ];
    startm=startt(cindex)+1;
    endm=endt(cindex);   
    LINE=[];
    linefeat1=[];
    linefeat2=[]; 
    for m=startm:endm
        mTOTLINout=[matPATH '/LINE/beta' PUTDAT num2str(m,'%02i') '.mat'];
        mLINEout=[matPATH '/LINE/Z/linez' PUTDAT num2str(m,'%02i') '.mat'];    
        load(mLINEout, 'linez');
        mLINEoutd=[matPATH '/LINE/DELZ/linedelz' PUTDAT num2str(m,'%02i') '.mat'];
        load(mLINEoutd, 'linedelz');
        mLINEout2=[matPATH '/LINE/CELLZ/widecellz' PUTDAT num2str(m,'%02i') '.mat'];
        load(mLINEout2, 'widecellz');
        mcartout=[matPATH '/CART/cart' PUTDAT num2str(m,'%02i') '.mat']; 
        load(mcartout, 'PARITP'); 
        nanind=isnan(PARITP(:,:,1));
        pbeta=(linez+linedelz)./2;
        pbeta(nanind==1)=nan;
        beta=pbeta-widecellz; 
        beta(beta<0)=0;
        save(mTOTLINout,'beta');   
        
%         figure(m)
%         pcolor(xi2,yi2,double(pbeta))
%         shading flat
        
        
        clear linez linedelz widecellz pbeta beta;
   end
end





   