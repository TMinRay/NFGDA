% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % Yunsung Hwang 2015 
% % % % % % % step 13: obtaining evaluation boxes
% % now evaluating the performace of newly obtained ANFIS
% % evaluation box is obtained as
% % 1. get the local maxima in the hand picked region
% % 2. based on the local maxima, the region 2.5 km enlarged from that pixel
% either x or y direction is considered as the evaluation box
% % as the final outputs are only considered in x or y axis, the evaluation
% boxes are also considered only x or y direction
% The width of the evaluation box is 5 km and it is consistant with
% previous studies
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

NF00_header;
se = strel('ball',5,1);

for cindex=1:numel(ttable(:,1));

    PUTDAT=ttable(cindex,:);
    startm=startt(cindex)+1;
    endm=endt(cindex);
  
    for m=startm:endm
     
        mcartout=[matPATH '/CART/cart' PUTDAT num2str(m,'%02i') '.mat']; 
        load(mcartout, 'PARITP'); 
        REF=PARITP(:,:,1); 
        for iii=1:401
            for jjj=1:401
                if isnan(REF(iii,jjj))==1
                    REF(iii,jjj)=0;
                end
            end
        end

        mhandpick=[ matPATH '/HANDPICK/handpick' PUTDAT num2str(m,'%02i') '.mat']; 
        load(mhandpick,'handpick');
        pevalbox = bwmorph(double(handpick), 'skel', inf);
        evalbox = double(imdilate(double(pevalbox),se)>1);
        mevalbox=[ matPATH '/EVAL/newevalbox' PUTDAT num2str(m,'%02i') '.mat']; 
        
%         figure(m)
%         pcolor(xi2,yi2,evalbox)
%         shading flat;
%         hold on;
%         colorbar
%         contour(xi2,yi2,evalbox,'y');
%         title('REF + evalbox')
        
        save(mevalbox,'evalbox');   
    end
end

