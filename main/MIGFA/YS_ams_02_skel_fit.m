% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % Yunsung Hwang 2015 
% % % % % % % step 12: curve fitting: final output
% % 2nd order polynomail curve fitting only in x and y direction
% % it is assumed to represent GF enough but there are various methods to
% represent GF as final form
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
run '../NF00_header';

for cindex=1:numel(ttable(:,1));
% for cindex=7:9;

    PUTDAT=ttable(cindex,:);
    startm=startt(cindex)+1;
    endm=endt(cindex);
  
    for m=startm:endm
        
    mindxy=['../../mat/OUT/MIGFAfinal' PUTDAT num2str(m,'%02i') '.mat'];

    Mresultf2=[ '../../mat/MIGFAOUT/cart' PUTDAT  ...
    num2str(m,'%02i')  '.mat'];

    load(Mresultf2, 'migfa'); 

%     locf=ppresult;
%     h = fspecial('gaussian', 11,1);
%     smoothedlocf = imfilter(locf,h,'replicate');
    indDTD3=migfa>0;
    se = strel('ball',3,0);
    pskel_nfout = double(imdilate(double(indDTD3),se)>0); 
    pskel=pskel_nfout;
    skel_nfout = bwmorph(double(pskel), 'thin', inf);
% % %         mindxy=[matPATH '/OUT/newfinal' PUTDAT num2str(m,'%02i') '.mat'];
        save(mindxy,'skel_nfout');  
%     figure(m)
% %     plot(xi2(logical(skel_nfout)),yi2(logical(skel_nfout)),'r.')
%     pcolor(xi4,yi4,double(skel_nfout))
%     shading flat
%     xlim([-70 70])
%     ylim([-70 70])
%     title('final output')
    end
end
  


