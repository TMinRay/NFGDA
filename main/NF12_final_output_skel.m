% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % Yunsung Hwang 2015 
% % % % % % % step 12: curve fitting: final output
% % 2nd order polynomail curve fitting only in x and y direction
% % it is assumed to represent GF enough but there are various methods to
% represent GF as final form
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
NF00_header;

for cindex=1:numel(ttable(:,1));
% for cindex=1:9;

    PUTDAT=ttable(cindex,:);
    startm=startt(cindex)+1;
    endm=endt(cindex);
  
    for m=startm:endm
    mindxy=[matPATH '/OUT/final' PUTDAT num2str(m,'%02i') '.mat'];
    
    mGSTout=[matPATH '/NFRESULT/nf' PUTDAT num2str(m,'%02i') '.mat'];
    load(mGSTout,'hh','hGST','smoothedhGST','skel_nfout');
        
	mppout=[matPATH '/POSTPROCESS/pp' PUTDAT num2str(m,'%02i') '.mat'];
    load(mppout, 'ppresult');
    locf=ppresult;
    
    indDTD3=locf>=0.6;
    se = strel('ball',5,1);
    pskel_nfout = double(imdilate(double(indDTD3),se)>1); 
    pskel=pskel_nfout.*hh;
    skel_nfout = bwmorph(double(pskel), 'skel', inf);
    skel_nfout2 = bwareaopen(skel_nfout,10);
    save(mindxy,'skel_nfout','skel_nfout2');  

    figure(m)
    pcolor(xi2,yi2,double(skel_nfout2))
    shading flat
    xlim([-100 100])
    ylim([-100 100])
    title('final output')
    end
end
  


