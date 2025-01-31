% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % Yunsung Hwang 2015 
% % % % % % % step 7: hand-picking possible regions of GF
% % % % % % % hand-picked region will be used to estimate evaluation boxes
% and histograms can be obtained based on those evaluation boxes
% %  based on preprocesses 6 variables
% %  PAR(0~1) 1.REF 2.BETA 3. DIF 4. RHO 5. SDV 6. SDP
% % possible region will be hand picked
% % rough thresholds can be applied.
% % 1. clicking when Z is ploted using (pcolor)
% % 2. click a certain point where it can be part of the GF
% % 3. make a possible closed shape linking the points clicked
% % - the clicked point are represented as "red o"
% % 4. after the clicking is done, type "y" hit enter to move to the next
% time step or different cases

NF00_header;

 
for cindex=1:numel(ttable(:,1));
    
    

    PUTDAT=ttable(cindex,:);
    startm=startt(cindex)+1;
    endm=endt(cindex);
    
    
     for m=startm:endm
         
        
    numpt=1000;
    minputNF=[matPATH '/CART/inputNF' PUTDAT num2str(m,'%02i') '.mat']; 
    load(minputNF,'inputNF')
%     FREF=double(inputNF(:,:,1)>0 & inputNF(:,:,1)<50).*1;  
%     FRHO=double(inputNF(:,:,5)<1.0).*1;  
%     FDIF=double(inputNF(:,:,6)>1).*1;  
%     FLINE=double(inputNF(:,:,8)>0).*1;
%     FSTDV=double(inputNF(:,:,7)>0).*1;
%     FSTDP=double(inputNF(:,:,9)>10).*1;
%     indh=FREF+FRHO+FDIF+FLINE+FSTDV+FSTDP;
    figure(m);
    figdiagam = figure(m);
    set(figdiagam,'Position',[10 200 700 700]);
    set(figdiagam, 'PaperPositionMode','auto') 
    pcolor(xi2,yi2,inputNF(:,:,1));
    hold on;
    caxis([5 30])
% % %     caxis are focused on Z close to GF values
    shading flat;   

    handpick = false(size(xi2));
    status = 'y';

    while status == 'y'
    'Press "enter" to click the points for the gust'
    pause;
    'take a look at a figure and resize the figure to start Rpick'
    'after clicking type "y" at asking to stop'
    for p = 1:numpt
        count = 0;
        gustin=[];
        [x1,y1] = ginput(1);
        xloc(p)  = x1 ;
        yloc(p)  = y1;
        text(x1,y1,['o'],'Color','r', 'Fontsize',6);     
        reply = input('Do you want stop? Y/N [Y]: ', 's');
        if reply=='y'    
        break;
        end
    end
    'p='
    p
    xloc2=xloc(1:p)';
    yloc2=yloc(1:p)';
    xax2 = [xloc2 ; xloc2(1)]; yax2 = [yloc2 ; yloc2(1)];
    [curr_handpick]= inpolygon(xi2,yi2,xax2,yax2);
    contour(xi2,yi2,double(curr_handpick),'Color','y','Linewidth',1);

    handpick = curr_handpick | handpick;
    status = input('Do you want to label another GF? Y/N: ','s');
    end

    mhandpick=[ matPATH '/HANDPICK/handpick' PUTDAT num2str(m,'%02i') '.mat']; 
    save(mhandpick,'handpick');
    clear xloc yloc xloc2 yloc2 xax2 yax2 x1 y1 p FREF;
    close(m);
%     clf
    end 
end