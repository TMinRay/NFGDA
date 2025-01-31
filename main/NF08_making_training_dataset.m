% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % Yunsung Hwang 2015 
% % % % % % % step 8: making a training dataset based on hand-picked region
% % % % % % % hand-picked region will be converted 0 and 1 and represent
% ground truth
% The region for the ground truth can be changed to that of "evaluation
% box"
% % but for this example, hand-picked region is assumed to be more trustful
% % only several representative cases will be used here and those cases
% should be eliminated while evaluating the performance
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Again, the cases used in making training data should not be evaluated for
% the performance since it will erroneously boost up the performace
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % training date and time should be determined in header file
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

NF00_header;
xsum=0;
totgstIN=[];
totothIN=[];
for cindex=trncasest:trncasend;
      
    PUTDAT=[ ttable(cindex,:)];
    startm=trnstartt(cindex);
    endm=trnendt(cindex);   
    for m=startm:endm
    % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % put fis in is for making training for FIS training %
    % putfisin=1 for "on" fpufisin=0 "off"               %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    putfisin=1;
    NF08_reading_files_FIS 
    [x1 x2]=size(gustin);
    xsum=xsum+x1;
               
    dgstIN=double(gustin);
    dothIN=double(nongustin);
    totgstIN=[totgstIN; dgstIN];
    totothIN=[totothIN; dothIN];
    
    clear ROWREF ROWstdREF ROWstdVEL ROWWID ROWGST ...
    ROWRHO ROWDIF ROWPHI dgstIN dothIN ROWINPUT preIN1;
    end
end


    [totgstROW totgstCOL]=size(totgstIN);
    [totothROW totothCOL]=size(totothIN);
    
    othrand=rand(totothROW,1);
    totothIN2=[totothIN, othrand];
    clear totothIN;
    preothsort=sortrows(totothIN2,totothCOL+1);
    totothIN=squeeze(preothsort(:,1:totothCOL));
     
    xsum;
    interxsum=round(xsum/6);
     
     
    trothnumst=1;
    trothnumen=interxsum*16;
    ttothnumst=trothnumen+1;
    ttothnumen=ttothnumst+interxsum;
    tvothnumst=ttothnumen+1;
    tvothnumen=tvothnumst+interxsum-5;
               
    trothIN=totothIN(trothnumst:trothnumen,:);
    ttothIN=totothIN(ttothnumst:ttothnumen,:);
    tvothIN=totothIN(tvothnumst:tvothnumen,:);
     
    tottr=[totgstIN; trothIN;];
    
    tottrout=['./NF09_new_training.dat'];
    save( tottrout, 'tottr' ,'-ascii');
