minputNF=[matPATH '/CART/inputNF' PUTDAT num2str(m,'%02i') '.mat']; 
load(minputNF,'inputNF')
ro1=401*401;
ROWREF = reshape(inputNF(:,:,1),ro1,1);
ROWSCR = reshape(inputNF(:,:,2),ro1,1); 
% % ROWSCR=beta
ROWDIF = reshape(inputNF(:,:,3),ro1,1);
ROWRHO = reshape(inputNF(:,:,4),ro1,1);


ROWSVEL = reshape(inputNF(:,:,5),ro1,1);
ROWSPHI = reshape(inputNF(:,:,6),ro1,1);

% ROWSVEL = reshape(GRADITP(:,:,2),ro1,1);
% ROWSPHI = reshape(GRADITP(:,:,4),ro1,1);


Vinput=[ ROWSCR ROWREF ROWRHO ROWDIF ROWSVEL ROWSPHI];
[inputrow inputnum]=size(Vinput);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % %              
% inputnum is num of variable excluding gust t/f value
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %               
% gstt_f is the value 1 for true(GF) % % % % % % % % % % %  
% gstt_f is the value 0 for false(non_GF) % % % % % % % % % 
% gstt_f is located last row of the variables row (+1) % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

gstt_f=inputnum+1;
sortnum=inputnum+2;
              
if putfisin==0;
ROWINPUT=Vinput;
end

if putfisin==1;

    mhandpick=[ matPATH '/HANDPICK/handpick' PUTDAT num2str(m,'%02i') '.mat']; 
%     mevalbox=[ matPATH '/EVAL/newevalbox' PUTDAT num2str(m,'%02i') '.mat']; 
            load(mhandpick,'handpick');
                    
%         load(mevalbox,'evalbox');   
% %     load(mhandpick,'handpick');
% % handpick=evalbo
    ROWGST = reshape(handpick,ro1,1);    
    ROWINPUT=[ Vinput ROWGST ];
    
    [numROW,numIN]=size(ROWINPUT);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%making other_area(the points are false) first%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
    indnan=zeros(numROW,1)+1;
    for i=1:inputnum
        indROW=(~isnan(ROWINPUT(:,i)));
        for j=1:numROW
            indnan(j)=indROW(j)*indnan(j);
        end
        clear indROW;
    end
    
    gustin=[];
    nongustin=[];
    ii1=1; ii2=1;
    for iii=1:numROW
        if (indnan(iii)==1 & ROWINPUT(iii,gstt_f)==1)
        gustin(ii1,:)=ROWINPUT(iii,:);
        ii1=ii1+1;
        end
    end
 
    for iii=1:numROW
        if  (indnan(iii)==1 & ROWINPUT(iii,gstt_f)==0)
        nongustin(ii2,:)=ROWINPUT(iii,:);
        ii2=ii2+1;
        end
    end
  
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % %  
    % % % sorting true_training and false training % % % 
    % % % to grab the points randomly % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % %


%     pregtsort=sortrows(pregtIN2,sortnum);
%     preothsort=sortrows(preothIN2,sortnum);
%     pregtIN1=squeeze(pregtsort(:,1:gstt_f));
%     preothIN1=squeeze(preothsort(:,1:gstt_f));
end


 
