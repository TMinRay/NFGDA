% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % Yunsung Hwang 2015 
% % % % % % % step 6: combining 6 varialbes to be used as inputs in NF
% system
% % % % % % % % % % % % % % % % % % % % % % % %
% %  previous variables: PAR(0~1) 1.REF 2.VEL 3.WID 4.PHI 5.RHO 6.DIF 
% %  are changed to 
% %  new PAR(0~1) 1.REF 2.BETA 3. DIF 4. RHO 5. SDV 6. SDP

NF00_header;

for cindex=1:numel(ttable(:,1));


    PUTDAT=ttable(cindex,:);
    startm=startt(cindex)+1;
    endm=endt(cindex);
    
    
    for m=startm:endm
    
    pPARITP=zeros(401,401,12);
    pnansum=zeros(401,401);
    
    minputNF=[matPATH '/CART/inputNF' PUTDAT num2str(m,'%02i') '.mat']; 
 
    mcartout=[matPATH '/CART/cart' PUTDAT num2str(m,'%02i') '.mat']; 
    load(mcartout, 'PARITP'); 
    mTOTLINout=[matPATH '/LINE/beta' PUTDAT num2str(m,'%02i') '.mat'];
    load(mTOTLINout,'beta');    
    mstdoutv=[matPATH '/SD/std' PUTDAT num2str(m,'%02i') '.mat']; 
    load(mstdoutv,'stda');    
    pPARITP(:,:,1)=PARITP(:,:,1);    
    pPARITP(:,:,2)=beta;
    pPARITP(:,:,3)=PARITP(:,:,6);    
    pPARITP(:,:,4)=PARITP(:,:,5);    
    pPARITP(:,:,5)=stda(:,:,2);
    pPARITP(:,:,6)=PARITP(:,:,4);
    pnan=(isnan(pPARITP));
    for ii=1:401
        nanval=[];
        for jj=1:401
            nanval=(pnan(ii,jj,:));
            pnansum(ii,jj)=sum(nanval(:));    
            clear nanval;
        end
    end
    for i=1:6
        for ki=1:401
            for kj=1:401 
                if(pnansum(ki,kj)>0)
                pPARITP(ki,kj,i)=nan;    
                end
            end
        end
    end
    inputNF=pPARITP; 
    save(minputNF,'inputNF')
    clear pPARTITP PARITP pnansum pnan PARtot;
    end
end

