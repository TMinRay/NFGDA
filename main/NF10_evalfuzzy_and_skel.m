
NF00_header;
h = fspecial('gaussian', 11,1);

for cindex=1:numel(ttable(:,1));
    

    PUTDAT=ttable(cindex,:);
    startm=startt(cindex)+1;
    endm=endt(cindex);
    

%     fuzzfile=['./new_comb.fis'];
% % % Reference for the FIS based on Y.Hwang's until 2015 study    
% % % using 6 variables of Z, beta, Z_DR, rho_hv, SD(v_r) and SD(phi_dp)
% % % can be used as comparison but new FIS should do better than the
% refernce FIS
    fuzzfile=['./NF00ref_YHWANG.fis'];

    fuzzGST=readfis(fuzzfile);
    for m=startm:endm
        % % % % % % % % % % % % % % % % % % % % % % % 0% % % % 
        % put fis in is for making training for FIS training %
        % putfisin=1 for "on" fpufisin=0 "off"               %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        putfisin=0;
        NF08_reading_files_FIS;            
        dINPUT=double(ROWINPUT);
        PREGST=evalfis(dINPUT,fuzzGST);
        hGST=reshape(PREGST,401,401);
        % % % % % % numthr below determines Y/N so the number should be
        % based on statistical results!!!!!
        hGST(hGST<=0)=nan;
         numthr=0.7;
         
        for iii=1:401
            for jjj=1:401
                if (hGST(iii,jjj)>=numthr) 
                    hGST(iii,jjj)=1;
                else
                    hGST(iii,jjj)=0;
                end
                if (hGST(iii,jjj)>1) 
                hGST(iii,jjj)=1;
                end
                if (isnan(hGST(iii,jjj))==1) 
                hGST(iii,jjj)=0;
                end
            end
        end
        
         hh=hGST;
        
        clear hGST;
        hGST=medfilt2(hh,[3 3]);
        
        smoothedhGST = imfilter(hGST,h,'replicate');
%         se = strel('ball',5,1);
%         
%         pskel_nfout = double(imdilate(double(hGST),se)>1);

        smoothedhGST(smoothedhGST<0.1)=nan;
%         
        skel_nfout = bwmorph(double(smoothedhGST>0.3), 'skel', inf);
        
        mGSTout=[matPATH '/NFRESULT/nf' PUTDAT num2str(m,'%02i') '.mat'];
        save(mGSTout,'hh','hGST','smoothedhGST','skel_nfout');


%         

%         
%         figure(m)
%         smoothedhGST(smoothedhGST<0.1)=nan;
%         pcolor(xi2,yi2,smoothedhGST)
%         shading flat
%         hold on;
%                 plot(xi2(logical(skel_nfout)),yi2(logical(skel_nfout)),'k.')
%         caxis([0 1])
%         xlim([-100 100])
%         ylim([-100 100])
%         title('right after NFsystem (before postprocessing)')
%         
% 
        figure(m+10)
        plot(xi2(logical(hGST)),yi2(logical(hGST)),'b.')
        xlim([-100 100])
        ylim([-100 100])
        title('Postprocessing_1 NFsystem (3 by 3 window median filter)')
% 
%         
%         pcolor(xi2,yi2,hGST)
%         shading flat
        
        
        clear ROWREF ROWSHR ROWWID ROWRHO ROWDIF ROWPHI dINPUT ...
        PREGST hGST ROWIN ROWindlin;
        clear sGSTout hGST;
    end 
end

