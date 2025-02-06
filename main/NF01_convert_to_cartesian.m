% % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % Yunsung Hwang 2015 
% % % % % % % % step 2: convert polar to cartesian
% % % % % % % % file names can be changed
% % % % % % % % (401,401,6) sized polar files are in ../mat/polar/
% % % % % % % % mat files will be written in ../mat/cart/
% % % % % % % % ttable contains info about "RADAR,YY,MM,DD
% % % % % % % % % % % % % % % % % % % % % % % % %
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         
% % % % % % % this is an ex version 
% % % % % % % after converting NCDC to netcdf (from V06 to netcdf)
% % % % % % % the folder /POLAR/ should contain files
% % % % % % % also "header.m" should be changed!!!!
% % %  PAR(0~1) 1.REF 2.VEL 3.WID 4.PHI 5.RHO 6.DIF 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

NF00_header;
%%%% NF01_convert_to_cartesian
for cindex=1:numel(ttable(:,1));

    PUTDAT=ttable(cindex,:);
    startm=startt(cindex);
    endm=endt(cindex);
    PARITP=[];
    for m=startm:endm
        mpolarout=[matPATH '/POLAR/polar' PUTDAT num2str(m,'%02i') '.mat'];
        load(mpolarout, 'PARROT'); 
        mcartout=[matPATH '/CART/cart' PUTDAT num2str(m,'%02i') '.mat'];
        PARITP=zeros(401,401,6);
        varnum=[1 2 3 4 5 6];
        sdphi=zeros(400,720);
        phi=double(PARROT(:,:,4));
        phi(phi<0)=nan;
        phi(phi>360)=nan;
        buf =zeros(400,720,5);
        % for displaceR=-2:2
        %     buf(3:end-2,:,displaceR+3)=phi(3+displaceR:end-2+displaceR,:);
        % end
        % sdphi(3:end-2,:) = std(buf(3:end-2,:,:),0,3,"omitmissing");
        for displaceR=-2:2
            buf(5:Gate2-2,:,displaceR+3)=phi(5+displaceR:Gate2-2+displaceR,:);
        end
        sdphi(5:Gate2-2,:) = std(buf(5:Gate2-2,:,:),0,3,"omitmissing");
        PARROT(:,:,4)=sdphi;
% % % % %         new way to obtain SD(phi_DP) only in radial direction
        for i=1:6
            a=double(PARROT(:,:,i));
            PARITP(:,:,i) = griddata(x,y,a,xi2,yi2);
        end
        save(mcartout, 'PARITP');
        clear PARITP PARITP;
    end
end
