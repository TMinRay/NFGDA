% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % Yunsung Hwang 2015 
% % % % % % % step 4: calculating Delta Z (time difference of Z)
% % % % % % % file names can be changed
% % % % % % % mat files will be written in ../mat/DELZ/
% % % % % % % ttable contains info about "RADAR,YY,MM,DD
% % % % % % % % % % % % % % % % % % % % % % % %

NF00_header;
     
for cindex=1:numel(ttable(:,1));
    

    PUTDAT=ttable(cindex,:);
    startm=startt(cindex)+1;
    endm=endt(cindex);
 
    for m=startm:endm
        t=m;
        pre=m-1;
        dout=[matPATH '/DELZ/delz' PUTDAT num2str(m,'%02i') '.mat'];
        mrout2=[matPATH '/POLAR/polar' PUTDAT num2str(m,'%02i') '.mat']; 
        mrout1=[matPATH '/POLAR/polar' PUTDAT num2str(pre,'%02i') '.mat'];
        a1=zeros(400,720);
        a2=a1;
        load(mrout1, 'PARROT'); 
        a1(:,:)=PARROT(:,:,1);
        a1(isnan(a1))=0;
        load(mrout2, 'PARROT'); 
        a2(:,:)=PARROT(:,:,1);  
        a2(isnan(a2))=0;
        dif2=a2-a1;   
        save(dout,'dif2');

%         mrout2=[matPATH '/CART/cart' PUTDAT num2str(m,'%02i') '.mat'];
%         mrout1=[matPATH '/CART/cart' PUTDAT num2str(pre,'%02i') '.mat'];
%         a1=zeros(401,401);
%         a2=a1;
%         load(mrout1, 'PARITP'); 
%         a1(:,:)=PARITP(:,:,1);
%         load(mrout2, 'PARITP'); 
%         a2(:,:)=PARITP(:,:,1);  
%         dif2=a2-a1;   
%         save(dout,'dif2');
%         clear a dif2 ;
% figure(m)
% pcolor(x,y,double(dif2))
% shading flat;
        clear a dif2 ;
    end
end 


