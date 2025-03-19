clear all;
close all;

AZM=720;
Gate2=100*4;

   r=(0:400-1)*0.25;  
   rotaz=0:0.5:359.5;
   for j=1:Gate2
   for k=1:AZM
   x(j,k)=r(j)*sin(rotaz(k)*pi/180);
   y(j,k)=r(j)*cos(rotaz(k)*pi/180);
   end
   end
   
   [xi2,yi2] = meshgrid(-100:0.5:100,-100:0.5:100);
  [xi4,yi4] = meshgrid(-100:0.5:100,-100:0.5:100);

matPATH=['../mat'];


case_name = 'KABX20200705_21';
% case_name = 'KABX20200715_23';

no_eval = false;

startt=1;
endt=8;

         
% % for image-processing FTC % %
% % rotdegree determine degree increase i.e., 20 deg here
% % angint is angle interval where how many angle should be shifted
% % i.e., 0.5 means super resolution so 40 radius should be rotated

rotdegree=180/9;
angint=0.5;
rotAZ=round(rotdegree/0.5);
rotnum=round(180/rotdegree);
thrREF=5;
rotbackrad=deg2rad(rotdegree);
cellcsrthresh=0.5;
idcellscrthresh=0.5;
% ref cbox threshold ref sbox threshold
thrdREF=0.3; 
cellthresh=5;
cbcellthrsh=0.8;
