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

% clc

oriPATH=['../necdf'];
matPATH=['../mat'];

% ttable=char( 'KABR20140621_21');
ttable=char( 'KGGW20210722_23');
         
         
         
%   1       ttable=char( 'KOUN100611', ...
%   2           'KOUN100707', ...
%   3           'KOUN100712', ...             
%   4           'KAMX120927', ...
%   5           'KSJT130606', ...
%   6           'KTLX130601', ...
%   7           'KMQT120609', ... 
%   8           'KTLX130309', ...
%   9           'KDDC120330', ...
%   10           'KABR140621', ...
% % % %   11           'KDMX140630', ...
%   12           'KILX140621', ...             
%   13           'KLSX140621');
% startt=[ 2 12 2 28 12  9  26 33 22   1 11  1 ];
%   endt=[11 18 10 39 21 15  47 45 39 10 20 10 ];

startt=[1];
endt=[2];


% % % % % this part is for choosing specific date for making training
% dataset
% % set for training dataset
% trncasest=the first case to be used as training data
% trncasend=the last case to be used as training data
% trncasest=2;
% trncasend=4;
% trnstartt=[5 16 7 30 20 15 11 ];
% trnendt=[6 17 8 31 21 16 12];
% this is extremely important to train membership functions
% be careful for selecting training dataset!!

trncasest=1;
trncasend=1;
trnstartt=[2];
trnendt=[2];

         
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
% % % Uncomment this and take a look at MIGFA folder 
% % % You might contact Darrel M. Kingfield in order to extract information
% of MIGFA(operational) from Level II data
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % Otherwise, you might think of comparing dataset,
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% This Mlath, Mlonh is for converting MIGFA results to 0.5 by 0.5 km
% Cartesina Grid, the Mlat, Mlon represents the center(location) of radar
% in United States, those should be provided case by case
% For example cases, lats and lons are provided as follows
% % Mlath=[-97.4708 -97.4708 -97.4708 -80.4131 ...
% %     -100.4925 -97.2775 -98.413 -89.337 -90.683 -92.4274 -97.2775 -80.0344 -80.4131 ];
% % Mlonh=[35.2435 35.2435 35.2435 25.6106 ...
% %     31.3714 35.3331 45.456 40.151 38.6991 46.5339 35.3331 37.7631 25.6106 ];
% % This is the location of "KABR" only for example
Mlonh=[ -98.4224444 ];
Mlath=[45.4483056 ];





% % % ccase2=[1 2 3 5 6 7 8 9 10 11]; 
% % % ccasenum2=numel(ccase2);
% % % 
% % % spcase=[2 3 5 6 11]; 
% % % spcasenum=numel(spcase);
% % % 
% % % casenumnum=numel(casenumtotal);
% % % casenumnum2=numel(casenumtotal2);
% % % stcasenum=1;
% % % ndcasenum=10;
% % % 
% % % % startt=[1 12 2 1 28 12 9 1 11 1 28 ];
% % % % endt=[11 18 10 12 40 22 15 10 20 10 40];
% % % 
% % % startt=[1 12 2 1 28 12 9 1 11 1 28 ];
% % % endt=[11 18 10 12 39 21 15 10 20 10 39];
% % % 


%   1       ttable=char( 'KOUN100611', ...
%   2           'KOUN100707', ...
%   3           'KOUN100712', ...             
%   4           'KAMX120927', ...
%   5           'KSJT130606', ...
%   6           'KTLX130601', ...
%   7           'KMQT120609', ... 
%   8           'KTLX130309', ...
%   9           'KDDC120330', ...
%   10           'KABR140621', ...
% % % %   11           'KDMX140630', ...
%   12           'KILX140621', ...             
%   13           'KLSX140621');

% % % % % % % % % % % % % % % % % % % % % % % % % % % originally

% % % % ttable=char( 'KOUN100611', ...
% % % %              'KOUN100707', ...
% % % %              'KOUN100712', ...
% % % %              'KOUN120403', ...
% % % %              'KAMX120927', ...
% % % %              'KSJT130606', ...
% % % %              'KTLX130601', ...
% % % %              'KMQT120609', ... 
% % % %              'KTLX130309', ...
% % % %              'KDDC120330', ...
% % % %              'KAMX120927' ...    
% % % %              );



   







