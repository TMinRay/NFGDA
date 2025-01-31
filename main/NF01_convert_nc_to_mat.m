% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % Yunsung Hwang 2015 
% % % % % % % step 1: convert mat from netcdf
% % % % % % % file names can be changed
% % % % % % % netcdf files are in ../netcdf/KOUN100611/
% % % % % % % mat files will be written in ../mat/polar/
% % % % % % % ttable contains info about "RADAR,YY,MM,DD
% % % % % % % ttable example
% % % % % % % ttable=char( 'KOUN100611', ...
% % % % % % %              'KOUN100707', ...
% % % % % % %              'KOUN100712', ...
% % % % % % %              'KOUN120403', ...
% % % % % % %              'KAMX120927', ...
% % % % % % %              'KSJT130606', ...
% % % % % % %              'KTLX130601', ...
% % % % % % %              'KMQT120609', ... 
% % % % % % %              'KTLX130309', ...
% % % % % % %              'KDDC120330', ...
% % % % % % %              'KAMX120927' ...    
% % % % % % %              );
% % % % % % % % % % % % % % % % % % % % % % % % 

% clear all;
% close all;
% % % % % % or header can be modified, here only header is available for
% image processing for FTC in ../IMG


% ttable=char( 'KTLX130601', ...
%              'KMQT120609');
% startt=[1 12];
% endt=[11 18];

NF00_header;

oriPATH=['../netcdf'];
matPATH=['../mat'];

% % numel(ttable(:,1)) give the number of cases   
% % cindex=case index
% % startm=starting time if different cases starts different time,
% % they can represented in an array
% % endm=ending time

for cindex=1:numel(ttable(:,1));
    
     PATH=[ oriPATH '/' ttable(cindex,:) '/'];
     PUTDAT=ttable(cindex,:);
%      startm=1;
%      endm=2;

    startm=startt(cindex);
    endm=endt(cindex);

     PARROT=[];
     INDROT=[];     

WIDPATH=[PATH 'SpectrumWidth/00.50/'];
REFPATH=[PATH 'Reflectivity/00.50/'];
VEL1PATH=[PATH 'Velocity/00.50/'];
RHOPATH=[PATH 'RhoHV/00.50/'];
DIFPATH=[PATH 'Differential_Reflectivity/00.50/'];
PHIPATH=[PATH 'PhiDP/00.50/'];

WIDF=[WIDPATH '*.netcdf'];
REFF=[REFPATH '*.netcdf'];
RHOF=[RHOPATH '*.netcdf'];
PHIF=[PHIPATH '*.netcdf'];
DIFF=[DIFPATH '*.netcdf'];
VEL1F=[VEL1PATH '*.netcdf'];

WIDfile=dir(WIDF);
REFfile=dir(REFF);
VEL1file=dir(VEL1F);
RHOfile=dir(RHOF);
PHIfile=dir(PHIF);
DIFfile=dir(DIFF);

for m=startm:endm
    
    t=m;
    
    mpolarout=[matPATH '/POLAR/polar' PUTDAT num2str(m,'%02i') '.mat']; 

%  load(mRAWout, 'PARROT'); 
% %  PAR(0~1) 1.REF 2.VEL 3.WID 4.PHI 
% %  5.RHO 6.DIF GRID (400*720)


WIDFULL=[WIDPATH WIDfile(m).name];
REFFULL=[REFPATH REFfile(m).name];
VEL1FULL=[VEL1PATH VEL1file(m).name];
RHOFULL=[RHOPATH RHOfile(m).name];
PHIFULL=[PHIPATH PHIfile(m).name];
DIFFULL=[DIFPATH DIFfile(m).name];

ncid=netcdf.open(WIDFULL,'nowrite');
ncid=netcdf.open(REFFULL,'nowrite');
ncid=netcdf.open(RHOFULL,'nowrite');
ncid=netcdf.open(PHIFULL,'nowrite');
ncid=netcdf.open(DIFFULL,'nowrite');
ncid=netcdf.open(VEL1FULL,'nowrite');

azR=ncread(REFFULL,'Azimuth');
azV1=ncread(VEL1FULL,'Azimuth');

Ref=ncread(REFFULL,'Reflectivity');
Width=ncread(WIDFULL,'SpectrumWidth');
Vel1=ncread(VEL1FULL,'Velocity');
Phi=ncread(PHIFULL,'PhiDP');
Rho=ncread(RHOFULL,'RhoHV');
Dif=ncread(DIFFULL,'Differential_Reflectivity');     
 netcdf.close(ncid)
[RAzimuth] = size(azR);
[VAzimuth1] = size(azV1);

[GateR,RAzimuth] = size(Ref);
[GateV,VAzimuth1] = size(Width);
[GateV,VAzimuth1] = size(Vel1);
[GateR,RAzimuth] = size(Phi);
[GateR,RAzimuth] = size(Rho);
[GateR,RAzimuth] = size(Dif);

Gate2=100*4;
SRef=Ref(1:Gate2,1:RAzimuth);
SRho=Phi(1:Gate2,1:RAzimuth);
SPhi=Rho(1:Gate2,1:RAzimuth);
SDif=Dif(1:Gate2,1:RAzimuth);
SWidth=Width(1:Gate2,1:VAzimuth1);
SVel1=Vel1(1:Gate2,1:VAzimuth1);

AZM=720;

kidR=find(azR<0.4);
kidV1=find(azV1<0.4);
% kidR

if numel(kidR(:))>1
    kidR=max(kidR);
end

if numel(kidV1(:))>1
    kidV1=max(kidV1);
end



azRh=zeros(size(azR));
azV1h=zeros(size(azV1));


for i = 1:AZM
    azRh(i)=azR(kidR+i-1);
    if kidR+i-1==720
        kidR=1-i;
    end ; end ;

for i = 1:AZM
    azV1h(i)=azV1(kidV1+i-1);
    if kidV1+i-1==720
        kidV1=1-i;
    end ; end ;
 
RidR=find(azR<0.4);
if numel(RidR(:))>1
    RidR=max(RidR);
end


RidV1=find(azV1<0.4);

% RidV1=find(azR<0.4);
if numel(RidV1(:))>1
    RidV1=max(RidV1);
end



for i = 1:AZM
PARROT(1:Gate2,i,1)=SRef(1:Gate2,RidR+i-1);
PARROT(1:Gate2,i,4)=SPhi(1:Gate2,RidR+i-1);
PARROT(1:Gate2,i,5)=SRho(1:Gate2,RidR+i-1);
PARROT(1:Gate2,i,6)=SDif(1:Gate2,RidR+i-1);
    if RidR+i-1==720
        RidR=1-i;
    end
end

for i = 1:AZM
PARROT(1:Gate2,i,3)=SWidth(1:Gate2,RidV1+i-1);
PARROT(1:Gate2,i,2)=SVel1(1:Gate2,RidV1+i-1);
    if RidV1+i-1==720
        RidV1=1-i;
    end
end

clear KidR KidV1 KidV2 RidR RidV1 RidV2 ;
    
INDROT(:,:,1)=(PARROT(:,:,1)>-999 & PARROT(:,:,1)<100);   
INDROT(:,:,2)=(PARROT(:,:,2)>-100 & PARROT(:,:,2)<100);
INDROT(:,:,6)=(PARROT(:,:,6)>-5 & PARROT(:,:,6)<8);
INDROT(:,:,3)=(PARROT(:,:,3)<999 &  PARROT(:,:,3)>0);
INDROT(:,:,4)=(PARROT(:,:,4)<=360 & PARROT(:,:,4) >=0);
INDROT(:,:,5)=(PARROT(:,:,5)<=1 & PARROT(:,:,5) >=0);

idx=6;
 for kk=1:idx
    for gg=1:Gate2
        for aa=1:AZM
        if( INDROT(gg,aa,kk) ==0 )
        PARROT(gg,aa,kk)=nan;
        end
        end
    end
 end 
 save(mpolarout, 'PARROT');
 clear  PARROT INDROT;
end %%%%%%%% for m==startm:endm
end %%%%%%%% for cindex==1:5


