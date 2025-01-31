% tic
clear all;
close all;
NF00_header;
totnumag=0;
totnumao=0;
gf=cell(1,7);
st=gf;
ca=gf;

for cindex=1:numel(ttable(:,1));

    PUTDAT=ttable(cindex,:);
    startm=startt(cindex)+1;
    endm=endt(cindex);

     
for m=startm:endm
    t=m;
    tot=zeros(401,401,6);
        minputNF=[matPATH '/CART/inputNF' PUTDAT num2str(m,'%02i') '.mat']; 
        load(minputNF,'inputNF')
                nokink=inputNF(:,:,3);
        nokink(nokink>=8)=nan;
        inputNF(:,:,3)=nokink;
        clear nokink;
        nokink=inputNF(:,:,5);
        nokink(nokink==0)=nan;
        inputNF(:,:,5)=nokink;
        clear nokink;
        tot=inputNF;
        clear inputNF;
        mevalbox=[ matPATH '/EVAL/newevalbox' PUTDAT num2str(m,'%02i') '.mat']; 
        load(mevalbox,'evalbox');   
      
      gustin=evalbox;
      ist=zeros(401,401);
      ica=ist;
      oth=logical(gustin.*(-1)+1);
      ref=tot(:,:,1);
      for ii=1:401
          for ij=1:401
              if double(oth(ii,ij))>0 && ref(ii,ij)>=25
                  ist(ii,ij)=1;
              end
              if double(oth(ii,ij))>0 && ref(ii,ij)<25
                  ica(ii,ij)=1;
              end
          end
      end
      gst=logical(gustin);
      ist=logical(ist);
      ica=logical(ica);
      clear gustin;
for k=1:6

      aa=tot(:,:,k);
      agst=aa(gst);
      ast=aa(ist);
      aca=aa(ica);
      
      agnidx=~isnan(agst);
      if k==2
      agnidx2=(agst~=0);  
      for kkk=1:numel(agnidx2)
      agnidx(kkk)=agnidx(kkk)*agnidx2(kkk);
          
      end
     
      end
      
      astidx=~isnan(ast);
      acaidx=~isnan(aca);
      
%       agst(agnidx)
      gf{k}=[gf{k}; agst(agnidx);];
      st{k}=[st{k}; ast(astidx);];
      ca{k}=[ca{k}; aca(acaidx);];
      clear aa agst ast aca;
end
    clear tot;
%             toc

end
% save(mANALGSTout,'gf','st');
end
mANALGSTout=[matPATH '/HIST/stat3regions.mat'];
save(mANALGSTout,'gf','st','ca');
