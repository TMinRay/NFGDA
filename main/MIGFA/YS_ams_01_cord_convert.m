run '../NF00_header';
     [xi3, yi3]=meshgrid(-100:0.1:100,-100:0.1:100);


for cindex=1:numel(ttable(:,1));
% for cindex=7:7;
    

    PUTDAT=ttable(cindex,:);
    startm=startt(cindex)+1;
    endm=endt(cindex);
for m=startm:endm
% for m=startm:startm

    t=m;
    
    Mresultf=[ '../../mat/MIGFAOUT/new' PUTDAT  ...
    num2str(t,'%02i')  '.mat'];
    load(Mresultf,'Mresult');
    
    Mresultf2=[ '../../mat/MIGFAOUT/cart' PUTDAT  ...
    num2str(t,'%02i')  '.mat'];

    migfa=zeros(401,401);
%     migfa= griddata(xi3,yi3,Mresult,xi4,yi4);
    migfa= interp2(xi3,yi3,Mresult,xi4,yi4);
    

        
%         figure(1)
%         pcolor(xi4,yi4,migfa)
%         shading flat
%         caxis([0 1])
%         colorbar
        
        save(Mresultf2, 'migfa'); 
        clear migfa;
    end 
end 
