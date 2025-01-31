NF00_header;

dtdout=['./NFFIG/'];


%     for cindex=1:numel(ttable(:,1));
    for cindex=1:1
        
    PUTDAT=ttable(cindex,:);
    startm=startt(cindex)+1;
    endm=endt(cindex);    
%     for m=startm:endm   
    for m=3:3   
    minputNF=[matPATH '/CART/inputNF' PUTDAT num2str(m,'%02i') '.mat']; 
    load(minputNF, 'inputNF'); 
% %  PAR(0~1) 1.REF 2.BETA 3.Z_DR 4.rho_hv 
% %  5.SD(v_r) 6.SD(phi_DP) // GRID (401*401)
         mevalbox=[ matPATH '/EVAL/newevalbox' PUTDAT num2str(m,'%02i') '.mat']; 
        
        load(mevalbox,'evalbox');   

    xx2=xi2(logical(evalbox));
    yy2=yi2(logical(evalbox));
    maxx=max([xx2;]);
    minx=min([xx2;]);
    maxy=max([yy2;]);
    miny=min([yy2;]);
    delx=abs(maxx-minx)/2;
    dely=abs(maxy-miny)/2;
    delgrd=max(delx,dely)+5;
    xgrd1=(maxx+minx)/2-delgrd;
    xgrd2=(maxx+minx)/2+delgrd;
    ygrd1=(maxy+miny)/2-delgrd;
    ygrd2=(maxy+miny)/2+delgrd;
%     xlim([xgrd1 xgrd2])
%     ylim([ygrd1 ygrd2])
    xlim([-70 70])
    ylim([-70 70])
    
    
    
    hout{1}=['{Case' num2str(cindex) ', t=' num2str(m-(startm-1)) ' \medskip, $Z$}'];
    hout{2}=['{Case' num2str(cindex) ', t=' num2str(m-(startm-1)) ' \medskip, $\beta$}'];
    hout{3}=['{Case' num2str(cindex) ', t=' num2str(m-(startm-1)) ' \medskip, $Z_{DR}$}'];
    hout{4}=['{Case' num2str(cindex) ', t=' num2str(m-(startm-1)) ' \medskip, $\rho_{hv}$}'];
    hout{5}=['{Case' num2str(cindex) ', t=' num2str(m-(startm-1)) ' \medskip, SD $(v_{r})$}'];
    hout{6}=['{Case' num2str(cindex) ', t=' num2str(m-(startm-1)) ' \medskip, SD $(\phi_{DP})$}'];
       
    REF=inputNF(:,:,1); 
    DIF=inputNF(:,:,3); 
    RHO=inputNF(:,:,4); 
    LIN=inputNF(:,:,2); 
    SDV=inputNF(:,:,5); 
    SDP=inputNF(:,:,6); 
   
     fig=figure(m);
    set(fig,'Position',[100 100 630 800]);
    ha = tight_subplot(3,2,[.04 .01],[.03 .06],[.05 .03]);
    plotout=[ dtdout '/fig3.png'];   
    
    axes(ha(1))    % similar with subplot(221)
    axis square;
    pcolor(xi2,yi2,REF);
    hold on;
    clim = [-5 65];
    cmap1 = boonlib('zmap3',35);
    colormap(cmap1);
    shading flat; 
    caxis([-5 60])
    shading flat;
    g1=colorbar;
%     ylabel('{Meridional}(km)','interpreter', 'latex','Fontsize',12);
%     ylabel('Meridional(km)');
%     set(gca,'tick','out','box','on','TickLength'  , [.01 .01], 'LineWidth'   , 1.5);
    contour(xi2,yi2,evalbox,'y-','linewidth',1); hold on;
%     title(hout{1},'interpreter', 'latex','Fontsize',13);
    set(gca,'XTickLabel',[]);
    cbfreeze(g1);
    freezeColors;
    xlim([-70 70])
    ylim([-70 70])
    

    axes(ha(2))    % similar with subplot(221)
    axis square;
    pcolor(xi2,yi2,LIN);
    caxis([0 0.5]); 
    hold on;
    shading flat;
    g1=colorbar;
    cmap2 = boonlib('rbmap',51);
    colormap(cmap2);
    shading flat;
    g1=colorbar;
%     set(g1,'tick','in','box','on', 'LineWidth'  , 1,'TickLength'  , [.012 .012]);
%     set(gca,'tick','out','box','on','TickLength'  , [.01 .01], 'LineWidth'   , 1.5);
    contour(xi2,yi2,evalbox,'y-','linewidth',1); hold on;
%     title(hout{2},'interpreter', 'latex','Fontsize',13);
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    cbfreeze(g1)
    freezeColors;
    xlim([-70 70])
    ylim([-70 70])
    
    
    axes(ha(3))    % similar with subplot(221)
    pcolor(xi2,yi2,DIF);
    caxis([-5 8]); 
    cmap3 = hsv(32);
    cmap3(1,:) = [0.900    0.3        0];
    cmap3(2,:) = [0.9000    0.4         0];
    cmap3(3,:) = [0.90000    0.55         0];
    colormap(cmap3);
    hold on;
    shading flat;
    g1=colorbar;
%     set(g1,'tick','in','box','on', 'LineWidth'  , 1,'TickLength'  , [.012 .012]);
%     ylabel('Meridional(km)');
%     set(gca,'tick','out','box','on','TickLength'  , [.01 .01], 'LineWidth'   , 1.5);
    contour(xi2,yi2,evalbox,'y-','linewidth',1); hold on;
%     title(hout{3},'interpreter', 'latex','Fontsize',13);
    cbfreeze(g1)
    set(gca,'XTickLabel',[]);
    freezeColors;
    xlim([-70 70])
    ylim([-70 70])
    

    axes(ha(4))    % similar with subplot(221)
    axis square;
    r_HV=RHO;
    thres = [0.3 0.42 0.55 0.75 0.89 0.95 0.98 0.99 0.995 0.999 1];
    nlev = length(thres);
    clim = [0 nlev+1/12];
    xx=zeros(size(r_HV));
    xx(r_HV<0.3 ) =1;
    xx(r_HV>0.3 & r_HV<=0.42) =0.2;
    xx(r_HV>0.42 & r_HV<=0.55) =0.3;
    xx(r_HV>0.55 & r_HV<=0.75) =0.4;
    xx(r_HV>0.75 & r_HV<=0.89) =0.5;
    xx(r_HV>0.89 & r_HV<=0.95) =0.6;
    xx(r_HV>0.95 & r_HV<=0.98)=0.7;
    xx(r_HV>0.98 & r_HV<=0.99)=0.8;
    xx(r_HV>0.99 & r_HV<=0.995)=0.9;
    xx(r_HV>0.995 & r_HV<=0.999)=1.0;
    xx(r_HV>0.999 & r_HV<=1 )=1.1;
    cmap = boonlib('czmap',17);
    cmap = [0 0 0; 0.4 0.4 0.4; cmap(2:nlev-3,:);  0.92  0.4  0.67];
    cmap = [cmap(1:end-1,:); 0.92  0.4  0.67 ];
    colormap(cmap);  
    xx(isnan(RHO)==1)=nan;
    pcolor(xi2,yi2,xx);
    hold on;
    shading flat;
    g1=colorbar;
%     set(g1,'tick','in','box','on', 'LineWidth'  , 1,'TickLength'  , [.012 .012]);
%     set(gca,'tick','out','box','on','TickLength'  , [.01 .01], 'LineWidth'   , 1.5);
    contour(xi2,yi2,evalbox,'y-','linewidth',1); hold on;
%     title(hout{4},'interpreter', 'latex','Fontsize',13);
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    cbfreeze(g1)
    freezeColors;
    xlim([-70 70])
    ylim([-70 70])
    
    axes(ha(5))    % similar with subplot(221)
    axis square;
    cmap2 = boonlib('rbmap',51);
    colormap(cmap2); 
    pcolor(xi2,yi2,SDV);
    hold on;
    shading flat;
    caxis([0 10]); 
    g1=colorbar;
%     set(g1,'tick','in','box','on', 'LineWidth'  , 1,'TickLength'  , [.012 .012]);
%     xlabel('Zonal(km)');
%     ylabel('Meridional(km)');
%     set(gca,'tick','out','box','on','TickLength'  , [.01 .01], 'LineWidth'   , 1.5);
    contour(xi2,yi2,evalbox,'y-','linewidth',1); hold on;
%     title(hout{5},'interpreter', 'latex','Fontsize',13);
    cbfreeze(g1)
    freezeColors;
    xlim([-70 70])
    ylim([-70 70])
    

    axes(ha(6))    % similar with subplot(221)
    axis square;
    cmap2 = boonlib('rbmap',51);
    colormap(cmap2); 
    pcolor(xi2,yi2,SDP);
    hold on;
    shading flat;
    caxis([0 100]); 
    g1=colorbar;
%     set(g1,'tick','in','box','on', 'LineWidth'  , 1,'TickLength'  , [.012 .012]);
%     xlabel('Zonal (km)');
%     set(gca,'tick','out','box','on','TickLength'  , [.01 .01], 'LineWidth'   , 1.5);
    contour(xi2,yi2,evalbox,'y-','linewidth',1); hold on;
%     title(hout{6},'interpreter', 'latex','Fontsize',13);
    cbfreeze(g1)
    set(gca,'YTickLabel',[]);
    freezeColors;
    xlim([-70 70])
    ylim([-70 70])
    
    
    set(gcf,'color','w');
    set(gcf, 'PaperPositionMode','auto');
%     set(gcf,'render','painter');
%     set(gcf, 'Renderer', 'zbuffer')
    frame = getframe(figure(m));
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,plotout,'png');
%     close(fig)
    end
    end