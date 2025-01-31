NF00_header;

dtdout=['../NFFIG/'];


%     for cindex=1:numel(ttable(:,1));    
    for cindex=1:1
        
    PUTDAT=ttable(cindex,:);
%     for m=startm:endm   
    for m=3:3   
    minputNF=[matPATH '/CART/inputNF' PUTDAT num2str(m,'%02i') '.mat']; 
    load(minputNF, 'inputNF'); 
% %  PAR(0~1) 1.REF 2.BETA 3.Z_DR 4.rho_hv 
% %  5.SD(v_r) 6.SD(phi_DP) // GRID (401*401)
        mevalbox=[ matPATH '/EVAL/newevalbox' PUTDAT num2str(m,'%02i') '.mat']; 
        load(mevalbox,'evalbox');   
        mhandpick=[ matPATH '/HANDPICK/handpick' PUTDAT num2str(m,'%02i') '.mat']; 
        load(mhandpick,'handpick');
        
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
    xlim([xgrd1 xgrd2])
    ylim([ygrd1 ygrd2])
    
    
    
    hout{1}=['{ .}'];
    hout{2}=['{ .}'];
    
    REF=inputNF(:,:,1); 
    
   
    fig=figure(m);
    set(fig,'Position',[100 100 630 270]);
    ha = tight_subplot(1,2,[.05 .03],[.09 .06],[.09 .03]);
    plotout=[ dtdout '/fig1.png'];   
    
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
    xlabel('Meridional(km)');
    set(gca,'TickDir','out','box','on','TickLength'  , [.01 .01], 'LineWidth'   , 1.5);
%     contour(xi2,yi2,evalbox,'y-','linewidth',1); hold on;
    title(hout{1},'interpreter', 'latex','Fontsize',13);
    cbfreeze(g1);
    freezeColors;
    xlim([xgrd1 xgrd2])
    ylim([ygrd1 ygrd2])
    

    axes(ha(2))    % similar with subplot(221)
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
    xlabel('Meridional(km)');
    set(gca,'TickDir','out','box','on','TickLength'  , [.01 .01], 'LineWidth'   , 1.5);
    contour(xi2,yi2,evalbox,'y-','linewidth',1); hold on;
    contour(xi2,yi2,handpick,'b-','linewidth',1); hold on;
    title(hout{1},'interpreter', 'latex','Fontsize',13);
    cbfreeze(g1);
    freezeColors;
    xlim([xgrd1 xgrd2])
    ylim([ygrd1 ygrd2])
        
    
    
    
    set(gcf,'color','w');
    set(gcf, 'PaperPositionMode','auto');
%     set(gcf,'render','painter');
%     set(gcf, 'Renderer', 'zbuffer')
    frame = getframe(figure(m));
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,plotout,'png');
    close(fig)
    end
end