NF00_header;

fig_dir=['./ieee2016/out/'];
export_preds_dir = ['./tracking_points/nf_preds'];

% indcirc=zeros(401,401);
% for i=1:401
%     for j=1:401
%         if sqrt(xi2(i,j)^2+yi2(i,j)^2)<=70
%             indcirc(i,j)=1;
%         else
%             indcirc(i,j)=0;
%         end
%     end
% end
indcirc = double(sqrt(xi2.^2+yi2.^2)<=70);

n=1;

for cindex=1:numel(ttable(:,1));
%     for cindex=7:numel(ttable(:,1));
        
        
    PUTDAT=ttable(cindex,:);

    dtdout = fullfile(fig_dir, PUTDAT);    
    if not(isfolder(dtdout))
        mkdir(dtdout);
    end

    exp_preds_event = fullfile(export_preds_dir,PUTDAT);
    if not(isfolder(exp_preds_event))
        mkdir(exp_preds_event);
    end

    startm=startt(cindex)+1;
    endm=endt(cindex);    
%     for m=startm:startm
    for m=startm:endm   
        
    minputNF=[matPATH '/CART/inputNF' PUTDAT num2str(m,'%02i') '.mat']; 
    load(minputNF, 'inputNF'); 
% %  PAR(0~1) 1.REF 2.BETA 3.Z_DR 4.rho_hv 
% %  5.SD(v_r) 6.SD(phi_DP) // GRID (401*401)

    REF=inputNF(:,:,1);
    
%     if GFYN(cindex)==1
%     mevalbox=[ matPATH '/EVAL/evalbox' PUTDAT num2str(m,'%02i') '.mat']; 
    mevalbox=[ matPATH '/EVAL/newevalbox' PUTDAT num2str(m,'%02i') '.mat']; 
    if exist(mevalbox)>0
        load(mevalbox,'evalbox');     
        else
            evalbox=zeros(401,401);
    end
    
% % % % % % % % % % % % % % % % % %     else
% % % % % % % % % % % % % % % % % %     evalbox=zeros(401,401);
% % % % % % % % % % % % % % % % % %     end
        mindxy=[matPATH '/OUT/final' PUTDAT num2str(m,'%02i') '.mat'];
         if exist(mindxy)>0
        load(mindxy,'skel_nfout','skel_nfout2');   
        else
            skel_nfout=zeros(401,401);
        end
          
%     mindxyall=[matPATH '/OUT/ysfinal' PUTDAT num2str(m,'%02i') '.mat'];
%     load(mindxyall,'farray');   
    nfout=skel_nfout2;
    clear skel_nfout;

    matout = [exp_preds_event '\' num2str(m,'%02i') '.mat']; 
    save(matout,"xi2","yi2","REF","nfout","evalbox");
    
%         mindxy2=['../mat/OUT/MIGFAfinal' PUTDAT num2str(m,'%02i') '.mat'];
%         if exist(mindxy2)>0
%         load(mindxy2,'skel_nfout');   
%         else
%             skel_nfout=zeros(401,401);
%         end
%     mfout=skel_nfout;
%     clear farray;


    truescr=evalbox;
    falsescr=(truescr-1).*-1;

    tp=zeros(401,401);
    fp=tp;
%     tm=tp;
%     fm=tp;
    
% 
%     for ti=1:401
%         for tj=1:401
%            tp(ti,tj)=nfout(ti,tj)*truescr(ti,tj)*indcirc(ti,tj);   
%            fp(ti,tj)=nfout(ti,tj)*falsescr(ti,tj)*indcirc(ti,tj);
% %            tm(ti,tj)=mfout(ti,tj)*truescr(ti,tj)*indcirc(ti,tj);   
% %            fm(ti,tj)=mfout(ti,tj)*falsescr(ti,tj)*indcirc(ti,tj);
%            REF(ti,tj)=REF(ti,tj)*indcirc(ti,tj);
%            evalbox(ti,tj)=evalbox(ti,tj)*indcirc(ti,tj);
%         end
%     end
       tp=nfout.*truescr.*indcirc;   
       fp=nfout.*falsescr.*indcirc;
       REF=REF.*indcirc;
       evalbox=evalbox.*indcirc;

  
% %     fig=figure(m);
% %     set(fig,'Position',[100 100 500 480]);
% %     ha = tight_subplot(1,1,[.05 .05],[.05 .05],[.05 .05]);
% %     REF(indcirc<1)=nan;
% %     pcolor(xi2,yi2,REF); 
% %     shading flat;
% %     hold on;
% %     clim = [-5 65];
% %     cmap = boonlib('zmap3',35);
% %     colormap((cmap));
% %     
% % 
% % 
% %     plot(xi2(logical(tp)),yi2(logical(tp)),'k.','linewidth',2); hold on;
% %     plot(xi2(logical(fp)),yi2(logical(fp)),'r.','linewidth',2); hold on;
% %     contour(xi2,yi2,evalbox,'y-','linewidth',1); hold on;
% %     xlim([-70 70])
% %     ylim([-70 70])
% %     
% %     % plot(xi2(logical(tm)),yi2(logical(tm)),'c.','linewidth',3); hold on;
% %     % plot(xi2(logical(fm)),yi2(logical(fm)),'m.','linewidth',2); hold on;
% %     % plot(xt,yt,'c--','linewidth',2); hold on;   
% % 
% %     
% % %     PNGgridlim4;
% %     
% %     
% %     xlabel('Zonal (km)','Fontsize',14);
% %     ylabel('Meridional (km)','Fontsize',15);
% %     set(gca,'tick','out','box','on','TickLength'  , [.01 .01], 'LineWidth' , 2);
% %     
% %     text_t=['{Case' num2str(cindex) ', t=' num2str(m-(startm-1)) ' \medskip, NFout}'];
% %     title(text_t,'interpreter', 'latex','Fontsize',15);
% %     set(gcf,'color','w');
% %     set(gcf, 'PaperPositionMode','auto');
% %     set(gcf,'render','painter')
% %     plotout=[ '../ieee2016/out/nf_' num2str(cindex) '_t_' num2str(m,'%02i') '.png'];       
% %     
% %     frame = getframe(figure(m));
% %     im = frame2im(frame);
% %     [imind,cm] = rgb2ind(im,256);
% %     imwrite(imind,cm,plotout,'png');
% %     close(figure(m))
% %     
    
    fig=figure(m);
    set(fig,'Position',[100 100 500 480]);
    ha = tight_subplot(1,1,[.05 .05],[.05 .05],[.05 .05]);
    REF(indcirc<1)=nan;
    pcolor(xi2,yi2,REF); 
    shading flat;
    hold on;
    clim = [-5 65];
    cmap = boonlib('zmap3',35);
    colormap((cmap));
    


    plot(xi2(logical(tp)),yi2(logical(tp)),'k.','linewidth',2); hold on;
    plot(xi2(logical(fp)),yi2(logical(fp)),'r.','linewidth',2); hold on;

%     plot(xi2(logical(tm)),yi2(logical(tm)),'k.','linewidth',3); hold on;
%     plot(xi2(logical(fm)),yi2(logical(fm)),'m.','linewidth',2); hold on;
        contour(xi2,yi2,evalbox,'y-','linewidth',1); hold on;
    % plot(xt,yt,'c--','linewidth',2); hold on;   
    xlim([-70 70])
    ylim([-70 70])
    

    
%     PNGgridlim4;
    
    
    xlabel('Zonal (km)','Fontsize',14);
    ylabel('Meridional (km)','Fontsize',15);
    set(gca,'TickDir','out','box','on','TickLength'  , [.01 .01], 'LineWidth' , 2);
    
    text_t2=['{Case' num2str(cindex) ', t=' num2str(m-(startm-1)) ' \medskip, MFout}'];
    title(text_t2,'interpreter', 'latex','Fontsize',15);
    set(gcf,'color','w');
    set(gcf, 'PaperPositionMode','auto');
    set(gcf,'render','painter')
    plotout=[ dtdout '\' 'mf_' num2str(cindex) '_t_' num2str(m,'%02i') '.png'];       
    
    frame = getframe(figure(m));
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,plotout,'png');
    clear a11 truescr falsescr trueregion PARITP;
    close(figure(m))
    
    end 
end

