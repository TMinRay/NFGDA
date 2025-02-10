NF00_header;

fig_dir=['./tracking_points/FIG'];
export_preds_dir = ['./tracking_points/nf_preds'];

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
        matout = [exp_preds_event '\' num2str(m,'%02i') '.mat']; 
        load(matout,"xi2","yi2","REF","nfout","evalbox");

        truescr=evalbox;
        falsescr=(truescr-1).*-1;

        tp=zeros(401,401);
        fp=tp;

       tp=nfout.*truescr.*indcirc;   
       fp=nfout.*falsescr.*indcirc;
       REF=REF.*indcirc;
       evalbox=evalbox.*indcirc;

        
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
        

        plot(xi2(logical(nfout)&logical(indcirc)),yi2(logical(nfout)&logical(indcirc)),'k.','linewidth',2); hold on;
        % plot(xi2(logical(tp)),yi2(logical(tp)),'k.','linewidth',2); hold on;
        % plot(xi2(logical(fp)),yi2(logical(fp)),'r.','linewidth',2); hold on;

    %     plot(xi2(logical(tm)),yi2(logical(tm)),'k.','linewidth',3); hold on;
    %     plot(xi2(logical(fm)),yi2(logical(fm)),'m.','linewidth',2); hold on;
            contour(xi2,yi2,evalbox,'y-','linewidth',1); hold on;
        % plot(xt,yt,'c--','linewidth',2); hold on;   
        xlim([-70 70])
        ylim([-70 70])

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
        % clear a11 truescr falsescr trueregion PARITP;
        close(figure(m))
    
    end
end

