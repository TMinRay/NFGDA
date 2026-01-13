NF00_header;

header = ini2struct('NFGDA.ini');
case_name = header.settings.case_name;
nf_pred_dir = header.settings.export_preds_dir;
fig_dir = header.settings.fig_dir;

[xi2,yi2] = meshgrid(-100:0.5:100,-100:0.5:100);

% indcirc = double(sqrt(xi2.^2+yi2.^2)<=70);
indcirc = double(sqrt(xi2.^2+yi2.^2)<=100);

dtdout = fullfile(fig_dir, case_name);
if not(isfolder(dtdout))
    mkdir(dtdout);
end

data_path = fullfile(nf_pred_dir, case_name);
nf_preds = {dir(fullfile(data_path,'nf_pred*')).name};
nf_preds = sort(nf_preds);

for m=1:length(nf_preds)
    ppi_id = nf_preds{m};
    load(fullfile(data_path,ppi_id),"xi2","yi2","REF","nfout","evalbox");

    % matout = [exp_preds_event '\' num2str(m,'%02i') '.mat']; 
    % load(matout,"xi2","yi2","REF","nfout","evalbox");

    truescr = evalbox;
    falsescr = -(truescr-1);

    tp=zeros(401,401);
    fp=tp;

   tp=nfout.*truescr.*indcirc;
   fp=nfout.*falsescr.*indcirc;
   REF=REF.*indcirc;
   evalbox=evalbox.*indcirc;

    fig=figure(m);
    set(fig,'Position',[100 100 500 480]);
    % ha = tight_subplot(1,1,[.05 .05],[.05 .05],[.05 .05]);
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
    % xlim([-70 70])
    % ylim([-70 70])
    xlim([-100 100])
    ylim([-100 100])

    % xlabel('Zonal (km)','Fontsize',14);
    % ylabel('Meridional (km)','Fontsize',15);
    % set(gca,'TickDir','out','box','on','TickLength'  , [.01 .01], 'LineWidth' , 2);
    set(gca,'TickDir','out','box','on','TickLength'  , [.01 .01], ...
        'LineWidth', 2, 'XTick', -100:20:100,'YTick', -100:20:100);

    ppi_name = ppi_id(12:end);
    tstamp = datetime(ppi_name(5:12), 'InputFormat', 'yyyyMMdd', 'TimeZone', 'UTC') + ...
            timeofday(datetime(ppi_name(14:19), 'InputFormat', 'HHmmss'));
    radar_id = ppi_name(1:4);
    ppi_desc = [radar_id ', ' char(tstamp, 'MM/dd/yyyy, HH:mm:ss z')];

    title([ppi_desc '  NFGF'],'Fontsize',13);
    % set(gcf,'color','w');
    % set(gcf, 'PaperPositionMode','auto');
    % set(gcf,'render','painter')

    plotout=[ dtdout '\' ppi_id(1:end-4) '.png'];
    print(plotout, '-dpng');
    % frame = getframe(figure(m));
    % im = frame2im(frame);
    % [imind,cm] = rgb2ind(im,256);
    % imwrite(imind,cm,plotout,'png');
    % clear a11 truescr falsescr trueregion PARITP;
    close(figure(m))

end
% end

