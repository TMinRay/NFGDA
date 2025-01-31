NF00_header
dtdout=['./NFFIG/stat'] ;                                                                                      

mANALGSTout=[matPATH '/HIST/stat3regions.mat'];
load(mANALGSTout,'gf','st','ca');
mANALGSTout=[matPATH '/HIST/stat2regions.mat'];
load(mANALGSTout,'totg','toto');

% toto={st; ca;};

for i=1:6
    stc=st{i};
    cac=ca{i};
    gfc=gf{i};
    % rmin=min([nanmin(stc(:)), nanmin(cac(:)), nanmin(gfc(:)) ]);
    % rmax=max([nanmax(stc(:)), nanmax(cac(:)), nanmax(gfc(:)) ]);
    rmin=min([min(stc(:)),min(cac(:)),min(gfc(:))]);
    rmax=max([max(stc(:)),max(cac(:)),max(gfc(:))]);
    
    srange{i}=[rmin:(rmax-rmin)/50:rmax];
    srange2{i}=[rmin; rmax];
    
        totoc=toto{i};
    % ssrange{i}=[nanmin(totoc(:)):(nanmax(totoc(:))-nanmin(totoc(:)))/30:nanmax(totoc(:))];
    % ssrange2{i}=[nanmin(totoc(:)); nanmax(totoc(:))*0.9];
    min_tot=min(totoc(:),[],"omitmissing");
    max_tot=max(totoc(:),[],"omitmissing");
    ssrange{i}=[min_tot:(max_tot-min_tot)/30:max_tot];
    ssrange2{i}=[min_tot; max_tot*0.9];
    
end
    srange{4}=[0:0.03:1];
    srange2{2}=[0; 0.8];
    srange{6}=[0:(360-(0))/300:360];
    srange2{6}=[0; 180];
    srange{5}=[0:(15-(-0))/100:15];
    srange2{5}=[-0.1; 5];
    srange2{2}=[0; 0.4];
    srange2{3}=[-4; 8.5];
    ssrange{4}=[0:0.02:1.1];
    ssrange2{4}=[0; 0.99];
        ssrange2{3}=[-4; 8];
    ssrange{7}=[-30:(35-(-30))/30:35];
    ssrange2{7}=[-30; 35];
    ssrange{6}=[0:(360-(0))/200:360];
    ssrange2{6}=[0; 180];
    ssrange{5}=[0:(15-(-0))/50:15];
    ssrange2{5}=[-0.1; 10];
    
    ssrange2{2}=[0; 1];
    ssrange2{3}=[-4; 8.1];


  aout{1}=[dtdout 'sdphiHIST6.eps'];
  
  xa{1}=['{$Z$ (dBZ)}'];
  xa{4}=['{$\rho_{hv}$}'];
  xa{3}=['{$Z_{DR}$ (dB)}'];
  xa{2}=['{ $\beta$}'];
  xa{5}=['{ SD($v_{r}$) (ms$^{-1}$)}'];
  xa{6}=['{ SD($\phi_{DP}$) ($^{\circ}$)}'];
  hout{7}=['{Normalized \medskip histogram \medskip of \medskip $\Delta \, Z$}'];
  
    

 lsize=20; 
 charsize=15;

 fig=figure(1); 
 set(fig,'Position',[100 100 900 850]);
%  ha = tight_subplot(3,2,[.05 .05],[.07 .07],[.07 .1]);
  ha = tight_subplot(3,2,[.07 .05],[.07 .03],[.07 .03]);


for i=1:6
%      nanmean(totg{i})
%     nanstd(gf{i})
%      nanmean(toto{i})
%     nanstd(st{i})
%     nanmean(ca{i})
%     nanstd(ca{i})
[snum_gt,sall_gt]=hist(totg{i},ssrange{i});
[snum_oth,sall_oth]=hist(toto{i},ssrange{i});
slengt=sum(snum_gt);
slenoth=sum(snum_oth);
    
[num_gt,all_gt]=hist(gf{i},srange{i});
[num_st,all_st]=hist(st{i},srange{i});
[num_ca,all_ca]=hist(ca{i},srange{i});

 lengt=sum(num_gt);
 lenst=sum(num_st);
 lenca=sum(num_ca);

 axes(ha(i))
 plot(sall_gt,snum_gt./slengt,'r','linewidth',3);
 hold on;
 plot(sall_oth,snum_oth./slenoth,'b','linewidth',3);
 xlabel(xa{i},'interpreter', 'latex','Fontsize',charsize);
 
%  plot(all_gt,num_gt./lengt,'r');
 hold on;
 plot(all_st,num_st./lenst,'color',[0 .5 0],'linewidth',2);
 plot(all_ca,num_ca./lenca,'k--','linewidth',2);

 xlim([ srange2{i}]);
   if i==1
 xlim([ -5 50]);
   end

     if i==3
 xlim([ -4 8]);
  end
 
 
  if i==4
 xlim([ 0.5 1.0]);
  ylim([ 0 0.3]);
  end
 
   if i==6
 xlim([ -0.1 40]);
   ylim([ 0 0.31]);
  end
 
  
 if i==2
 ylim([ 0 0.5]);
 xlim([ -0 0.5]);
 h_legend=legend('GF','non-GF', 'Storm','Clear-air');
 legend boxoff
 set(h_legend,'interpreter', 'latex','Fontsize',charsize,'location','north');
 
 end    

%  title(hout{i},'interpreter', 'latex','Fontsize',lsize);
 if i==2 || i==4 || i==6
%  set(gca, 'YTicklabel', [ ])
 else
      ylabel('Normalized \medskip Occurence','interpreter', 'latex','Fontsize',charsize);
 
 end 

end
set(gcf,'color','w');
set(gcf, 'PaperPositionMode','auto');
 set(gcf,'render','painter')
%  frame = getframe(figure(1));
% im = frame2im(frame);
% [imind,cm] = rgb2ind(im,256);
% imwrite(imind,cm,aout{1},'png');

 print(gcf,'-depsc',aout{1});
%  close(fig)

