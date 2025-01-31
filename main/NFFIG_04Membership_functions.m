run NF00_header;
charsize=15;
% tmpout=[pngPATH2 '/MF/REF.eps'];
tmpout=['./NFFIG/'];
Range1=0:0.01:1  ;
Range2=-20:0.2:70  ;
Range3=0:0.001:1  ;
Range4=-4:0.13:8  ;
Range5=0:0.1:10  ;
Range6=0:1.5:150  ;

% gf1=[0.308714780257392 0.443732960021408];
% gf2=[4.11626734773716 15.0816697718829];
% gf3=[0.322820389562963 0.607429069609535];
% gf4=[2.83468941228872 5.11404048855749];
% gf5=[3.04779895757988 3.42115689385854];
% gf6=[20.8206753838925 32.0004366895764];

gf1=[0.343011474971548 0.343013205344144];
gf2=[4.33829177506464 15.0283112304065];
gf3=[0.270037916522748 0.570877577406138];
gf4=[2.86421357321656 5.13168272670849];
gf5=[2.98574819999131 3.44033892630897];
gf6=[20.8241488634125 32.00808009822];


% ngf1=[0.0557963275716467 -0.0127740952571925];
% ngf2=[16.5007809112174 33.6995101991228];
% ngf3=[0.26409745431802 0.951207741805216];
% ngf4=[2.13655172522198 1.2630295937188];
% ngf5=[2.9717920154568 1.19864850662713];
% ngf6=[20.6900928328566 3.55112530444254];

ngf1=[0.0777363037870736 0.0012529062358556];
ngf2=[16.5008839478689 33.6994186164538];
ngf3=[0.269854461486406 0.946117835249615];
ngf4=[2.13691930702597 1.26307356423026];
ngf5=[2.97222823213022 1.19934594258286];
ngf6=[20.6901423649177 3.55116611243056];

GGF1=gaussmf(Range1,gf1);
GGF2=gaussmf(Range2,gf2);
GGF3=gaussmf(Range3,gf3);
GGF4=gaussmf(Range4,gf4);
GGF5=gaussmf(Range5,gf5);
GGF6=gaussmf(Range6,gf6);

NGF1=gaussmf(Range1,ngf1);
NGF2=gaussmf(Range2,ngf2);
NGF3=gaussmf(Range3,ngf3);
NGF4=gaussmf(Range4,ngf4);
NGF5=gaussmf(Range5,ngf5);
NGF6=gaussmf(Range6,ngf6);


% del=abs(AGGF3-ANGGF3);
% % Range3=0:0.001:1  ;
% % for i=1:numel(Range3)
%     
% [kk]=find(AGGF3-ANGGF3==min(del(:)));
% Range3(kk)

xa{1}=['{$Z$ (dBZ)}'];
  xa{4}=['{$\rho_{hv}$}'];
  xa{3}=['{$Z_{DR}$ (dB)}'];
  xa{2}=['{ $\beta$}'];
  xa{5}=['{ SD($v_{r}$) (ms$^{-1}$)}'];
  xa{6}=['{ SD($\phi_{DP}$) ($^{\circ}$)}'];


VARIL=[1:6];
% VARIL=[2];
fig=figure(1);
set(fig,'Position',[10 10 800 700]);
ha = tight_subplot(3,2,[.07 .05],[.08 .03],[.07 .05]);
for kkk=1:numel(VARIL)

if VARIL(kkk)==1
axes(ha(2))
% text_title2 = [ 'Membership  \medskip functions \medskip of \medskip {$\beta$}'];

plot(Range1,GGF1,'r','Linewidth',3) ; hold on;
plot(Range1,NGF1,'b','Linewidth',2)
 xlabel(xa{2},'interpreter', 'latex','Fontsize',charsize);

% xlbl='{$\beta$}';
elseif VARIL(kkk)==2
axes(ha(1))
% text_title1 = [ 'Membership  \medskip functions \medskip of \medskip {$Z$ (dB$Z$)}'];
plot(Range2,GGF2,'r','Linewidth',3) ; hold on;
plot(Range2,NGF2,'b','Linewidth',2)
xlabel(xa{1},'interpreter', 'latex','Fontsize',charsize);

% xlbl='{$Z$} (dBZ)';
xlim([-20 70])

elseif VARIL(kkk)==3
axes(ha(4))
plot(Range3,GGF3,'r','Linewidth',3) ; hold on;
plot(Range3,NGF3,'b','Linewidth',2)
xlabel(xa{4},'interpreter', 'latex','Fontsize',charsize);

elseif VARIL(kkk)==4
axes(ha(3))
% text_title3 = [ 'Membership  \medskip functions \medskip of \medskip {$Z_{DR}$ (dB)}'];
plot(Range4,GGF4,'r','Linewidth',3) ; hold on;
plot(Range4,NGF4,'b','Linewidth',2)
% xlbl='{$Z_{DR}$ (dB)}';
 xlabel(xa{3},'interpreter', 'latex','Fontsize',charsize);

elseif VARIL(kkk)==5
axes(ha(5))
% text_title5 = [ 'Membership  \medskip functions \medskip of \medskip {SD$(v_{r})$ (ms$^{-1})$}'];
plot(Range5,GGF5,'r','Linewidth',3) ; hold on;
plot(Range5,NGF5,'b','Linewidth',2)
% xlbl='{SD{$(v_{r})$} (ms{$^{-1})$}}';
 xlabel(xa{5},'interpreter', 'latex','Fontsize',charsize);

elseif VARIL(kkk)==6
axes(ha(6))
% text_title6 = [ 'Membership  \medskip functions \medskip of \medskip {SD$(\Phi_{DP})$ ($^\circ$)}'];
plot(Range6,GGF6,'r','Linewidth',3) ; hold on;
plot(Range6,NGF6,'b','Linewidth',2)
 xlabel(xa{6},'interpreter', 'latex','Fontsize',charsize);

% xlbl='{SD$(\phi_{DP})$ ($^\circ$)}';


end

ylim([0 1.1])
% grid on;
if VARIL(kkk)==1
h_legend=legend('$\mu^{G}_{2}$', '$\mu^{N}_{2}$');
set(h_legend,'interpreter', 'latex','Fontsize',15);
% title(text_title2,'interpreter', 'latex','Fontsize',15);
elseif VARIL(kkk)==2
h_legend=legend('$\mu^{G}_{1}$', '$\mu^{N}_{1}$');
set(h_legend,'interpreter', 'latex','Fontsize',15);
% title(text_title1,'interpreter', 'latex','Fontsize',15);
ylabel('{Membership degree}','interpreter', 'latex','Fontsize',14);
% title(text_title1,'interpreter', 'latex','Fontsize',15);
elseif VARIL(kkk)==3
h_legend=legend('$\mu^{G}_{4}$', '$\mu^{N}_{4}$');
set(h_legend,'interpreter', 'latex','Fontsize',15,'Location','best'); 
% title(text_title4,'interpreter', 'latex','Fontsize',15);
elseif VARIL(kkk)==4
h_legend=legend('$\mu^{G}_{3}$', '$\mu^{N}_{3}$');
set(h_legend,'interpreter', 'latex','Fontsize',15,'Location','best');    
ylabel('{Membership degree}','interpreter', 'latex','Fontsize',14);
% title(text_title3,'interpreter', 'latex','Fontsize',15);
elseif VARIL(kkk)==5
h_legend=legend('$\mu^{G}_{5}$', '$\mu^{N}_{5}$');
set(h_legend,'interpreter', 'latex','Fontsize',15);
ylabel('{Membership degree}','interpreter', 'latex','Fontsize',14);
% title(text_title5,'interpreter', 'latex','Fontsize',15);
elseif VARIL(kkk)==6
h_legend=legend('$\mu^{G}_{6}$', '$\mu^{N}_{6}$');
set(h_legend,'interpreter', 'latex','Fontsize',15);
% title(text_title6,'interpreter', 'latex','Fontsize',15);
end
% end
% grid on;
% h_legend=legend('GF', 'None-GF');
% 
% if VARIL(kkk)==3 | VARIL(kkk)==4
% set(h_legend,'interpreter', 'latex','Fontsize',10,'Location','NorthWest');    
% else
% set(h_legend,'interpreter', 'latex','Fontsize',10);
% end



end
mfpout=[tmpout 'MFS.eps'];
set(gcf,'color','w');
set(gcf, 'PaperPositionMode','auto');
set(gcf,'render','painter')

print(gcf,'-depsc',mfpout,'-loose');
% close(fig)

 