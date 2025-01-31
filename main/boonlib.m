% BOONLIB - Boon Leng's library of functions
%
% 5/31/2010 - Code clean up for MATLAB 2010a
%           - Added blackngold().
%           - Added ezmap().
%           - Added ezmap2().
%           - Updated cspace() to have shade counts
%
% 8/16/2009 - Code clean up.
%           - Replaced some depreciating function calls.
%           - Added zmapn().
%           - Added zmapnx().
%           - Added rgmapn().
%
% 3/7/2009  - Added czmapx().
%
% 4/18/2008 - Added rcmap().
%           - Added bjetmapx().
%
% 4/13/2008 - Updated czmap().
%           - Added rgmap2().
%           - Added spidergrid.
%
% ------------------------------ Figures -------------------------------
% BMOVEFIG     - Move figure/axes without changing the size
%                boonlib('bmovefig',HANDLES,[X_OFFSET Y_OFFSET])
%
% BMOVEFIGANIM - Move figure/axes without changing the size
%                boonlib('bmovefiganim',HANDLES,[X_OFFSET Y_OFFSET])
%
% PLACEFIG     - Place figure/axes in an absolute position
%                boonlib('placefig',HANDLES,[X Y])
%
% PLACEFIGUL   - Place figure/axes in an absolute position wrt upper-left
%                boonlib('placefigul',HANDLES,[X Y])
%
% BZOOM        - Scale up/down a figure.
%                boonlib('bzoom',FIGURE_HANDLE,ZOOM_FACTOR)
%
% BSIZEWINDOW  - Resize the figure/axes and keep upper left position
%                boonlib('bsizewindow',FIGURE_HANDLE,[WIDTH HEIGHT])
%
% ----------------------------- Colormaps ------------------------------
%
% RBMAP        - Red to White to Blue Map.  Good for symmetric +/- range
%                (e.g. radial velocity)
%                cmap = boonlib('rbmap',NUMBER_OF_COLORS)
%
% GOMAP        - Green to Yellow to Orange, symmetric.
%                cmap = boonlib('gomap',NUMBER_OF_COLORS)
%
% OGMAP        - Orange to White to Green, symmetric.
%                cmap = boonlib('ogmap',NUMBER_OF_COLORS)
%
% RGMAP        - Red to Gray to Green Map.  Common for radial velocity
%                cmap = boonlib('rgmap',NUMBER_OF_COLORS)
%
% RGMAP2       - Red to Gray to Green Map 2.  Similar to one in NCDC Viewer.
%                Best for odd number of bins.
%                cmap = boonlib('rgmap2',NUMBER_OF_COLORS)
%
% RGMAPN       - Red to Gray to Green Map similar to NCDC Viewer.
%                Best for ven number of bins.
%                cmap = boonlib('rgmapn',NUMBER_OF_COLORS)
%
% RGMAPINV     - Using the above map, flipped around the center, Gray stays.
%                Good for highlighting middle shade
%                cmap = boonlib('rgmapinv',NUMBER_OF_COLORS)
%
% RAPMAP       - Colormap that prints well on both color and grayscale.
%                Black --> Blue --> Magenta --> Red --> Yellow --> White
%                Generated using ideas from Rappaport's paper.
%                cmap = boonlib('rapmap',NUMBER_OF_COLORS)
%
% ASYMAP       - Asymmetric Blue to White to Red Map.
%                (e.g. radial velocity with more emphasis on negative)
%                cmap = boonlib('asymap',[MIN MAX])
%                cmap = boonlib('asymap',[MIN MAX MINCLIP MAXCLIP])
%
% ZMAP         - A fairly standard reflectivity 16-shade colormap for radar images
%                cmap = boonlib('zmap') creates a colormap for [0 75] dBZ; or
%                cmap = boonlib('zmap',[MIN MAX]) for [MIN MAX] dBZ
%
% ZMAP2        - Reflectivity colormap used by Plymouth State Weather Center
%                cmap = boonlib('zmap2') creates a colormap for [0 75] dBZ; or
%                cmap = boonlib('zmap2',[MIN MAX]) for [MIN MAX] dBZ
%
% ZMAP3        - Reflectivity colormap used by WDSS-II
%                cmap = boonlib('zmap3') creates a colormap for [0 75] dBZ; or
%                cmap = boonlib('zmap3',[MIN MAX]) for [MIN MAX] dBZ
%
% ZMAPN        - Reflectivity colormap used by NOAA Weather and Climate Toolkit
%                cmap = boonlib('zmapn') creates a colormap for [-2.5 77.5] dBZ; or
%
% CARBMAP      - Carbone map. Suitable for refractivity.
%                cmap = boonlib('carbmap',NUMBER_OF_COLORS)
%
% CZMAP        - Similar to zmap with flexible number of colors.
%                cmap = boonlib('czmap',NUMBER_OF_COLORS)
%                Requirement: Minimum NUMBER_OF_COLORS = 16
%
% T7MAP        - Dark blue to Green at 0.7, Magenta at 0.7+ to yellow then white
%                cmap = boonlib('t7map',NUMBER_OF_COLORS)
%
% T6MAP        - Dark blue to Green at 0.6, Magenta at 0.6+ to yellow then white
%                cmap = boonlib('t6map',NUMBER_OF_COLORS)
%
% RCMAP        - Red, Orange, Yellow, White, Green, Dark Cyan, Cyan.
%                cmap = boonlib('rcmap',NUMBER_OF_COLORS)
%
% BJETMAP      - A cross between RGMAP and JET
%                cmap = boonlib('bjetmap',NUMBER_OF_COLORS)
%
% BJETMAPX     - Extended BJETMAP
%                cmap = boonlib('bjetmapx',NUMBER_OF_COLORS)
%
% MAPINHEX     - Conver colormap from 0..1 representations to 0 to 255 int8
%                in HEX representation
%                boonlib('mapinhex',COLOR_MAP)
%
% ----------------------------- Drawings -------------------------------
%
% DRAW3DBOX    - Draw a 3D Cube at (x y z) = (0 0 0) and with a size of [1 1 1].
%                h = boonlib('draw3dbox',[X_MIN X_MAX Y_MIN Y_MAX Z_MIN Z_MAX])
%
% DRAWARROW    - Draw a solid arrow using fill() function.
%                Position parameters are in normalized unit of the gcf figure.
%                h = boonlib('drawarrow',[X Y LENGTH ANGLE_RAD])
%                h = boonlib('drawarrow',[X Y LENGTH ANGLE_RAD THICKNESS])
%
% CIRCLE       - Draw circles
%                h = boonlib('circle',[X Y R],NPTS) draws concentric circles at
%                [X Y] with radii of R.  Each circle is drawn with NPTS points.
%                Note: [X Y Z] must be in the size of Nx3 for drawing N circles.
%                boonlib('circle',[X Y Z]) draws circles with default NPTS=90.
%
% CSPACE       - Demonstrate RGB and YUV components of a colormap
%                boonlib('cspace', COLORMAP)
%
% GRAYME       - Set colormap to a Gray version (Preview of Monochrome printout)
%                boonlib('grayme')
%                boonlib('grayme',FIGURE_HANDLE)
%
% SPIDERGRID   - Draws spider-like grids at specified ring spacing and
%                number of lines
%                boonlib('spidergrid')
%                boonlib('spidergrid',AXES,RING_SPACING,NUMBER_OF_LINES)
%
% ------------------------------- Misc ---------------------------------
%
% CHOOSEFILE   - Interface to choose file using dir()
%                boonlib('choosefile',FILEPATH,FILE_CRITERIA)
%                boonlib('choosefile',FILEPATH_FUNCTION,FILE_CRITERIA)
%
%

function [varargout] = boonlib(fn,varargin)
if ~exist('fn','var')
	eval('help boonlib')
    if isempty(whos), feval('demo'); return; end
else
    if nargout==0; feval(fn,varargin{:});
    elseif nargout==1, x = feval(fn,varargin{:}); varargout(1) = {x};
    end
end
return

% ------------------------------- Demo ---------------------------------
function [x] = demo()
feval('mydefault')
h = feval('assignfig','demo');
figure(h)
cm = {'bjetmap', 'bjetmapinv', 'bjetmapxinv',...
    'carbonemap', 'zmap','zmap2', 'zmap3', 'czmap', 'czmapx', 'zmapn', 'zmapnx',...
    'gomap', 'ogmap', 'refmap',...
    'rgmapinv', 'rgmap2', 'rgmapw', 'rgmapn',...
    't7map', 't6map', 't5map',...
    'bgraymap', 'wmap', 'tempmapinv',...
    'brmap', 'rcmap',...
    'asymap0', 'asymap',...
    'vrmap','blackngold','ezmap','ezmap2',...
	'carbmap2',...
    };
idx = 1;
while (idx<length(cm))
    fprintf('Map: %s\n',cm{idx});
    cmap = feval(cm{idx});
    cspace(cmap);
    if (idx==1)
        bzoom(gcf,1.0); drawnow;
        bxwin();
    end
    drawnow;
    idx = idx+1;
end
clf
set(gcf,'Visible','Off')
bsizewin(gcf,[200 200])
placefig(gcf,[200 200])
placefigul(gcf,[200 10]);
cm = {'mapinhexq', 'mapinhex',...
    'draw3dbox','drawarrow','circle','spidergrid',...
    'grayme','greenme','sepiame'};
idx = 1;
while (idx<length(cm))
    feval(cm{idx});
    idx = idx+1;
end
close(gcf);
feval('aintersect',1);
feval('choosefile','.','randomFilenameThatDoesNotExist');
x = feval('isoverlap',[0 0 200 200],[100 100 200 200]);
feval('mydefault2')
return;
    
% ------------------------------ Figures -------------------------------

function bmovefig(arg1,arg2)
if nargin<2, error('Must have two input.'); end
if size(arg2,1) ~= length(arg1), error('Dimension of arg1 and arg2 does not match.'); end
nfig = length(arg1);
tmp = zeros([nfig 4]); for iarg = 1:length(arg1), tmp(iarg,1:4) = get(arg1(iarg),'Position'); end
for iarg = 1:length(arg1)
    set(arg1(iarg),'Position',[tmp(iarg,1:2)+arg2(iarg,1:2) tmp(iarg,3:4)]);
end


function bmovefiganim(arg1,arg2)
if nargin<2, error('Must have two input.'); end
if size(arg2,1) ~= length(arg1), error('Dimension of arg1 and arg2 does not match.'); end
del_step = 5.95943165613589.\arg2;
nfig = length(arg1);
tmp = zeros([nfig 4]); for iarg = 1:length(arg1), tmp(iarg,1:4) = get(arg1(iarg),'Position'); end
for idx = cumsum(1-sin(0:0.5*pi/15:0.49*pi))
    for iarg = 1:length(arg1)
        set(arg1(iarg),'Position',[round(tmp(iarg,1:2)+idx*del_step(iarg,:)) tmp(iarg,3:4)]);
        pause(0.01);
    end
end
return


function placefig(arg1,arg2)
if nargin<2, error('Must have two input.'); end
if size(arg2,1) ~= length(arg1), error('Dimension of arg1 and arg2 does not match.'); end
nfig = length(arg1);
tmp = zeros([nfig 4]); for iarg = 1:length(arg1), tmp(iarg,1:4) = get(arg1(iarg),'Position'); end
for iarg = 1:length(arg1)
    %set(arg1(iarg),'Position',[tmp(iarg,1:2)+arg2(iarg,1:2) tmp(iarg,3:4)]);
	set(arg1(iarg),'Position',[arg2(iarg,1:2) tmp(iarg,3:4)]);
end
return


function placefigul(arg1,arg2)
if nargin<2, error('Must have two input.'); end
if size(arg2,1) ~= length(arg1), error('Dimension of arg1 and arg2 does not match.'); end
% Convert arg2 into UR coordinate system
tmp = get(0,'ScreenSize');
pos = get(arg1,'Position');
if iscell(pos), pos = cat(1,pos{:}); end
newpos = arg2;
newpos(:,2) = tmp(4)-20-arg2(:,2)-pos(:,4);
placefig(arg1,newpos)
return


function bsizewin(h,newsize)
if nargin < 1; fprintf('Please put something.\n'); end
win_position = feval('get',h,'Position');
new_lower_right = win_position(1:2) + [0,win_position(4)-newsize(2)];
feval('set',h,'Position',[new_lower_right,newsize]);
return


function bzoom(h,scale)
if nargin==0; h = feval('gcf'); scale = 2; end
win_position = feval('get',h,'Position');
newsize = floor(win_position(3:4)*scale);
new_lower_right = win_position(1:2) + [0,win_position(4)-newsize(2)];
feval('set',h,'Position',[new_lower_right,newsize]);
return


function bxwin(h,cn,reorder)
if (nargin<1)||isempty(h); h = get(0,'Children'); end
scnsize = get(0,'ScreenSize');
% If there are two screen (==2690 pixels), use my office setting, otherwise single screen settings
if scnsize(3)>2560, x_off = 1280; if (nargin<2)||isempty(cn), cn = [50 80]; end
else x_off = 0;
	if (nargin<2)||isempty(cn)
        cn = [10 40];
	end
end
%scnspace = repmat(0,[scnsize(4) scnsize(3)]);
% Y offset for the dock bar
y_off = 50;
% Buffer space among windows
bf = [10 10];

% Only get figures that still exist
h = h(ishandle(h)&(h~=0));
nfig = numel(h);
if ~exist('reorder','var'), reorder = 1; end
if (reorder)
	% Sort the figure so that they are from big to small
    hsiz = zeros(1,nfig);
	for ifig = 1:nfig
		tmp = get(h(ifig),'Position');
		hsiz(ifig) = tmp(3)*tmp(4);
	end
	[tmp I] = sort(hsiz,2,'descend'); %#ok<ASGLU>
	h = h(I);
end
% 1-st available space is always at the upper-right corner
avspace = scnsize(3:4)-cn;
% Define it first
target = avspace;
f = zeros([nfig 4]);
for ifig = 1:nfig
    f(ifig,1:4) = get(h(ifig),'Position');
    % Number of available space
	nspace = size(avspace,1);
    % Fudge factor to measure remaining area
    ff = zeros([nspace 1]);
    % If there is menubar, penalize more
    if ~strcmpi(get(h(ifig),'Menubar'),'figure'), ybar = 21;
    else ybar = 80;
    end
    % Go through the available space to see if current window fits without overlaps
    for ispace = 1:nspace
    	ff(ispace) = (avspace(ispace,1)-f(ifig,3)-x_off)*(avspace(ispace,2)-f(ifig,4)-y_off-ybar);
        tmp = (avspace(ispace,1)-f(ifig,3)     >=(target(1:ifig-1,1)-f(1:ifig-1,3)) & ...
               avspace(ispace,1)-f(ifig,3)     <=(target(1:ifig-1,1)              ) & ...
               avspace(ispace,2)-f(ifig,4)-ybar>=(target(1:ifig-1,2)-f(1:ifig-1,4)) & ...
               avspace(ispace,2)-f(ifig,4)-ybar<=(target(1:ifig-1,2)              )) | ...
              (avspace(ispace,1)-f(ifig,3)     >=(target(1:ifig-1,1)-f(1:ifig-1,3)) & ...
               avspace(ispace,1)-f(ifig,3)     <=(target(1:ifig-1,1)              ) & ...
               avspace(ispace,2)               >=(target(1:ifig-1,2)-f(1:ifig-1,4)) & ...
               avspace(ispace,2)               <=(target(1:ifig-1,2)              )) | ...
              (avspace(ispace,1)               >=(target(1:ifig-1,1)-f(1:ifig-1,3)) & ...
               avspace(ispace,1)               <=(target(1:ifig-1,1)              ) & ...
               avspace(ispace,2)-f(ifig,4)-ybar>=(target(1:ifig-1,2)-f(1:ifig-1,4)) & ...
               avspace(ispace,2)-f(ifig,4)-ybar<=(target(1:ifig-1,2)              )) | ...
              (avspace(ispace,1)               >=(target(1:ifig-1,1)-f(1:ifig-1,3)) & ...
               avspace(ispace,1)               <=(target(1:ifig-1,1)              ) & ...
               avspace(ispace,2)               >=(target(1:ifig-1,2)-f(1:ifig-1,4)) & ...
               avspace(ispace,2)               <=(target(1:ifig-1,2)              )) | ...
              (target(1:ifig-1,1)-f(1:ifig-1,3)>=(avspace(ispace,1)-f(ifig,3)     ) & ...
               target(1:ifig-1,1)-f(1:ifig-1,3)<=(avspace(ispace,1)               ) & ...
               target(1:ifig-1,2)-f(1:ifig-1,4)>=(avspace(ispace,2)-f(ifig,4)-ybar) & ...
               target(1:ifig-1,2)-f(1:ifig-1,4)<=(avspace(ispace,2)               )) | ...
              (target(1:ifig-1,1)-f(1:ifig-1,3)>=(avspace(ispace,1)-f(ifig,3)     ) & ...
               target(1:ifig-1,1)-f(1:ifig-1,3)<=(avspace(ispace,1)               ) & ...
               target(1:ifig-1,2)              >=(avspace(ispace,2)-f(ifig,4)-ybar) & ...
               target(1:ifig-1,2)              <=(avspace(ispace,2)               )) | ...
              (target(1:ifig-1,1)              >=(avspace(ispace,1)-f(ifig,3)     ) & ...
               target(1:ifig-1,1)              <=(avspace(ispace,1)               ) & ...
               target(1:ifig-1,2)-f(1:ifig-1,4)>=(avspace(ispace,2)-f(ifig,4)-ybar) & ...
               target(1:ifig-1,2)-f(1:ifig-1,4)<=(avspace(ispace,2)               )) | ...
              (target(1:ifig-1,1)              >=(avspace(ispace,1)-f(ifig,3)     ) & ...
               target(1:ifig-1,1)              <=(avspace(ispace,1)               ) & ...
               target(1:ifig-1,2)              >=(avspace(ispace,2)-f(ifig,4)-ybar) & ...
               target(1:ifig-1,2)              <=(avspace(ispace,2)               ));
        %if ~isempty(find(tmp,1)), ff(ispace) = -999999; end
        if any(tmp), ff(ispace) = -99999; end
    end
    % Check the position scores (debug only)
    % disp([avspace repmat([f(ifig,3) f(ifig,4)],[size(avspace,1) 1]) ff])
    % Find out which space can fit it better
    [tmp ispace] = max(ff); %#ok<ASGLU>
    target(ifig,1:2) = avspace(ispace,:)-[0 ybar];
    %scnspace(max(target(ifig,2)-f(ifig,4)+1,1):target(ifig,2),...
    %         max(target(ifig,1)-f(ifig,3)+1,1):target(ifig,1)) = ifig;
    if ispace==nspace
        avspace = [avspace(1:ispace-1,:); ...
                   avspace(ispace,:)-[f(ifig,3)+bf(1) 0]; ...
                   avspace(ispace,:)-[0 f(ifig,4)+bf(2)+ybar]];
    else
        avspace = [avspace(1:ispace-1,:); ...
                   avspace(ispace,:)-[f(ifig,3)+bf(1) 0]; ...
                   avspace(ispace,:)-[0 f(ifig,4)+bf(2)+ybar]; ...
                   avspace(ispace+1:nspace,:)];
    end
end
offset = target-f(:,1:2)-f(:,3:4);
% figure(1)
% imagesc(scnspace);
% set(gca,'YDir','Normal')
if any(offset(:)~=0), 
    if verLessThan('matlab','6.0.0') || isempty(javachk('jvm'))
        feval('bmovefig',h,offset);
    else
        feval('bmovefiganim',h,offset);
    end
end
return

% ----------------------------- Colormaps ------------------------------

function [x] = mybgcolor()
tmp = get(0,'defaultTextColor');
tmp = [0.3 0.59 0.11]*tmp(:);
if tmp>0.5
	x = [0.27 0.23 0.20];
else
	% x = [0.925 0.900 0.850];
	x = [0.89 0.87 0.83];
end
return


function [x] = bjetmap(num)
if nargin < 1,num = size(colormap,1); end
pt = [0.00 0.00 0.00 0.50;...
      0.13 0.00 0.10 1.00;...
      0.37 0.00 1.00 1.00;...
      0.50 1.00 1.00 1.00;...
      0.57 1.00 1.00 0.00;...
      0.88 1.00 0.00 0.00;...
      1.00 0.50 0.00 0.00];
x = feval('fleximap',num,pt);
return


function [x] = bjetmapinv(num)
if nargin < 1,num = size(colormap,1); end
x = bjetmap(num);
x = flipud(x);
return


function [x] = bjetmapx(num)
if nargin < 1,num = size(colormap,1); end
pt = [0.00 0.00 0.80 0.00;...
      0.15 0.00 0.30 0.00;...
      0.20 0.00 0.00 0.50;...
      0.27 0.00 0.10 1.00;...
      0.44 0.00 1.00 1.00;...
      0.50 1.00 1.00 1.00;...
      0.54 1.00 1.00 0.00;...
      0.73 1.00 0.00 0.00;...
      0.80 0.50 0.00 0.00;...
      0.85 0.50 0.00 0.50;...
      1.00 0.90 0.00 0.85];
x = feval('fleximap',num,pt);
return


function [x] = bjetmapxinv(num)
if nargin < 1,num = size(colormap,1); end
x = bjetmapx(num);
x = flipud(x);
return 


function [x] = carbonemap(num)
if ~exist('num','var'), num = 16; end
x = carbmap(num);
return


function [x] = carbmap(num)
if nargin<1,num = size(colormap,1); end
if num<16, fprintf('Cannot generate map with %d elements\n',num); x = []; return; end
num = num-1;
pt = [0.000                  mybgcolor(); ...
      1/num                  0.10 0.00 0.30; ...  % dark blue
      0.125                  0.33 0.06 0.73; ...  % blue-purple
      round(0.25*num-1)/num  0.50 0.45 1.00; ...  % light blue
      round(0.25*num)/num    0.00 0.40 0.00; ...  % green
      0.325                  0.00 0.65 0.00; ...  % mid-green
      round(0.475*num-1)/num 0.70 0.90 0.70; ...  % light green
      round(0.475*num)/num   0.90 0.90 0.90; ...  % gray
      0.575                  1.00 1.00 0.20; ...  % light yellow
      0.675                  0.95 0.60 0.10; ...  % yellowish orage
      round(0.825*num)/num   0.50 0.30 0.25; ...  % brown
      round(0.825*num+1)/num 1.00 0.30 0.51; ...  % 
      1.000                  0.60 0.00 0.05];
x = feval('fleximap',num+1,pt);
return


function [x] = carbmap2(num)
if nargin<1,num = size(colormap,1); end
if num<16, fprintf('Cannot generate map with %d elements\n',num); x = []; return; end
num = num-1;
pt = [0.000                  mybgcolor(); ...
      1/num                  0.3059    0.0980    0.6196; ...  % dark blue
      0.125                  0.4863    0.3608    0.8902; ...  % blue-purple
      round(0.25*num-1)/num  0.4627    0.6902    0.8157; ...  % light blue
      round(0.25*num)/num    0.0588    0.5451    0.0510; ...  % green
      round(0.475*num-1)/num 0.6431    0.8196    0.6392; ...  % light green
      %round(0.475*num)/num   0.9098    0.9059    0.9098; ...  % gray
      round(0.475*num)/num   0.9686    0.8627    0.1137; ...  % gray
      %0.575                  0.9686    0.8627    0.1137; ...  % light yellow
      0.675                  0.8706    0.6471    0.2196; ...  % yellowish orage
      round(0.825*num)/num   0.6275    0.4431    0.3451; ...  % brown
      round(0.825*num+1)/num 0.8588    0.3765    0.4745; ...  % 
      1.000                  0.8000    0.1490    0.2510];
x = feval('fleximap',num+1,pt);
return


function [x] = czmap(num)
if nargin<1,num = 16; end
if num<12, fprintf('Cannot generate map with %d elements\n',num); num = 12; end
num = num-1;
% pt = [0.000                  mybgcolor(); ...
%       1/num                  0.10 1.00 1.00; ...  % cyan blue
%       0.118                  0.00 0.40 1.00; ...  % ocean blue
%       round(0.261*num-1)/num 0.00 0.00 0.75; ...  % dark blue
%       round(0.261*num)/num   0.20 1.00 0.40; ...  % light green
%       0.332                  0.00 0.80 0.00; ...  % green
%       round(0.451*num-1)/num 0.00 0.50 0.00; ...  % dark green
%       round(0.451*num)/num   1.00 1.00 0.20; ...  % light yellow
%       0.5462                 1.00 0.60 0.00; ...  % orange
%       0.6418                 1.00 0.40 0.40; ...  % pink
%       round(0.737*num)/num   0.95 0.00 0.00; ...  % torch read
%       round(0.785*num)/num   0.67 0.00 0.00; ...  % dark red
%       round(0.785*num+1)/num 1.00 0.00 1.00; ...  % magenta
%       (num-1)/num            0.40 0.10 0.60; ...  % dark purple
%       1.000                  1.00 1.00 1.00];
c = (num-1)/num;
d = c/(5-7/num);
m = num-1;
pt = [0.000                          mybgcolor(); ...
      1/num                          0.10 1.00 1.00; ...  % cyan blue
      1/num+(round(d*m)-1)/num       0.00 0.00 0.75; ...  % dark blue
      1/num+(round(d*m))/num         0.00 1.00 0.00; ...  % light green
      1/num+(round(2*d*m)-1)/num     0.00 0.50 0.00; ...  % dark green
      1/num+(round(2*d*m))/num       1.00 1.00 0.00; ...  % yellow
      1/num+(round(3*d*m)-1)/num     1.00 0.50 0.00; ...  % orange
      1/num+(round(3*d*m))/num       1.00 0.00 0.00; ...  % torch red
      1/num+(round(4*d*m)-1)/num     0.50 0.00 0.00; ...  % dark red
      1/num+(round(4*d*m))/num       1.00 0.00 1.00; ...  % magenta
      (num-1)/num                    0.60 0.20 0.80; ...  % purple
      1.000                          1.00 1.00 1.00];
x = feval('fleximap',num+1,pt);
return


function [x] = czmapx(num)
if nargin<1,num = 23; end
if num<16, fprintf('Cannot generate map with %d elements\n',num); num = 16; end
num = num-1;
c = (num-1)/num;
d = c/(7-4/num);
m = num-1;
pt = [0.000                          mybgcolor(); ...
      1/num                          0.80 1.00 1.00; ...
      2/num                          0.80 0.60 0.80; ...  % light purple
      2/num+(round(d*m)-1)/num       0.40 0.20 0.40; ...  % dark purple
      2/num+(round(d*m))/num         0.80 0.80 0.60; ...  % light dirty
      2/num+(round(2*d*m)-1)/num     0.40 0.40 0.40; ...  % dark gray
      2/num+(round(2*d*m))/num       0.10 1.00 1.00; ...  % cyan blue
      2/num+(round(3*d*m)-1)/num     0.00 0.00 0.75; ...  % dark blue
      2/num+(round(3*d*m))/num       0.00 1.00 0.00; ...  % light green
      2/num+(round(4*d*m)-1)/num     0.00 0.50 0.00; ...  % dark green
      2/num+(round(4*d*m))/num       1.00 1.00 0.00; ...  % yellow
      2/num+(round(5*d*m)-1)/num     1.00 0.50 0.00; ...  % orange
      2/num+(round(5*d*m))/num       1.00 0.00 0.00; ...  % torch red
      2/num+(round(6*d*m)-1)/num     0.50 0.00 0.00; ...  % dark red
      2/num+(round(6*d*m))/num       1.00 0.00 1.00; ...  % magenta
      num/(num+1)                    0.60 0.30 0.80; ...  % purple
      1.000                          1.00 1.00 1.00];
x = feval('fleximap',num+1,pt);
return


function [x] = ezmap(num)
if nargin<1,num = 23; end
if num<16, fprintf('Cannot generate map with %d elements\n',num); num = 16; end
num = num-1;
c = (num-1)/num;
d = c/(7-4/num);
m = num-1;
pt = [0.000                          0.00 0.00 0.00; ...  % black
      0.20                           0.25 0.20 0.30; ...
      0.30                           0.35 0.30 0.43; ...
      2/num+(round(3*d*m)-1)/num     0.65 0.65 0.65; ...  % light gray
      2/num+(round(3*d*m))/num       0.00 1.00 0.00; ...  % light green
      2/num+(round(4*d*m)-1)/num     0.00 0.50 0.00; ...  % dark green
      2/num+(round(4*d*m))/num       1.00 1.00 0.00; ...  % yellow
      2/num+(round(5*d*m)-1)/num     1.00 0.50 0.00; ...  % orange
      2/num+(round(5*d*m))/num       1.00 0.00 0.00; ...  % torch red
      2/num+(round(6*d*m)-1)/num     0.50 0.00 0.00; ...  % dark red
      2/num+(round(6*d*m))/num       1.00 0.00 1.00; ...  % magenta
      num/(num+1)                    0.60 0.33 0.75; ...  % purple
      1.000                          1.00 1.00 1.00];
x = feval('fleximap',num+1,pt);
return


function [x] = ezmap2(num)
if nargin<1,num = 23; end
if num<16, fprintf('Cannot generate map with %d elements\n',num); num = 16; end
num = num-1;
c = (num-1)/num;
d = c/(7-4/num);
m = num-1;
pt = [0.000                          0.00 0.00 0.00; ...  % black
      0.20                           0.25 0.20 0.30; ...
      2/num+(round(2*d*m)-1)/num     0.35 0.30 0.43; ...
      2/num+(round(2*d*m))/num       0.10 1.00 1.00; ...  % cyan blue
      2/num+(round(3*d*m)-1)/num     0.00 0.00 0.75; ...  % dark blue
      2/num+(round(3*d*m))/num       0.00 1.00 0.00; ...  % light green
      2/num+(round(4*d*m)-1)/num     0.00 0.50 0.00; ...  % dark green
      2/num+(round(4*d*m))/num       1.00 1.00 0.00; ...  % yellow
      2/num+(round(5*d*m)-1)/num     1.00 0.50 0.00; ...  % orange
      2/num+(round(5*d*m))/num       1.00 0.00 0.00; ...  % torch red
      2/num+(round(6*d*m)-1)/num     0.50 0.00 0.00; ...  % dark red
      2/num+(round(6*d*m))/num       1.00 0.00 1.00; ...  % magenta
      num/(num+1)                    0.60 0.33 0.75; ...  % purple
      1.000                          1.00 1.00 1.00];
x = feval('fleximap',num+1,pt);
return


function [x] = gomap(num)
if nargin < 1,num = size(colormap,1); end
pt = [0.00 0.25 0.10 0.00;...
      0.20 0.85 0.45 0.05;...
      0.50 1.00 1.00 0.65;...
      0.80 0.00 0.85 0.00;...
      1.00 0.00 0.25 0.00];
x = feval('fleximap',num,pt);
return


function [x] = ogmap(num)
if nargin < 1,num = size(colormap,1); end
pt = [0    0.00 0.25 0.00;...
      0.20 0.00 0.85 0.00;...
      0.50 1.00 1.00 1.00;...
      0.80 0.85 0.45 0.05;...
      1.00 0.25 0.10 0.00];
x = feval('fleximap',num,pt);
return


function [x] = refmap(num)
if nargin<1,num = size(colormap,1); end
pt = [0.00 1.0 0.5 0.0;...  % 
      0.45 1.0 1.0 0.0;...  %    orange
      0.50 1.0 1.0 1.0;...  %    white
      0.55 0.0 1.0 0.0;...  %    green
      1.00 0.0 0.5 0.0];    % 
x = feval('fleximap',num,pt);
return 


function [x] = rbmap(num)
if nargin < 1,num = size(colormap,1); end
pt = [0    0.0 0.0 0.5;...
      0.20 0.0 0.0 1.0;...
      0.50 1.0 1.0 1.0;...
      0.80 1.0 0.0 0.0;...
      1.00 0.5 0.0 0.0];
x = feval('fleximap',num,pt);
return


function [x] = rgmap(num)
if nargin<1,num = size(colormap,1); end
pt = [0    0.00 0.20 0.00;...
      0.30 0.00 0.80 0.00;...
      0.50 0.85 0.85 0.85;...
      0.70 0.80 0.00 0.00;...
      1.00 0.20 0.00 0.00];
x = feval('fleximap',num,pt);
return


function [x] = rgmapinv(num)
if nargin<1,num = size(colormap,1); end
x = feval('rgmap',num);
idx1 = floor(num/2); idx2 = num-idx1+1;
x = x([idx1:-1:1 idx1+1:idx2-1 num:-1:idx2],:);
return


function [x] = rgmap2(num)
if nargin<1,num = size(colormap,1); end
c = floor(num/2);
pt = [0             0.00 1.00 0.00;...
      (c-2)/(num-1) 0.00 0.40 0.00;...
      (c-1)/(num-1) 0.22 0.33 0.22;...
      (c)/(num-1)   0.40 0.40 0.40;...
      (c+1)/(num-1) 0.33 0.22 0.22;...
      (c+2)/(num-1) 0.45 0.00 0.00;...
      1.000         1.00 0.00 0.00];
x = feval('fleximap',num,pt);
return


function [x] = rgmapw(num)
if nargin<1,num = size(colormap,1); end
c = floor(num/2);
pt = [0             0.00 1.00 0.00;...
      (c-2)/(num-1) 0.00 0.40 0.00;...
      (c-1)/(num-1) 0.45 0.58 0.45;...
      (c)/(num-1)   0.58 0.45 0.45;...
      (c+1)/(num-1) 0.45 0.00 0.00;...
      1.000         1.00 0.00 0.00];
x = feval('fleximap',num,pt);
return


function [x] = rgmapn(num)
if nargin<1,num = size(colormap,1); end
c = floor(num/2);
pt = [0             0.00 1.00 0.00;...
      (c-2)/(num-1) 0.00 0.45 0.00;...
      (c-1)/(num-1) 0.15 0.27 0.15;...
      (c)/(num-1)   0.27 0.15 0.15;...
      (c+1)/(num-1) 0.45 0.00 0.00;...
      1.000         1.00 0.00 0.00];
x = feval('fleximap',num,pt);
return


function [x] = t7map(num)
if nargin<1,num = size(colormap,1); end
pt = [0     0.00 0.20 0.30;...
      0.350 0.00 0.47 0.35;...
      0.699 0.00 0.75 0.20;...
      0.700 1.00 0.25 1.00;...
      0.950 1.00 1.00 0.00;...
      1.000 1.00 1.00 1.00];
x = feval('fleximap',num,pt);
return


function [x] = t6map(num)
if nargin<1,num = size(colormap,1); end
pt = [0     0.00 0.10 0.25;...
      0.300 0.00 0.32 0.38;...
      0.599 0.00 0.65 0.20;...
      0.600 1.00 0.30 0.90;...
      0.920 1.00 1.00 0.00;...
      1.000 1.00 1.00 1.00];
x = feval('fleximap',num,pt);
return


function [x] = t5map(num)
if nargin<1,num = size(colormap,1); end
pt = [0     0.00 0.10 0.25;...
      0.250 0.00 0.32 0.38;...
      0.499 0.00 0.65 0.20;...
      0.500 1.00 0.20 1.00;...
      0.850 1.00 1.00 0.00;...
      1.000 1.00 1.00 1.00];
x = feval('fleximap',num,pt);
return


function [x] = bgraymap(num)
if nargin<1,num = size(colormap,1); end
pt = [0.000 0.00 0.90 0.10
      0.499 1.00 1.00 1.00;...
      0.500 0.50 0.50 0.50;...
      0.501 0.00 0.00 0.00;...
      1.000 1.00 0.10 0.10];
x = feval('fleximap',num,pt);
return


function [x] = wmap(num)
if nargin<1,num = size(colormap,1); end
pt = [0.00 1.0 1.0 1.0;...  % 
      0.25 0.0 1.0 1.0;...  %    cyan
      0.80 0.0 0.0 1.0;...  %    blue
      1.00 0.0 0.0 0.5];    % 
x = feval('fleximap',num,pt);
return 


function [x] = tempmap(num)
if nargin<1,num = size(colormap,1); end
pt = [0.000 0.0 0.0 0.8;...  %    blue
	  0.100 0.2 0.5 1.0;...  %    mid blue
      0.200 0.0 1.0 1.0;...  %    fluoroscent cyan
      0.333 0.0 0.7 0.1;...  %    forest green
      0.467 0.5 1.0 0.0;...  %    lemon green
      0.533 1.0 1.0 0.0;...  %    pure yellow
      0.700 1.0 0.3 0.0;...  %    red-orange
      0.733 1.0 0.0 0.0;...  %    husker red
      0.900 0.3 0.0 0.0;...  %    devil's red
      1.000 0.5 0.0 0.6];    %    
x = feval('fleximap',num,pt);
return 


function [x] = tempmapinv(num)
if nargin<1,num = size(colormap,1); end
x = tempmap(num);
x = flipud(x);
return 


function [x] = brmap(num)
if nargin<1,num = size(colormap,1); end
midcolor = 0.7*[1 1 1];
c = 2;
pt = [0.00 0.0 0.0 0.4;...  % Dark blue
      (0.5*num-c-0.01)/num 0.2 0.5 1.0;...  % Light blue
      (0.5*num-c)/num midcolor;...  % Gray
      (0.5*num+c)/num midcolor;...  % Gray
      (0.5*num+c+0.01)/num 1.0 0.4 0.4;...  % Light red
      1.00 0.5 0.0 0.0];    % Dark red
x = feval('fleximap',num,pt);
return


function [x] = rcmap(num)
if nargin<1,num = size(colormap,1); end
pt = [0.00 0.7 0.0 0.3;...  % 
      0.25 1.0 0.5 0.0;...  %    orange
      0.45 1.0 1.0 0.0;...  %    yellow
      0.50 1.0 1.0 1.0;...  %    white
      0.58 0.0 1.0 0.0;...  %    green
      0.78 0.0 0.5 0.3;...  %    dark green
      1.00 0.0 1.0 1.0];    %    cyan
x = feval('fleximap',num,pt);
return 


function [x] = mapinhexq(map)
if nargin<1, map = rgmap(); end
if size(map,2)~=3
	fprintf('Colormap must be in N x 3');
	x = [];
	return;
end
x = reshape(dec2hex(round(255*map'))',[6 size(map,1)])';
idx = 1;
while(idx<=size(map,1))
	if rem(idx,6)==1
		fprintf('\n');
	end
	fprintf('"%s",',x(idx,1:6));
	idx=idx+1;
end
fprintf('\n');
return; 


function [x] = mapinhex(map)
if nargin<1, map = rgmap(); end
if size(map,2)~=3
	fprintf('Colormap must be in N x 3');
	x = [];
	return;
end
x = reshape(dec2hex(round(255*map'))',[6 size(map,1)])';
idx = 1;
while(idx<=size(map,1))
	if rem(idx,6)==1
		fprintf('\n');
	end
	fprintf('0x%s,',x(idx,1:6));
	idx=idx+1;
end
fprintf('\n');
return; 


function [x] = rapmap(npts,mode)
if nargin < 1,npts = size(colormap,1); mode='lin';
elseif nargin < 2,mode='lin'; end
cromIQ = [0.00, 0.00; ...
         -0.07, 0.13; ...
         -0.11, 0.20; ...
         -0.03, 0.21; ...
          0.23, 0.19; ...
          0.39, 0.11; ...
          0.41,-0.01; ...
          0.39,-0.12; ...
          0.32,-0.19; ...
          0.20,-0.15; ...
          0.00, 0.00];
oldx = 1:11;
newx = linspace(1,11,npts).';
switch (mode)
    case 'lin'
        yiqmap(1:npts,1) = linspace(0,1,npts).';
    case 'exu'
        tmp = linspace(0,1,npts).';
        yiqmap(1:npts,1) = 0.45*tmp+0.55*(exp(1)-exp(1-tmp))/(exp(1)-1);
    case 'exd'
        tmp = linspace(0,1,npts).';
        yiqmap(1:npts,1) = 0.45*tmp+0.55*(exp(tmp)-1)/(exp(1)-1);
    case 'rap',
        Raw = [0.00 0.00 0.00; ...
               0.15 0.15 0.50; ...
               0.30 0.15 0.75; ...
               0.60 0.20 0.50; ...
               1.00 0.25 0.15; ...
               0.90 0.50 0.00; ...
               0.90 0.75 0.10; ...
               0.90 0.90 0.50; ...
               1.00 1.00 1.00];
        newx = linspace(1,9,npts).';
        oldx = 1:9;
        cromIQ = Raw*[0.596,-0.274,-0.322; ...
                      0.211,-0.523, 0.312].';
        yiqmap(1:npts,1) = interp1(oldx,Raw*[0.299,0.587,0.114].',newx,'cubic');
end
yiqmap(1:npts,2) = interp1(oldx,cromIQ(:,1),newx,'cubic');
yiqmap(1:npts,3) = interp1(oldx,cromIQ(:,2),newx,'cubic');
x=yiqmap*[1.0000, 1.0000, 1.0000; ...
          0.9562,-0.2727,-1.1037; ...
          0.6214,-0.6468, 1.7006];
x(x<0)=0; x(x>1)=1;


function [x] = asymap0(npts,lims)
if nargin < 1,npts = 64; end
if nargin < 2,lims = [-9,+4,-9,+5]; end
if length(lims) < 2,fprintf('LIMS must be at least a vector of length 2.\n'); x = []; return; end
if length(lims) < 4,lims = [lims,-9,+5]; end
[minnum,maxnum,minclip,maxclip] = deal(lims(1),lims(2),lims(3),lims(4));
minnum = max(minnum,minclip); maxnum = min(maxnum,maxclip);
rang = maxnum - minnum;
if minnum > 0,fprintf('minnum should be negative.\n'); x = []; return; end
if minclip > 0,fprintf('minclip should be negative.\n'); x = []; return; end
if maxnum < 0,fprintf('maxnum should be positive.\n'); x = []; return; end
if maxclip < 0,fprintf('maxclip should be positive.\n'); x = []; return; end
map_num = maxclip - minclip;
long_map_len = ceil(npts*map_num/rang);
num_clr = ceil(-minclip/map_num*long_map_len);
end_sec = floor(num_clr/4);
clr1 = linspace(0,0.05,num_clr-3*end_sec); clr2 = linspace(0.05,0.4,2*end_sec+1); clr3 = linspace(0.4,1,end_sec+1);
rcurv = [clr1,clr2(2:end),clr3(2:end)];
clr1 = linspace(0,0.20,end_sec); clr2 = linspace(0.20,1,num_clr-2*end_sec+1); clr3 = ones(1,end_sec);
gcurv = [clr1,clr2(2:end),clr3];
clr1 = linspace(0.35,0.7,end_sec); clr2 = linspace(0.7,1,num_clr-end_sec+1);
bcurv = [clr1,clr2(2:end)];
tmp1 = [rcurv',gcurv',bcurv'];
tmp2 = feval('boonlib','rbmap',ceil(2*maxclip/map_num*long_map_len));
tmp2 = tmp2(ceil(size(tmp2,1)/2)+1:end,:);
asymap = [tmp1; tmp2];
stidx = round(long_map_len*(minnum-minclip)/(maxclip-minclip))+1;
x = asymap(stidx:stidx+npts-1,:);
return


function [x] = asymap(num,lims)
if nargin < 1,num = size(colormap,1); end
if nargin < 2,lims = [-9,+4,-9,+5]; end
if length(lims) < 2,fprintf('LIMS must be at least a vector of length 2.\n'); x = []; return; end
if length(lims) < 4,lims = [lims,lims]; end
if lims(1)<lims(3),lims(3) = lims(1); end
if lims(2)>lims(4),lims(4) = lims(2); end
tmp = -lims(3)/(lims(4)-lims(3));
pt = [0           0.00 0.00 0.45;...
      0.25*tmp    0.00 0.30 0.75;...
      0.75*tmp    0.40 1.00 0.90;...
      tmp         1.00 1.00 1.00;...
      0.6+0.4*tmp 1.00 0.00 0.00;...
      1.00        0.50 0.00 0.00];
x = feval('fleximap',ceil(num*(lims(4)-lims(3))/(lims(2)-lims(1))),pt);
stidx = max(floor(num*(lims(1)-lims(3))/(lims(2)-lims(1)))+1,1);
x = x(stidx:stidx+num-1,:);
return


function [x] = zmap(zlim)
if nargin<1,zlim = [0 75]; end
if numel(zlim)==1, zlim = [0 5*zlim]; end
if zlim(1)<0,zlim(1) = 0; end
if zlim(2)>75,zlim(2) = 75; end
%tmp = linspace(zlim(1),zlim(2),16);
tmp = zlim(1):5:zlim(2);
tmp = floor(tmp/5)+1;
std_map = [mybgcolor(); ...
           0.20 1.00 1.00; ...
           0.20 0.60 1.00; ...
           0.00 0.00 1.00; ...
           0.30 1.00 0.00; ...
           0.10 0.80 0.00; ...
           0.00 0.60 0.00; ...
           1.00 1.00 0.00; ...
           1.00 0.75 0.00; ...
           1.00 0.50 0.00; ...
           1.00 0.00 0.00; ...
           0.75 0.00 0.00; ...
           0.50 0.00 0.00; ...
           1.00 0.00 0.80; ...
           0.60 0.30 1.00; ...
           1.00 1.00 1.00];
x = std_map(tmp,:);
return




function [x] = zmap2(zlim)
if nargin<1,zlim = [7.5 72.5]; end
if zlim(1)<7.5,zlim(1) = 7.5; end
if zlim(2)>72.5,zlim(2) = 72.5; end
tmp = zlim(1):5:zlim(2);
tmp = floor(tmp/5);
std_map = [  0    0    0; ...
             5  143  143; ...
           128  225   80; ...
            99  185   63; ...
            72  143   48; ...
            44  104   33; ...
            16   63   15; ...
           241  191   16; ...
           240  127   33; ...
           240   15   33; ...
           143    1    2; ...
           177   32  127; ...
           202   63  161; ...
           255  255  255]/255;
x = std_map(tmp,:);
return


function [x] = zmap3(num)
if nargin<1,num = 26; end
if num<12, fprintf('Cannot generate map with %d elements\n',num); num = 12; end
num = num-1;
% pt = [0.00   0.00 0.00 0.00; ...
%       0.10   0.25 0.20 0.30; ...
%       0.20   0.35 0.30 0.43; ...
%       0.30   0.60 0.60 0.60; ...
%       0.35   0.00 0.75 0.00; ...
%       0.40   0.00 0.40 0.00; ...
%       0.50   0.80 0.80 0.00; ...
%       0.55   0.70 0.70 0.00; ...
% %       0.73   0.70 0.30 0.00; ...
%       0.73   1.00 0.00 0.00; ...
%       0.83   0.55 0.00 0.00; ...
%       0.90   0.70 0.00 0.70; ...
%       0.93   0.60 0.20 0.80; ...
%       1.00   1.00 1.00 1.00];
pt = [0.00   0.00 0.00 0.00; ...
%       0.10   0.25 0.20 0.30; ...
%       0.20   0.35 0.30 0.43; ...
      0.10   0.60 0.60 0.60; ...
      0.20   0.00 0.75 0.00; ...
      0.30   0.00 0.40 0.00; ...
      0.40   0.80 0.80 0.00; ...
      0.50   0.70 0.70 0.00; ...
      0.60   1.00 0.00 0.00; ...
      0.70   0.90 0.00 0.00; ...
      0.75   0.55 0.00 0.00; ...
      0.85   0.70 0.00 0.70; ...
      0.90   0.60 0.20 0.80; ...
      1.00   1.00 1.00 1.00];
x = feval('fleximap',num+1,pt);
return


function [x] = zmapn(zlim)
if nargin<1,zlim = [-5 75]; end
if numel(zlim)==1, zlim = [-5 5*zlim]; end
if zlim(1)<-5,zlim(1) = -5; end
if zlim(2)>75,zlim(2) = 75; end
%tmp = linspace(zlim(1),zlim(2),16);
tmp = zlim(1):5:zlim(2);
tmp = floor((tmp-zlim(1))/5)+1;
std_map = [0.00 0.00 0.00; ...
           0.20 0.30 0.30; ...
           0.20 1.00 1.00; ...
           0.20 0.60 1.00; ...
           0.00 0.00 1.00; ...
           0.30 1.00 0.00; ...
           0.10 0.80 0.00; ...
           0.00 0.60 0.00; ...
           1.00 1.00 0.00; ...
           1.00 0.75 0.00; ...
           1.00 0.50 0.00; ...
           1.00 0.00 0.00; ...
           0.75 0.00 0.00; ...
           0.50 0.00 0.00; ...
           1.00 0.00 0.80; ...
           0.60 0.30 1.00; ...
           1.00 1.00 1.00];
x = std_map(tmp,:);
return


function [x] = zmapnx(varargin)
	tmp = zmapn(varargin{:});
	x = [0.03 0.06 0.06; ...
	     0.05 0.08 0.08; ...
	     0.10 0.15 0.15; ...
	     tmp(2:end,:)];
return


function [x] = vrmap(lim)
if nargin<1,lim = [-32 32]; end
if lim(1)<-32,lim(1) = 32; end
if lim(2)>32,lim(2) = 32; end
tmp = lim(1):4.5:lim(2);
tmp = floor((tmp+36.5)/4.5);
std_map = [ 24  24 181; ...
            76  76 255; ...
             1 179 179; ...
            76 255 255; ...
             0 179   1; ...
            76 255  76; ...
           102 102 102; ...
            10  10  10; ...
           179 179 179; ...
           255 255  76; ...
           179 179   0; ...
           255  76  76; ...
           179   0   1; ...
           255  76 255; ...
           179   0 179]/255;
x = std_map(tmp,:);
return


function [x] = blackngold(num)
if nargin<1,num = size(colormap,1); end
c = floor(num/2);
pt = [0             0.10 0.10 0.10;...
      (c-1)/(num-1) 0.20 0.20 0.20;...
      (c)/(num-1)   1.00 0.90 0.00;...
      0.83          1.00 1.00 0.00;...
      1.000         1.00 1.00 1.00];
x = feval('fleximap',num,pt);
return


% The engine behind all the color maps
function [cmap] = fleximap(num,pt)
if nargin < 1,
    num = 15;
    pt = [0    0.5 0.0 0.0;...
          0.20 1.0 0.0 0.0;...
          0.50 1.0 1.0 1.0;...
          0.80 0.0 0.0 1.0;...
          1.00 0.0 0.0 0.5];
end
if any(pt>1),fprintf('PT Matrix cannot have value > 1.0.\n'); return; end
x = linspace(0,1,num); x = x(:);
cmap = interp1(pt(:,1),pt(:,2:4),x,'linear');
cmap(cmap<0) = 0;
cmap(cmap>1) = 1;
return

% ------------------------------ Drawings ------------------------------

function [x] = draw3dbox(lims)
if nargin<1,lims = [0 1 0 1 0 1]; end
% Vertices of each surface, ordered according to right hand rule
nor_sfc_x = [0 0 1 1]';   nor_sfc_y = [1 1 1 1]';    nor_sfc_z = [0 1 1 0]';
sou_sfc_x = [0 1 1 0]';   sou_sfc_y = [0 0 0 0]';    sou_sfc_z = [0 0 1 1]';
eas_sfc_x = [1 1 1 1]';   eas_sfc_y = [0 1 1 0]';    eas_sfc_z = [0 0 1 1]';
wes_sfc_x = [0 0 0 0]';   wes_sfc_y = [1 0 0 1]';    wes_sfc_z = [0 0 1 1]';
top_sfc_x = [0 1 1 0]';   top_sfc_y = [0 0 1 1]';    top_sfc_z = [1 1 1 1]';
bot_sfc_x = [1 1 0 0]';   bot_sfc_y = [1 0 0 1]';    bot_sfc_z = [0 0 0 0]';
% Composite surface
xx = [nor_sfc_x,sou_sfc_x,eas_sfc_x,wes_sfc_x,top_sfc_x,bot_sfc_x];
yy = [nor_sfc_y,sou_sfc_y,eas_sfc_y,wes_sfc_y,top_sfc_y,bot_sfc_y];
zz = [nor_sfc_z,sou_sfc_z,eas_sfc_z,wes_sfc_z,top_sfc_z,bot_sfc_z];
xx = xx*(lims(2)-lims(1))+lims(1);
yy = yy*(lims(4)-lims(3))+lims(3);
zz = zz*(lims(6)-lims(5))+lims(5);
x = patch(xx,yy,zz,'b');
return


function [x] = drawarrow(lims)
if nargin<1,lims = [0 0 0.5]; end
if length(lims)<4,t = pi/4; else t = lims(4); end
if length(lims)<5,a = 1; else a = lims(5); end
ha = findobj(gcf,'UserData','ArrowOverlay');
orig_unit = get(gcf,'Unit');
figsize = get(gcf,'Position');
if isempty(ha)
    ha = axes('Unit','Pixel','Position',[1 1 figsize(3:4)]);
elseif length(ha)>1
    delete(ha)
    ha = axes('Unit','Pixel','Position',[1 1 figsize(3:4)]);
else
    set(ha,'Unit','Pixel'); hold on
end
axsize = get(ha,'Position');
axaspect = axsize(3)/axsize(4);
a = a/50;
xx = [0 0 lims(3)-a lims(3)-a lims(3) lims(3)-a lims(3)-a 0 0];
yy = [1 -1 -1 -2 0 2 1 1 -1]*0.3*a;
xxp = xx*cos(t)+yy*sin(t)+lims(1);
yyp = (xx*sin(t)-yy*cos(t))*axaspect+lims(2);
tmp = fill(xxp,yyp,'g');
set(ha,'Visible','Off','Unit','Pixel',...
       'XLim',[-0.01 1.01],'XLimMode','Manual',...
       'YLim',[-0.01 1.01],'YLimMode','Manual','UserData','ArrowOverlay');
hold off
set(gcf,'Unit',orig_unit);
x = [ha tmp];
return


function [x] = cspace(cmap)
clf reset
bsizewin(gcf,[1000 500])
len=size(cmap,1);
yiqmap = cmap*[0.299, 0.587, 0.114; ...
               0.596,-0.274,-0.322; ...
               0.211,-0.523, 0.312].';
% Normalize it into [0,1] range
yiqmap(:,2) = (yiqmap(:,2)+0.596)/1.192;
yiqmap(:,3) = (yiqmap(:,3)+0.523)/1.046;

% Count number of shades per color, simple laplacian edge detection
laplacian = [ones(1,6); diff(cmap,1,1) diff(yiqmap,1,1)];
edges = max(abs(laplacian),[],2)>0.2;
epos = find(edges);
epossep = epos(2:end)-epos(1:end-1);
eposmid = 0.5*(epos(2:end)+epos(1:end-1)-1);
cmap2 = cmap;
for ie = 1:numel(epossep)
	cmap2(epos(ie):epos(ie+1)-1,1) = cmap(floor(eposmid(ie)),1);
	cmap2(epos(ie):epos(ie+1)-1,2) = cmap(floor(eposmid(ie)),2);
	cmap2(epos(ie):epos(ie+1)-1,3) = cmap(floor(eposmid(ie)),3);
end

% Modify edges to have 1 = nothing, 2 = start, 3 = end, 4 = both.
edges = edges+1;
edges(epos(1:end-1)+epossep-1) = edges(epos(1:end-1)+epossep-1)+2;

%subplot(2,1,1)
axes('Position',[0.05 0.55 0.9 0.4])
plot(linspace(0,1,len),':','Color',[1 1 1]*0.7);
hold on
x(1)=plot(cmap(:,1),'Color',[1 0 0]);
x(2)=plot(cmap(:,2),'Color',[0 1 0]);
x(3)=plot(cmap(:,3),'Color',[0 0 1]);
x(4)=plot(yiqmap(:,1),'Color',[0.5 0.5 0.5]);
x(5)=plot(yiqmap(:,2),'Color',[0.8 0.8 0.0]);
x(6)=plot(yiqmap(:,3),'Color',[0.4 0.7 1.0]);
% [legh objh outh outm] = legend(x,'Red Channel','Green Channel','Blue Channel',...
%          'Luminance Y','Chrominance I','Chrominance Q',-1);
[legh objh] = legend(x,'Red Channel','Green Channel','Blue Channel',...
	'Luminance Y','Chrominance I','Chrominance Q',-1);
tmp = findobj(objh,'Type','Text');
set(tmp,'Color',get(0,'defaultTextColor'));
tmp = get(0,'defaultAxesXColor');
set(legh,'XColor',tmp,'YColor',tmp,'ZColor',tmp);
hold off
axis([0.5 len+0.5 0 1]);
ylabel('Component Value');
pos=get(gca,'Position');
if len<50,set(x,'Marker','.'); end
set(gca,'XTickLabel',[]);

% subplot(2,1,2)
axes('Position',[0.05 0.1 0.9 0.4])
img=zeros(9,len,3);
img(1,1:len,1:3)=ind2rgb(1:len,cmap);
img(2,1:len,1:3)=ind2rgb(round(yiqmap(:,1).'*255)+1,gray(256));
img(3,1:len,1:3)=ind2rgb(round(yiqmap(:,2).'*255)+1,linspace(0,1,256).'*[1 1 0]);
img(4,1:len,1:3)=ind2rgb(round(yiqmap(:,3).'*255)+1,linspace(0,1,256).'*[0.2 0.8 1]);
img(5,1:len,1:3)=ind2rgb(round(cmap(:,1).'*255)+1,[linspace(0,1,256).',zeros(256,2)]);
img(6,1:len,1:3)=ind2rgb(round(cmap(:,2).'*255)+1,[zeros(256,1),linspace(0,1,256).',zeros(256,1)]);
img(7,1:len,1:3)=ind2rgb(round(cmap(:,3).'*255)+1,[zeros(256,2),linspace(0,1,256).']);
img(8,1:len,1:3)=ind2rgb(edges,[0 0 0; 1 0.3 0.3; 0.2 0.6 1; 1 0.5 1]);
img(9,1:len,1:3)=ind2rgb(1:len,cmap2);

image(img)
tmppos=get(gca,'Position');
set(gca,'Position',[tmppos(1:2),pos(3:4)],'YAxisLocation','Right','YTick',1:9,...
    'YTickLabel',{'Colors','Luminance Y','Chrominance I','Chrominance Q',...
                   'Red Channel','Green Channel','Blue Channel','Shade Counts','Mid Colors'})
hold on
text(eposmid,8*ones(1,numel(eposmid)),num2cell(epossep),'HorizontalAlignment','Center','Color','y')
hold off

xlabel('Colormap Index');
return


function [xx yy] = makecircle(params,npts)
if nargin<1,params = [0 0 1]; end
if nargin<2,npts = 90; end
if size(params,2)~=3
    fprintf('Parameter must be in Nx3,e.g. [0 0 1; 0 0 2]\n')
    return
end
x = params(:,1); y = params(:,2); radii = params(:,3);
ang_rad = linspace(0,2*pi,npts);
xx = radii*cos(ang_rad)+repmat(x,[1 npts]);
yy = radii*sin(ang_rad)+repmat(y,[1 npts]);
xx(:,npts+1) = nan;
yy(:,npts+1) = nan;
xx = xx.'; xx = xx(:);
yy = yy.'; yy = yy(:);
return


function [hcir] = circle(params,npts)
if nargin<1,params = [0 0 1]; end
if nargin<2,npts = 90; end
[xx yy] = makecircle(params,npts);
hold on
hcir = plot(xx,yy,'Color','k');
hold off
return


function [h] = spidergrid(a,sep,nsec)
if nargin<1 || isempty(a), a = gca; end
m = axis; m = max(m(:));
if nargin<2, sep = m/4; end
m = m+2*sep;
if nargin<3, nsec = 12; end
phi = (0:nsec-1)/nsec*2*pi;
b = sin(phi); b = round(b*1e4)*1e-4; % avoid some small misalign
c = cos(phi); c = round(c*1e4)*1e-4;
xx = kron(b,[nan 0 1])'*m;
xy = kron(c,[nan 0 1])'*m;
r = (sep:sep:m).';
n = length(r);
if (n>0)
	[cx cy] = makecircle([zeros(n,2) r]);
	xx = [cx; xx];
	yy = [cy; xy];
end
h = zeros(length(a),1);
for idx = 1:length(a)
	axes(a(idx)); %#ok<LAXES>
	c = get(gca,'Color');
	hold on
	h(idx) = plot(xx,yy,'Color',1-c);
	hold off
end
return;


% ------------------------------ Don't ask ------------------------------
function grayme(h,blend)
if nargin<1, h = gcf; end
if (nargin==1)&&~ishandle(h), blend = h; h = gcf; end
if ~exist('blend','var'), blend = 1; end
if ~ishandle(h), fprintf('Invalid handle object.\n'); return; end
feval('tintme',h,[0.5 0.5],blend);
return


function greenme(h,blend)
if nargin<1, h = gcf; end
if (nargin==1)&&~ishandle(h), blend = h; h = gcf; end
if ~exist('blend','var'), blend = 1; end
if ~ishandle(h), fprintf('Invalid handle object.\n'); return; end
feval('tintme',h,[0.41 0.37],blend);
return


function sepiame(h,blend)
if nargin<1, h = gcf; end
if (nargin==1)&&~ishandle(h), blend = h; h = gcf; end
if ~exist('blend','var'), blend = 1; end
if ~ishandle(h), fprintf('Invalid handle object.\n'); return; end
feval('tintme',h,[0.43 0.53],blend);
return


function tintme(h,adj,blend)
if nargin<1, h = gcf; end
if nargin<2||length(adj)~=2, adj = [0.6 0.4]; end
if ~exist('blend','var'), blend = 1; end
if ~ishandle(h), fprintf('Invalid handle object.\n'); return; end
cmap = get(h,'Colormap');
tmp = rgb2ntsc(cmap);
cmap = rgb2ycbcr(cmap);
cmap(:,1) = tmp(:,1);
cmap(:,2) = blend*adj(1)+(1-blend)*cmap(:,2);
cmap(:,3) = blend*adj(2)+(1-blend)*cmap(:,3);
cmap = ycbcr2rgb(cmap);
clr = get(h,'Color');
clr = rgb2ycbcr(clr); clr(:,2) = 0.4+0.2*adj(1); clr(:,3) = 0.4+0.2*adj(2); clr = ycbcr2rgb(clr);
set(h,'Colormap',cmap,'Color',clr);
%tags = {'Color','BackgroundColor','ForegroundColor','MarkerFaceColor','MarkerEdgeColor','XColor','YColor'};
tags = {'Color','MarkerFaceColor','MarkerEdgeColor','XColor','YColor','FaceColor'};
ha = get(h,'Children');
for idx=1:length(ha)
    for itag=1:length(tags)
        if isprop(ha(idx),tags{itag})
            clr = get(ha(idx),tags{itag});
            if isnumeric(clr),
                %fprintf('Prop = %s  Val = %s\n',tags{itag},num2str(tmp));
				clr = rgb2ycbcr(clr);
				clr(:,2) = blend*adj(1)+(1-blend)*clr(:,2);
				clr(:,3) = blend*adj(2)+(1-blend)*clr(:,3);
				clr = ycbcr2rgb(clr);
                set(ha(idx),tags{itag},clr);
            end
        end
    end
    % Second Level
    hb = get(ha(idx),'Children');
	if strcmp(get(ha(idx),'Type'),'axes')
    	hb = [hb; get(ha(idx),'XLabel');
        get(ha(idx),'YLabel');
        get(ha(idx),'Title')]; %#ok<AGROW>
	end
    for jdx=1:length(hb)
        for itag=1:length(tags)
            if isprop(hb(jdx),tags{itag}),
                clr = get(hb(jdx),tags{itag});
                if isnumeric(clr),
                    %fprintf('Prop = %s  Val = %s\n',tags{itag},num2str(tmp));
					clr = rgb2ycbcr(clr);
					clr(:,2) = blend*adj(1)+(1-blend)*clr(:,2);
					clr(:,3) = blend*adj(2)+(1-blend)*clr(:,3);
					clr = ycbcr2rgb(clr);
                    set(hb(jdx),tags{itag},clr);
                end
            end
        end
    end
end
% Special image type
ha = findobj(h,'Type','Image');
for idx=1:length(ha)
    tmp = get(ha(idx),'CData');
    if size(tmp,3)==3
        % The RGB type
		clr = rgb2ycbcr(tmp);
		clr(:,2) = blend*adj(1)+(1-blend)*clr(:,2);
		clr(:,3) = blend*adj(2)+(1-blend)*clr(:,3);
		clr = ycbcr2rgb(clr);
        set(ha(idx),'CData',clr);
    end
end
return


function [fname] = choosefile(arg1,arg2,arg3)
if nargin<2,arg2 = '*'; end
if nargin<1,arg1 = './'; end
if iscell(arg1), path_choice = arg1;
elseif exist(arg1,'dir')==7, path_choice = {arg1};
elseif exist(arg1,'file')==2, path_choice = feval(arg1);
elseif ~exist(arg1,'dir') && ~exist(arg1,'file'), fprintf('Path or function does not exist.\n'); fname = []; return;
end
if exist('arg3','var') && numel(arg3)~=2, clear arg3; end
if length(path_choice)>1
	if exist('arg3','var')&&((arg3(1)>0)&&(arg3(1)<=length(path_choice)))
		tmp = arg3(1);
		% fprintf('Pre-selected Folder #%d: %s\n',tmp,char(path_choice(tmp)));
	else
		for idx = 1:length(path_choice);
			fprintf(' %3d. %s\n',idx,char(path_choice{idx}));
		end
		fprintf('\n');
		tmp = input(['Selection (1-',num2str(idx),') : ']);
		if isempty(tmp) || ~ismember(tmp,1:idx)
			fprintf('No valid selection,action cancelled !\n\n');
			fname = [];
			return
		end
	end
	dirname = char(path_choice(tmp));
else
	% fprintf('Only one folder.  Automatically selected.\n');
    dirname = char(path_choice);
end
if ~exist(dirname,'dir')
    fprintf('Directory does not exists, double check the path(s).\n');
    fname = [];
    return;
end
if dirname(end)~='/', dirname = [dirname,'/']; end
flist = dir([dirname,arg2]);
% Take out . and ..
names = {flist.name}';
keep_idx = ~ismember(names,{'.','..','.DS_Store'});
flist = flist(keep_idx);
if isempty(flist)
	fprintf('No file matches the criteria.\n');
	fname = [];
	return;
end

if ~exist('arg3','var')
	flist = feval('showflist',flist);
else
	flist = feval('showflist',flist,1);
end

% Selecting the file
if exist('arg3','var')
	if ((arg3(2)>0)&&(arg3(2)<=length(flist)))
		tmp = arg3(2);
	elseif (arg3(2)<=0)
		tmp = max(length(flist)+arg3(2),1);
	else
		%fprintf('%s: Empty filename.\n',mfilename)
		fname = [];
		return;
	end
	fprintf('Pre-selected File #%d: %s\n',tmp,char(flist(tmp)));
else
	idx = length(flist);
	tmp = input(['\nSelection (1-',num2str(idx),') : ']);
	if isempty(tmp) || ~ismember(tmp,1:idx)
		fprintf('No valid selection,action cancelled !\n\n');
		fname = [];
		return
	end
	fprintf('\n');
end

% Construct the full path filename
if strcmp(dirname,'.'),
    fname = char(flist(tmp));
else
    fname = [dirname,char(flist(tmp))];
end
return


function x = showflist(flist,quiet)
if isempty(flist), x = []; return; end
if ~exist('quiet','var'), quiet = 0; end
if ~isfield(flist,'bytes')&&~isfield(flist,'name')
    fprintf('Input must be something returned by dir().\n')
    x = 0;
    return
end
% Collect files, display them.
fbyte = [flist.bytes]';
flist = {flist.name}';

tmp = regexp(flist,'(?<=\D*)\d+','match');
if ~isempty(tmp)
	% tmp = cat(1,tmp{:});   % <-- Array with numbers for each entry
	% Convert them into numbers
	tmp2 = zeros(size(tmp));
	for idx=1:numel(tmp2),
		tmp3 = str2double(char(tmp{idx}));
		tmp2(idx,1:length(tmp3)) = tmp3;
	end
	% Try to sort them in a smart way, using those numbers
	[tmp2 tmp] = sortrows(tmp2); %#ok<ASGLU>
	flist = flist(tmp);
% Try to find indexing scheme, this is just special case for my xpol project
% if ~isempty(strfind(flist{1},'_P'))|~isempty(strfind(flist{1},'_H'))
% 	dd = zeros(size(fbyte,1),1);
% 	for idx = 1:length(dd)
% 		ii = strfind(flist{idx},'_P');
% 		if isempty(ii), ii = strfind(flist{idx},'_H'); end
% 		tmp = flist{idx}(ii+2:length(flist{idx})-4);
% 		keyboard
% 		dd(idx) = str2num(tmp);
% 	end
% 	[dds tmp] = sort(dd);
% 	flist = flist(tmp);
else
 	[flist tmp] = sortrows(flist);
end
if ~quiet
	fbyte = fbyte(tmp);
	if length(fbyte)<1,msize = 0; else msize = mean(fbyte); end
	if msize<1e3,funits = 'B';
	elseif msize>1e3&&msize<1e6,fbyte = fbyte/1024; funits = 'KB';
	elseif msize>1e6&&msize<1e9,fbyte = fbyte/1024/1024; funits = 'MB';
	elseif msize>1e9&&msize<1e12,fbyte = fbyte/1024/1024/1024; funits = 'GB';
	else fprintf('Filesize is rediculously big at the year of 2005.\n'); x = []; return;
	end
	tmp = size(char(flist));
	msize = ceil(log10(max(fbyte+1)));
	ndig = ceil(log10(length(flist)));
	% Compact template
	ftemplate = ['%',num2str(ndig),'d. %',num2str(tmp(2)),'s %',num2str(msize+3),'.2f ',funits,'\n'];
	if ((tmp(2)+ndig+msize+3+5)<39)
		% Assume 79 screen width
		% ndig digits + 4 spaces; 
		% tmp(2) = filename's length; msize = filesize digits
		leftover = 79-2*(ndig+4+tmp(2)+msize+3+length(funits));
		pad = char(repmat(' ',[1 max(floor(leftover/3),1)]));
		ftemplate = [pad,ftemplate];
	else
		% Can't squeeze the info in two columns ==> space them out a bit
		ftemplate = ['%',num2str(ndig),'d.  %',num2str(tmp(2)),'s   %',num2str(msize+3),'.2f ',funits,'\n'];
		leftover = 79-(ndig+7+tmp(2)+msize+3+length(funits));
		pad = char(repmat(' ',[1 floor(leftover/2)]));
		ftemplate = [pad,ftemplate];
	end
	if length(flist)==1
		idx = 1;
		fprintf(ftemplate,idx,char(flist(idx)),fbyte(idx)) %#ok<PRTCAL>
	else
		if ((tmp(2)+ndig+msize+3+5)<39)
			half_idx = round(length(flist)/2);
			for idx=1:half_idx-1+double(mod(length(flist),2)==0)
				fprintf(ftemplate(1:end-2),idx,char(flist(idx)),fbyte(idx)); % Left column
				fprintf(ftemplate,half_idx+idx,char(flist(half_idx+idx)),fbyte(half_idx+idx)); % Right column
			end
			if mod(length(flist),2)~=0,
				fprintf(ftemplate,half_idx,char(flist(half_idx)),fbyte(half_idx));
			end
		else
			for idx=1:length(flist)
				fprintf(ftemplate,idx,char(flist(idx)),fbyte(idx));
			end
		end
	end
end
x = flist;
return


function [yesno] = aintersect(fh)
yesno = false;
if length(fh)<2,return; end
if nargin<1
    pos_set = get(get(0,'Children'),'Position');
    pos_set = reshape([pos_set{:}],4,length(pos_set)).';
else
    pos_set = get(fh,'Position');
    pos_set = reshape([pos_set{:}],4,length(pos_set)).';
end
xbuff = 1; ybuff = 1;
for iset = 1:size(pos_set,1)
    for iobj = 1:size(pos_set,1)
        if iset~=iobj
            tmp = [((pos_set(iobj,1)>=pos_set(iset,1)&...
                     pos_set(iobj,1)<=pos_set(iset,1)+pos_set(iset,3))&...
                    (pos_set(iobj,2)>=pos_set(iset,2)&...
                     pos_set(iobj,2)<=pos_set(iset,2)+pos_set(iset,4))),...
                   ((pos_set(iobj,1)+pos_set(iobj,3)+xbuff>=pos_set(iset,1)&...
                     pos_set(iobj,1)+pos_set(iobj,3)+xbuff<=pos_set(iset,1)+pos_set(iset,3))&...
                    (pos_set(iobj,2)>=pos_set(iset,2)&...
                     pos_set(iobj,2)<=pos_set(iset,2)+pos_set(iset,4))),...
                   ((pos_set(iobj,1)>=pos_set(iset,1)&...
                     pos_set(iobj,1)<=pos_set(iset,1)+pos_set(iset,3))&...
                    (pos_set(iobj,2)+pos_set(iobj,4)+ybuff>=pos_set(iset,2)&...
                     pos_set(iobj,2)+pos_set(iobj,4)+ybuff<=pos_set(iset,2)+pos_set(iset,4))),...  
                   ((pos_set(iobj,1)+pos_set(iobj,3)+xbuff>=pos_set(iset,1)&...
                     pos_set(iobj,1)+pos_set(iobj,3)+xbuff<=pos_set(iset,1)+pos_set(iset,3))&...
                    (pos_set(iobj,2)+pos_set(iobj,4)+ybuff>=pos_set(iset,2)&...
                     pos_set(iobj,2)+pos_set(iobj,4)+ybuff<=pos_set(iset,2)+pos_set(iset,4)))];
            if any(tmp), yesno = true; return; end
%         else
%             tmp = logical([0 0 0 0]);
        end
    end
end
return


function [x] = assignfig(tag,forcenew)
if nargin<2, forcenew=0; end
if nargin<1, fprintf('I cannot comply, no Tag for search.\n'); return; end
if ~ischar(tag), fprintf('Tag must be a string that was assigned to a figure.\n'); return; end
stidx = 3000;
interv = 20;
if (forcenew>1)
	% If a specific number was given, then just use that as figure handle
	x = forcenew;
	return
end
% Find the object with the match tag
tmp = findobj('Tag',tag);
if ~isempty(tmp), tmp = tmp(ishandle(tmp)); end
if ~isempty(tmp)
    % There is at least 1 Figure
	if (forcenew==1)
        % If forcenew, try to open up a new figure that does not overlap
        ifig = round(stidx+interv*rand(1));
        while ismember(ifig,tmp), ifig = ifig+1; end
        x = ifig;
        return
	elseif (forcenew==0)
        % Don't open up new figure, just use to existing one
        x = tmp(1);
        return
	end
else
    % If figure doesn't exist, just give any handle
    x = round(stidx+interv*rand(1));
    return
end


function [yn] = isoverlap(L,R)
% L = (xmin xmax ymin ymax)
if numel(L)~=4
    fprintf('The first input must be a 4-element vector only.\n');
    yn = -1;
    return
end
if size(R,1)~=4,
    if numel(R)==4,
        R = R(:);
    else
        fprintf('Sorry I can''t compute.\n');
        yn = -1;
        return
    end
end
yn = ( (L(1)>R(1,:))&(L(1)<R(2,:))&...
       (L(3)>R(3,:))&(L(3)<R(4,:)) )|...
     ( (L(1)>R(1,:))&(L(1)<R(2,:))&...
       (L(4)>R(3,:))&(L(4)<R(4,:)) )|...
     ( (L(2)>R(1,:))&(L(2)<R(2,:))&...
       (L(3)>R(3,:))&(L(3)<R(4,:)) )|...
     ( (L(2)>R(1,:))&(L(2)<R(2,:))&...
       (L(4)>R(3,:))&(L(4)<R(4,:)) )|...
     ( (R(1,:)>L(1))&(R(1,:)<L(2))&...
       (R(3,:)>L(3))&(R(3,:)<L(4)) )|...
     ( (R(1,:)>L(1))&(R(1,:)<L(2))&...
       (R(4,:)>L(3))&(R(4,:)<L(4)) )|...
     ( (R(2,:)>L(1))&(R(2,:)<L(2))&...
       (R(3,:)>L(3))&(R(3,:)<L(4)) )|...
     ( (R(2,:)>L(1))&(R(2,:)<L(2))&...
       (R(4,:)>L(3))&(R(4,:)<L(4)) );
return


function mydefault()
set(0,'defaultFigureColor',[0.15 0.15 0.15],'defaultAxesColor',[0 0 0],...
      'defaultAxesXColor',[1 1 1],'defaultAxesYColor',[1 1 1],'defaultAxesZColor',[1 1 1],...
      'defaultAxesYDir','Normal',...
      'defaultTextColor',[1 1 1],...
      'defaultTextFontSize',12,'defaultAxesFontSize',10,...
      'defaultAxesColorOrder',[1 0 0; 0 0.8 0; 0 0.3 1; 1 0.6 0; 0 0.8 0.7; 0.6 0.2 1; 0.6 0.6 0.6],...
      'defaultFigureColormap',rapmap(64),'defaultFigureInvertHardCopy','off');


function mydefault2()
set(0,'defaultFigureColor',[1 1 1],'defaultAxesColor',[1 1 1],...
      'defaultAxesXColor',[0 0 0],'defaultAxesYColor',[0 0 0],'defaultAxesZColor',[0 0 0],...
      'defaultAxesYDir','Normal',...
      'defaultTextColor',[0 0 0],...
      'defaultTextFontSize',12,'defaultAxesFontSize',10,...
      'defaultAxesColorOrder',[1 0 0; 0 0.8 0; 0 0.3 1; 1 0.6 0; 0 0.8 0.7; 0.6 0.2 1; 0.6 0.6 0.6],...
      'defaultFigureColormap',rapmap(64));