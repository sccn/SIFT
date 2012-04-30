function h = hlp_plotAtlas(Surface,hax,plotColors,color,theme)
% Author: Tim Mullen and Alejandro Ojeda, 2012, SCCN/INC UCSD

if nargin<2
    plotColors = true;
end
if nargin<3
    hax = gca;
end
if nargin<4
    color = [1,.75,.65];
end
if nargin<5
    % a theme from hlp_getBrainMovieTheme()
    theme = {};
end

% color = 'jet';

if iscell(color)
    color = color{1};
    % color is a standard colormap name (e.g. 'jet','hsv',...)
    % convert to N x 3 colortable
    eval(sprintf('color = %s(%d);',color,size(Surface.colortable.table,1)));
end

if size(color,1)>1 && size(color,2)==3
    % color is an N x 3 colortable
    if size(color,1)~=size(Surface.colortable.table,1)
        error('SIFT:hlp_plotAtlas','number of rows in colortable must equal number of rows in Surface.colortable.table');
    end
    Surface.colortable.table(:,1:3) = color;
    
end

if plotColors && isfield(Surface,'label')
    colors = zeros(length(Surface.label),3);

    for it=1:size(Surface.colortable.table,1)
        I = Surface.label == Surface.colortable.table(it,5);
        colors(I,1) = Surface.colortable.table(it,1);
        colors(I,2) = Surface.colortable.table(it,2);
        colors(I,3) = Surface.colortable.table(it,3);
    end
    colors = colors./(sqrt(sum(colors.^2,2))*[1 1 1]);
    
    h = patch('vertices',Surface.vertices,'faces',Surface.faces, ...
              'LineStyle','none','parent',hax,'FaceVertexCdata',colors,'facecolor','flat','edgecolor','flat',theme{:}); 
else
    % don't plot individual colors
    
    h = patch('vertices',Surface.vertices,'faces',Surface.faces,'facecolor', color, ...
              'LineStyle','none','parent',hax,theme{:});
end

% userData.h = [];
% userData.label = Surface.label;
% userData.colortable = Surface.colortable;

% colordef black

% axis off;