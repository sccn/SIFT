function [legend_h,object_h,plot_h,text_strings] = columnlegend(numcolumns, str, location)
%
%   columnlegend creates a legend with a specified number of columns.
%   
%   columnlegend(numcolumns, str, location)
%       numcolumns - number of columns in the legend
%       cellstr - cell array of strings for the legend
%       location - location variable for legend, default is 'NorthEast'
%                  possible values: 'NorthWest', 'NorthEast', 'SouthEast', 'SouthWest' 
%
%   example:
%      legend_str = []; 
%      for i=1:10, 
%           x = 1:i:(10*i); 
%           plot(x); hold on; 
%           legend_str = [legend_str; {num2str(i)}];
%      end
%      columnlegend(3, legend_str, 'NorthWest');
%
%
%   Author: Simon Henin <shenin@gc.cuny.edu>


if nargin < 3
   location = 'NorthEast'; 
end


%create the legend
[legend_h,object_h,plot_h,text_strings] = legend(str);

%some variables
numlines = length(str);
numpercolumn = ceil(numlines/numcolumns);

%get old width, new width and scale factor
pos = get(legend_h, 'position');
width = numcolumns*pos(3);
rescale = pos(3)/width;

%get some old values so we can scale everything later
xdata = get(object_h(numlines+1), 'xdata'); 
ydata1 = get(object_h(numlines+1), 'ydata');
ydata2 = get(object_h(numlines+3), 'ydata');

%we'll use these later to align things appropriately
sheight = ydata1(1)-ydata2(1);                  % height between data lines
height = ydata1(1);                             % height of the box. Used to top margin offset
line_width = (xdata(2)-xdata(1))*rescale;   % rescaled linewidth to match original
spacer = xdata(1)*rescale;                    % rescaled spacer used for margins


%put the legend on the upper left corner to make initial adjustments easier
loci = get(gca, 'position');
set(legend_h, 'position', [loci(1) pos(2) width pos(4)]);


col = -1;
for i=1:numlines,
    if mod(i,numpercolumn)==1,
        col = col+1;
    end
    
    if i==1
        linenum = i+numlines;
    else
        linenum = linenum+2;
    end
    labelnum = i;
    
    position = mod(i,numpercolumn);
    if position == 0,
         position = numpercolumn;
    end
    
    %realign the labels
    set(object_h(linenum), 'ydata', [(height-(position-1)*sheight) (height-(position-1)*sheight)]);
    set(object_h(linenum), 'xdata', [col/numcolumns+spacer col/numcolumns+spacer+line_width]);
    
    set(object_h(linenum+1), 'ydata', [height-(position-1)*sheight height-(position-1)*sheight]);
    set(object_h(linenum+1), 'xdata', [col/numcolumns+spacer*3.5 col/numcolumns+spacer*3.5]);
    
    set(object_h(labelnum), 'position', [col/numcolumns+spacer*2+line_width height-(position-1)*sheight]);
    
   
end

%unfortunately, it is not possible to force the box to be smaller than the
%original height, therefore, turn it off and set background color to none
%so that it no longer appears
set(legend_h, 'Color', 'None', 'Box', 'off');

%let's put it where you want it
pos = get(legend_h, 'position');
fig_pos = get(gca, 'position');
switch lower(location),
    case {'northeast'}
        set(legend_h, 'position', [pos(1)+fig_pos(3)-pos(3) pos(2) pos(3) pos(4)]);
    case {'southeast'}
        set(legend_h, 'position', [pos(1)+fig_pos(3)-pos(3) fig_pos(2)-pos(4)/2+pos(4)/4 pos(3) pos(4)]);
    case {'southwest'}
        set(legend_h, 'position', [fig_pos(1) fig_pos(2)-pos(4)/2+pos(4)/4 pos(3) pos(4)]);
end