function [bg bgh CollapsedConn GraphMeasure] = MakeConnectivityGraph(Conn,connmethod,TimeRangeToCollapse,FreqRangeToCollapse,UnivariateGraphMeasure,PermutedChanOrder,FreqCollapseMode,TimeCollapseMode,PlotFigures,NodeLabels,ChansToExclude,EdgeSizeLimits,NodeSizeLimits,LayoutType,MakeLegend)

% TimeRangeToCollapse can be a matrix of [K x 2] intervals indicating the
% [min max] time range for the kth segement
% UnivariateGraphMeasure can be one of the graph measures from
% hlp_computeGraphMeasure (e.g. 'causalflow' or 'outflow'). Leave empty for
% none
%
% LayoutType can be 'radial','hierarchical', or a [NumChans x Num_Coords]
% matrix where the ith row contains the [X Y] or [X Y Z] coordinates of the
% ith channel/source


if nargin<4
    error('you must supply 4 arguments');
end

if nargin<5
    UnivariateGraphMeasure = '';
end

if nargin<6
    PermutedChanOrder = [];
end

if nargin<7
    FreqCollapseMode = 'net';
end

if nargin<8
    TimeCollapseMode = 'mean';
end

if nargin<9
    PlotFigures = true;
end

if nargin<10
    NodeLabels = [];
end

if nargin<11
    ChansToExclude = [];
end

if nargin<12
    EdgeSizeLimits = 99;  % percentile range for EdgeSize
end

if nargin<13
    NodeSizeLimits = 99;  % percentile range for NodeSize
end

if nargin<14
    LayoutType = 'radial'; % hierarchical
end

if nargin<15
    MakeLegend = true;
end


%% defaults

bgh = {};
bg  = {};

mycolormap = [0.941176474094391,0.941176474094391,0.941176474094391;0.931764721870422,0.931764721870422,0.931764721870422;0.922352969646454,0.922352969646454,0.922352969646454;0.912941157817841,0.912941157817841,0.912941157817841;0.903529405593872,0.903529405593872,0.903529405593872;0.894117653369904,0.894117653369904,0.894117653369904;0.884705901145935,0.884705901145935,0.884705901145935;0.875294148921967,0.875294148921967,0.875294148921967;0.865882337093353,0.865882337093353,0.865882337093353;0.856470584869385,0.856470584869385,0.856470584869385;0.847058832645416,0.847058832645416,0.847058832645416;0.837647080421448,0.837647080421448,0.837647080421448;0.828235328197479,0.828235328197479,0.828235328197479;0.818823516368866,0.818823516368866,0.818823516368866;0.809411764144898,0.809411764144898,0.809411764144898;0.800000011920929,0.800000011920929,0.800000011920929;0.775000035762787,0.806250035762787,0.806250035762787;0.750000000000000,0.812500000000000,0.812500000000000;0.725000023841858,0.818750023841858,0.818750023841858;0.699999988079071,0.824999988079071,0.824999988079071;0.675000011920929,0.831250011920929,0.831250011920929;0.650000035762787,0.837500035762787,0.837500035762787;0.625000000000000,0.843750000000000,0.843750000000000;0.600000023841858,0.850000023841858,0.850000023841858;0.574999988079071,0.856249988079071,0.856249988079071;0.550000011920929,0.862500011920929,0.862500011920929;0.525000035762787,0.868750035762787,0.868750035762787;0.500000000000000,0.875000000000000,0.875000000000000;0.474999994039536,0.881250023841858,0.881250023841858;0.450000017881393,0.887499988079071,0.887499988079071;0.425000011920929,0.893750011920929,0.893750011920929;0.400000005960465,0.899999976158142,0.899999976158142;0.375000000000000,0.906250000000000,0.906250000000000;0.349999994039536,0.912500023841858,0.912500023841858;0.325000017881393,0.918749988079071,0.918749988079071;0.300000011920929,0.925000011920929,0.925000011920929;0.275000005960465,0.931249976158142,0.931249976158142;0.250000000000000,0.937500000000000,0.937500000000000;0.225000008940697,0.943750023841858,0.943750023841858;0.200000002980232,0.949999988079071,0.949999988079071;0.174999997019768,0.956250011920929,0.956250011920929;0.150000005960464,0.962499976158142,0.962499976158142;0.125000000000000,0.968750000000000,0.968750000000000;0.100000001490116,0.975000023841858,0.975000023841858;0.0750000029802322,0.981249988079071,0.981249988079071;0.0500000007450581,0.987500011920929,0.987500011920929;0.0250000003725290,0.993749976158142,0.993749976158142;0,1,1;0.0312500000000000,1,0.968750000000000;0.0625000000000000,1,0.937500000000000;0.0937500000000000,1,0.906250000000000;0.125000000000000,1,0.875000000000000;0.156250000000000,1,0.843750000000000;0.187500000000000,1,0.812500000000000;0.218750000000000,1,0.781250000000000;0.250000000000000,1,0.750000000000000;0.281250000000000,1,0.718750000000000;0.312500000000000,1,0.687500000000000;0.343750000000000,1,0.656250000000000;0.375000000000000,1,0.625000000000000;0.406250000000000,1,0.593750000000000;0.437500000000000,1,0.562500000000000;0.468750000000000,1,0.531250000000000;0.500000000000000,1,0.500000000000000;0.531250000000000,1,0.468750000000000;0.562500000000000,1,0.437500000000000;0.593750000000000,1,0.406250000000000;0.625000000000000,1,0.375000000000000;0.656250000000000,1,0.343750000000000;0.687500000000000,1,0.312500000000000;0.718750000000000,1,0.281250000000000;0.750000000000000,1,0.250000000000000;0.781250000000000,1,0.218750000000000;0.812500000000000,1,0.187500000000000;0.843750000000000,1,0.156250000000000;0.875000000000000,1,0.125000000000000;0.906250000000000,1,0.0937500000000000;0.937500000000000,1,0.0625000000000000;0.968750000000000,1,0.0312500000000000;1,1,0;1,0.968750000000000,0;1,0.937500000000000,0;1,0.906250000000000,0;1,0.875000000000000,0;1,0.843750000000000,0;1,0.812500000000000,0;1,0.781250000000000,0;1,0.750000000000000,0;1,0.718750000000000,0;1,0.687500000000000,0;1,0.656250000000000,0;1,0.625000000000000,0;1,0.593750000000000,0;1,0.562500000000000,0;1,0.531250000000000,0;1,0.500000000000000,0;1,0.468750000000000,0;1,0.437500000000000,0;1,0.406250000000000,0;1,0.375000000000000,0;1,0.343750000000000,0;1,0.312500000000000,0;1,0.281250000000000,0;1,0.250000000000000,0;1,0.218750000000000,0;1,0.187500000000000,0;1,0.156250000000000,0;1,0.125000000000000,0;1,0.0937500000000000,0;1,0.0625000000000000,0;1,0.0312500000000000,0;1,0,0;0.968750000000000,0,0;0.937500000000000,0,0;0.906250000000000,0,0;0.875000000000000,0,0;0.843750000000000,0,0;0.812500000000000,0,0;0.781250000000000,0,0;0.750000000000000,0,0;0.718750000000000,0,0;0.687500000000000,0,0;0.656250000000000,0,0;0.625000000000000,0,0;0.593750000000000,0,0;0.562500000000000,0,0;0.531250000000000,0,0;0.500000000000000,0,0;];

DEF_EDGESIZEMAX = 0.3;
DEF_NODESIZE    = 0.2;

% linewidth preferences
lMAX = 10;
lMIN = 0;

ComputeGraphMeasure = ~isempty(UnivariateGraphMeasure);
NumTimeWindows = size(TimeRangeToCollapse,1);
nchs = size(Conn.(connmethod),1)-length(ChansToExclude);

if isempty(PermutedChanOrder)
    PermutedChanOrder = 1:nchs;
end

if isempty(FreqRangeToCollapse) && length(Conn.freqs)>1
    FreqRangeToCollapse = [Conn.freqs(1) Conn.freqs(end)];
end

if isempty(TimeRangeToCollapse) && length(Conn.erWinCenterTimes)>1
    TimeRangeToCollapse = [Conn.erWinCenterTimes(1) Conn.erWinCenterTimes(end)];
end

%%


% obtain spectrum, if desired
if strcmpi(UnivariateGraphMeasure,'spectrum')
    if ~isfield(Conn,'S')
        fprintf('Warning: Spectrum is not present in Connectivity structure. Unable to use spectral information.\n');
        UnivariateGraphMeasure = '';
        ComputeGraphMeasure = false;
    else
        % extract only the 'S' connectivity matrix from Conn
        ConnSpect = rmfield(Conn,setdiff(hlp_getConnMethodNames(Conn),'S'));
    end
end

% remove all connectivity matrices that are not of interest
Conn = rmfield(Conn,setdiff(hlp_getConnMethodNames(Conn),connmethod));

CollapsedConn = Conn;
CollapsedConn.(connmethod) = zeros(nchs,nchs,NumTimeWindows);

% obtain labels for the nodes (applying permutation if desired)
GoodChanNumbers = setdiff(1:nchs,ChansToExclude);
if isempty(NodeLabels)
    NodeLabels = strtrim(cellstr(num2str(GoodChanNumbers(PermutedChanOrder)')))';
else
    NodeLabels = NodeLabels(GoodChanNumbers);
    NodeLabels = NodeLabels(PermutedChanOrder);
end


% initialize matrices
if ComputeGraphMeasure,
    GraphMeasure = zeros(nchs,NumTimeWindows);
else
    GraphMeasure = [];
end

for t=1:NumTimeWindows
    % for each time window, collapse the matrix
    
    tmp = hlp_filterConns(Conn,'connmethods',{connmethod},'method',{'freq',FreqCollapseMode,'time',TimeCollapseMode},'trange',TimeRangeToCollapse(t,:),'frange',FreqRangeToCollapse,'freqdim',3,'timedim',4,'badchans',ChansToExclude,'verb',false);
    
    % set the diagonals to zero (no auto-connections allowed)
    CollapsedConn.(connmethod)(:,:,t) = tmp.(connmethod).*~eye(nchs);
    
    % apply permutation of channel order if desired
    CollapsedConn.(connmethod)(:,:,t) = CollapsedConn.(connmethod)(PermutedChanOrder,PermutedChanOrder,t);
    
    if strcmpi(UnivariateGraphMeasure,'spectrum')
        % compute the spectral power for each channel
        % (e.g. for NodeSize modulation)
        
        tmp = hlp_filterConns(ConnSpect,'connmethods',{'S'},'method',{'freq',FreqCollapseMode,'time',TimeCollapseMode},'trange',TimeRangeToCollapse(t,:),'frange',FreqRangeToCollapse,'freqdim',3,'timedim',4,'badchans',ChansToExclude,'verb',false);
        GraphMeasure(:,t) = diag(tmp.('S'));
        
        % apply permutation of channel order if desired
        GraphMeasure(:,t) = GraphMeasure(PermutedChanOrder,t);
        
    elseif ComputeGraphMeasure
        % compute a univariate graph measure for each channel
        
        Ctmp = squeeze(CollapsedConn.(connmethod)(:,:,t));
        
        for ch=1:size(Ctmp,2)
            GraphMeasure(ch,t) = hlp_computeGraphMeasure(Ctmp,ch,1:nchs,lower(UnivariateGraphMeasure));
        end
    end
    
end

clear('tmp');

%% Compute percentiles for Node and Edge size limits
EdgeSizeMax = prctile(CollapsedConn.(connmethod)(:),EdgeSizeLimits);
EdgeSizeMin = prctile(CollapsedConn.(connmethod)(:),100-EdgeSizeLimits);

if ComputeGraphMeasure
    NodeSizeMax = prctile(GraphMeasure(:),NodeSizeLimits);
    NodeSizeMin = prctile(GraphMeasure(:),100-NodeSizeLimits);
end


%% Set up the NodeSize

if ComputeGraphMeasure
    % scale NodeSize according to graph measure
    
    % max and min radii (as fraction of radius of circlegraph)
    rMAX = 0.13;
    rMIN = 0.05;
    
    rmin = NodeSizeMin; %-EdgeSizeMax;
    rmax = NodeSizeMax; % EdgeSizeMax;
    
    % rescale to [0 1]
    NodeSize = rMIN + (rMAX-rMIN).*mat2gray(GraphMeasure,double([rmin,rmax]));
    
else
    NodeSize = DEF_NODESIZE;
end

% Plot graphs for each of the time windows

for t=1:NumTimeWindows
    
    % get the [nchs_to x nchs_from] connectivity matrix for time window t
    C = CollapsedConn.(connmethod)(:,:,t);
    
    if issym(C)
        csym = true;
        C = triu(C);
    else
        csym = false;
    end
    
    if all(C==0)
        fprintf('No non-zero connections at time [%s]. Graph will not be constructed\n', ...
                regexprep(num2str(TimeRangeToCollapse(t,:)),'\s*','  '));
        bg{t}   = [];
        bgh{t}  = [];
        continue;
    else
    
        % create the biograph
        bg{t} = biograph(C',NodeLabels,'ShowWeights','off',...
                                       'ShowTextInNodes','Label', ...
                                       'EdgeType','curved', ...
                                       'NodeAutoSize','off', ...
                                       'LayoutType',fastif(ischar(LayoutType),LayoutType,'radial'));
    end
    
    if PlotFigures
        % plot the graph
        if ismatrix(LayoutType) || strcmpi(LayoutType,'radial')
            bh = bgCircleGraph(bg{t},NodeSize); % view(bg{t});
        else
            bh = view(bg{t});
        end
        
        if csym % symmetric connectivity
            set(bh,'showArrows','off');
%             dolayout(bh,'pathsOnly',true);
        end
        
        % bh=view(biograph(C',ComponentNames,'ShowWeights','on','LayoutType','radial','ShowTextInNodes','Label','Scale',1.5)); %'hierarchical', 'radial'
        
        % if colorclusters
        %     nodecolors = distinguishable_colors(numclusters,[1 1 1; 0 0 0; 0 0 1]);
        %     % make each cluster a different color
        %     for i=1:numclusters
        %         set(bh.Nodes(myclustering(myclustering(:,2)==i,1)),'Color',nodecolors(i,:));
        %     end
        % end
        
        for i=1:length(bh.Nodes)
            set(bh.Nodes(i),'label',NodeLabels{i},'FontSize',12);
        end
        
        % set the node locations
        if ismatrix(LayoutType)
            for i=1:length(bh.Nodes)
                set(bh.Nodes(i),'Position',LayoutType(i,:));
            end
            
            % redraw the edges of the graph
            bh.dolayout('pathsonly',true);
            bh.hgUpdate
            bh.hgReDraw
        end
        
        cmaplen = 128;           %length(find(C_cluster>0))
        edgecolors = mycolormap; %jet(cmaplen);
        
        cmin = EdgeSizeMin;     %0;  min(C(C>0));
        cmax = EdgeSizeMax;     %    max(C(:));
        lmax = EdgeSizeMax;     %    max(C(:));
        lmin = EdgeSizeMin;     %0;  min(C(C>0));
        
        for c1=1:size(C,1)
            for c2=1:size(C,1)
                cval=C(c1,c2);
                
                if cval>0
                    % determine line color
                    index = min(cmaplen,fix((cval-cmin)/(cmax-cmin)*(cmaplen-1))+1);
                    LineColor = edgecolors(index,:);
                    
                    % determine linewidth
                    
                    LineWidth = min(lMAX,1+(cval-lmin)*(lMAX-lMIN)/(lmax-lmin+eps));
                    c2name = NodeLabels{c2};
                    c1name = NodeLabels{c1};
                    set(findobj(bh.Edges,'ID',sprintf('%s -> %s',c2name,c1name)),'LineColor',LineColor,'LineWidth',LineWidth);
                end
                
            end
        end
        
        % Add a few more things to the graph
        fighandle = get(bh.hgAxes,'parent');
        
        if MakeLegend
            
            % create legend
            
            timestr = regexprep(num2str(TimeRangeToCollapse(t,:)),'\s*','  ');
            freqstr = regexprep(num2str(FreqRangeToCollapse),'\s*','  ');
            
            edgeleg     = sprintf('EdgeSize:\t%s  ',connmethod);
            nodeleg     = fastif(isempty(UnivariateGraphMeasure),{},{sprintf('NodeSize:\t%s  ',UnivariateGraphMeasure)});
            timesleg    = fastif(isnan(TimeRangeToCollapse),{}, ...
                                 sprintf('Times:       [%s] sec ',fastif(isempty(timestr),Conn.erWinCenterTimes,num2str(timestr))));
            freqleg     = fastif(isnan(FreqRangeToCollapse),{}, ...   
                                sprintf('Freqs:       [%s] Hz ',  fastif(isempty(freqstr),num2str(Conn.freqs),freqstr)));
            legend = [...
                      edgeleg, ...
                      nodeleg, ...
                      timesleg, ...
                      freqleg ...
                    ];
            
            % TODO: need to set up callbacks to allow the legend to be moved
            % NOTE: if we want the legend to be preserved in 'print to figure',
            %       then we must use the 'text' version (no FaceAlpha available)
            %         ha = annotation(fighandle,'textbox',[0.6895, 0.8784,0.1,0.2],'FaceAlpha',0.5,'FitBoxToText','on','String',legend,'units','normalized','EdgeColor','k','BackGroundColor',[1 0.7 0.7],'FontSize',10);
            ht=text(0.6895, 0.8784,legend,'parent',bh.hgAxes,'units','normalized','EdgeColor','k','BackGroundColor',[1 0.7 0.7],'FontSize',10); % 'BackGroundColor',[1 0.7 0.7]
            %         pos = get(ht,'Extent');
            %         axespos = get(bh.hgAxes,'Position');
            % transform pos to data units
            %         pos = pos.*repmat(axespos(3:4),1,2);
            %         hp = patch([pos(1) pos(1) pos(1)+pos(3) pos(1)+pos(3)],[pos(2) pos(2)+pos(4) pos(2)+pos(4) pos(2)],[1 0.7 0.7],'FaceAlpha',0.3,'parent',bh.hgAxes);
        end
        
        set(fighandle,'Name',sprintf('Connectivity Graph %d',t));
        
        %         set(fighandle,'Name',sprintf('Times: [%0.3g %0.3g] sec, Freqs: [%0.4g %0.4g] Hz',...
        %             TimeRangeToCollapse(t,1),TimeRangeToCollapse(t,2),FreqRangeToCollapse(1),FreqRangeToCollapse(2)));
        if NumTimeWindows>1
            set(fighandle,'WindowStyle','docked');
        end
        
        bgh{t} = bh;
        
    end

end



end
