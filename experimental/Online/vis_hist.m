function vis_hist(eeg,varargin)
% vis_hist(EEG,Options...)
% browse through windowed histograms and/or signals of an EEG dataset
%
% In:
%   EEG        : continuous-time dataset
%   Options... : name-value pairs specifying the options, with names:
%                --- data selection ---
%                'wnd': window length to compute histograms / draw signal sections, in seconds (default: 10)
%                'cfg': list of foreground channels, blue (default: 1)
%                'cbg': list of background channels, grey (default: 'all'); excluding foreground channels
%                --- histogram plotting ---
%                'hfg': show histogram of foreground channels (default: 1)
%                'hbg': show histogram of background channels (defalt: 0)
%                'bins': number of bins to use for histogram (default: 100)
%                'xlimit': X plot limits
%                'ylimit': Y plot limits
%                --- signal plotting ---
%                'fgcol': foreground signal color (default: [0 0 0.5])
%                'bgcol': background signal color (default: [0.5 0.5 0.5])
%                'sfg': show foreground signals (default: 1)
%                'sbg': show background signals (default: 1)
%                'yrange': relative range of the histogram's Y axis that shall be occupied by signals (default [0.1 0.8])
%                'yscaling': signal Y scaling (in std. devs. between signal rows) (default: 3)
%
% Examples:
%  vis_hist(myeeg,'xlimit',[-200 200],'ylimit',[0 200],'cfg',7,'hfg',1) %- show channel 7 in blue, rest in grey, plus a histogram of channel 7 with some fixed map limits
%  vis_hist(myeeg,'cfg',[1:30],'sfg',1,'sbg',0,'hfg',0,'hbg',0) %- show the first 10 channels in blue
%  vis_hist(myeeg,'cfg',[1:10],'sfg',1,'sbg',1,'hfg',0,'hbg',0) %- show the first 10 channels in blue and all others in grey
%  vis_hist(myeeg,'cfg',7,'sfg',1,'sbg',1,'hfg',1,'hbg',0) %- show channel 7 in blue, others in grey, and channel 7's histogram
%  vis_hist(myeeg,'cfg',7,'sfg',1,'sbg',1,'hfg',1,'hbg',1) %- show channel 7 in blue, others in grey, histograms of chan 7 and the grey channels
%  vis_hist(myeeg,'cfg',7,'cbg',5:15,'sfg',1,'sbg',1,'hfg',1,'hbg',1) %- show channel 7 in blue, 5:15 in grey, histograms of chan 7 and the grey channels
%  vis_hist(myeeg,'cfg',7,'cbg',5:15,'sfg',1,'sbg',1,'hfg',1,'hbg',1,'xlimit',[-200 200],'ylimit',[0. 500]) % - set map limits for the histograms
% 
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-07-03

global gHist; gHist = struct('hSlider',[],'lastPos',0);

if iscell(eeg)
    gHist.eeg2 = eeg{2};
    eeg = eeg{1};
end

% create figure & slider
figure('ResizeFcn',@wResize); hold;
gHist.hSlider = uicontrol('style','slider'); slideResize();
jSlider = findjobj(gHist.hSlider);
jSlider.AdjustmentValueChangedCallback = @sliderUpdate;

if ndims(eeg.data) == 3
    eeg.data = eeg.data(:,:);
    [eeg.nbchan,eeg.pnts,eeg.trials] = size(eeg.data);
end

% populate with data
gHist.eeg = eeg;
gHist.opts = hlp_varargin2struct(varargin,'cfg',[],'cbg','all', 'hfg',0,'hbg',0, ...
    'wnd',10,'bins',100,'xlimit',[],'ylimit',[],'sfg',1,'sbg',1,'yrange',[0.1 0.8],'yscaling',3,'fgcol',[0 0 0.5],'bgcol',[0.5 0.5 0.5]);

% do the initial update
sliderUpdate;



% visualize the current eeg's histogram at the given relative position
function visualize(relPos,moved)
global gHist;

% read parameters
eeg = gHist.eeg;
opts = gHist.opts;
for fn = fieldnames(opts)'
    eval([fn{1} ' = opts.' fn{1} ';']); end

% if this happens, we are maxing out MATLAB's graphics pipeline: let it catch up
if relPos == gHist.lastPos && moved
    return; end

% find out the background channel list
[C,T] = size(eeg.data);
if strcmp(cbg,'all')
    cbg = 1:C; end
if strcmp(cfg,'all')
    cfg = 1:C; end
if ~isequal(cfg,cbg)
    % do a setdiff if they're not equal
    cbg = setdiff_bc(cbg,cfg); 
end
    

% title
title(sprintf('[%.1f - %.1f]',eeg.xmax*relPos,eeg.xmax*relPos+wnd(end)));

% sampling window
wnd = 1:wnd*eeg.srate;
% sparse index set to sample from background channels, for the grey histogram (otherwise too slow)
if ~isfield(gHist,'indices')
    gHist.indices = 1+floor((length(wnd)*length(cbg)-1)*rand(1,length(wnd))); end
% position into the data set
pos = floor((T-wnd(end))*relPos);

% axes
cla;
if ~isempty(xlimit)
    set(gca,'XLim',xlimit); end
if ~isempty(ylimit)
    set(gca,'YLim',ylimit); end

% background histogram in gray
if hbg &&~isempty(cbg)
    tmp = eeg.data(cbg,pos + wnd); 
    hist(tmp(gHist.indices),bins);
    h = findobj(gca,'Type','patch'); set(h,'FaceColor',[0.5 0.5 0.5])
end

% foreground histograms in blue
if hfg && ~isempty(cfg)
    for c = cfg
        hist(eeg.data(c,pos + wnd),bins); end
end

% compute some signal plotting parameters (map into the histogram's axes...)
xl = get(gca,'XLim');
yl = get(gca,'YLim');
fp = get(gcf,'position'); 
ap = get(gca,'position');
ylr = yl(1) + yrange*(yl(2)-yl(1));
if ~isfield(gHist,'datascale')
    gHist.datascale = iqr(eeg.data,2); end
scale = ((ylr(2)-ylr(1))/size(eeg.data,1)) ./ (yscaling*gHist.datascale); scale(~isfinite(scale)) = 0;
%scale = ((ylr(2)-ylr(1))/size(eeg.data,1)) ./ (yscaling*std(eeg.data(:,pos+wnd),0,2)); scale(~isfinite(scale)) = 0;
channel_y = (ylr(2):(ylr(1)-ylr(2))/(size(eeg.data,1)-1):ylr(1))';
pixels = (fp(3))*(ap(3)-ap(1));
wnd = 1+floor(0:wnd(end)/(pixels):(wnd(end)-1));
scale = repmat(scale,1,length(wnd)); 

% draw background signals
if sbg && ~isempty(cbg)
    if isfield(gHist,'eeg2')
        plot(xl(1):(xl(2)-xl(1))/(length(wnd)-1):xl(2), (repmat(channel_y(cbg),1,length(wnd)) + scale(cbg,:).*gHist.eeg2.data(cbg,pos+wnd))','Color',opts.bgcol);
    else
        plot(xl(1):(xl(2)-xl(1))/(length(wnd)-1):xl(2), (repmat(channel_y(cbg),1,length(wnd)) + scale(cbg,:).*eeg.data(cbg,pos+wnd))','Color',opts.bgcol);
    end
end
% draw foreground signals
if sfg && ~isempty(cfg)
    plot(xl(1):(xl(2)-xl(1))/(length(wnd)-1):xl(2), (repmat(channel_y(cfg),1,length(wnd)) + scale(cfg,:).*eeg.data(cfg,pos+wnd))','Color',opts.fgcol); end
drawnow;

gHist.lastPos = relPos;




% slider moved
function sliderUpdate(varargin)
global gHist;
visualize(get(gHist.hSlider,'Value'),~isempty(varargin))

% adapt/set the slider's size
function slideResize(varargin)
global gHist;
wPos = get(gcf,'Position'); 
if ~isempty(gHist.hSlider)
    try
        set(gHist.hSlider,'Position',[20,20,wPos(3)-40,20]); 
    catch,end    
    if isfield(gHist,'eeg')
        sliderUpdate; end
end

% window resized
function wResize(varargin)
slideResize();