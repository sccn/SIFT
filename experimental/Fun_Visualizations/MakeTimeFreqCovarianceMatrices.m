% compute the time- and frequency-based covariance matrices for each channel

EEG.CAT.Conn = hlp_absvalsq(EEG.CAT.Conn,hlp_getConnMethodNames(EEG.CAT.Conn),true,false,true);

ConnCov = EEG.CAT.Conn;
connmethod = 'RPDC';
method = 'cov';
Offset = 0; 0.980;  % time shift for event onset

FreqRangeToCollapse = [110 200]; [ConnCov.freqs(1) ConnCov.freqs(end)];
TimeRangeToCollapse = [0 0.5]; %[0 0.5];


ConnCov = hlp_filterConns(ConnCov,'connmethods',{connmethod},'method',{'freq','shrinkonly','time','shrinkonly'},'trange',TimeRangeToCollapse,'frange',FreqRangeToCollapse,'freqdim',3,'timedim',4,'verb',false);


C = ConnCov.(connmethod);

% remove all connectivity methods
ConnCov = rmfield(ConnCov,hlp_getConnMethodNames(ConnCov));

[nchs nchs nfreqs ntimes] = size(C);

ConnCov.(sprintf('%s_tcov',connmethod)) = zeros(nchs,nchs,ntimes,ntimes);
ConnCov.(sprintf('%s_fcov',connmethod)) = zeros(nchs,nchs,nfreqs,nfreqs);

StatCov.(sprintf('%s_tcov',connmethod)).pval = zeros(nchs,nchs,ntimes,ntimes);
StatCov.(sprintf('%s_fcov',connmethod)).pval = zeros(nchs,nchs,nfreqs,nfreqs);
StatCov.(sprintf('%s_tcov',connmethod)).ci = zeros(2,nchs,nchs,ntimes,ntimes);
StatCov.(sprintf('%s_fcov',connmethod)).ci = zeros(2,nchs,nchs,nfreqs,nfreqs);
StatCov.alpha = 0.00001;

% take the sqrt before transformation. This ensures that a pixel on the diagonal of the
% time-covariance matrix represents the total connectivity, integrated
% across freq. The off-diagonals represent covariance between integrated
% connectivities
C = sqrt(C);
for ch1=1:size(C,1)
    for ch2=1:size(C,2)
        
        switch method
            case 'corr'
                [Cor pval cilo cihi] = corrcoef(squeeze(C(ch1,ch2,:,:))','alpha',StatCov.alpha); % squeeze(C(ch1,ch2,:,:))*squeeze(C(ch1,ch2,:,:))';
                ConnCov.(sprintf('%s_fcov',connmethod))(ch1,ch2,:,:) = Cor;
                StatCov.(sprintf('%s_fcov',connmethod)).pval(ch1,ch2,:,:) = pval;
                StatCov.(sprintf('%s_fcov',connmethod)).ci(1,ch1,ch2,:,:) = cilo;
                StatCov.(sprintf('%s_fcov',connmethod)).ci(2,ch1,ch2,:,:) = cihi;
                [Cor pval cilo cihi] = corrcoef(squeeze(C(ch1,ch2,:,:)),'alpha',StatCov.alpha); %squeeze(C(ch1,ch2,:,:))'*squeeze(C(ch1,ch2,:,:));
                ConnCov.(sprintf('%s_tcov',connmethod))(ch1,ch2,:,:) = Cor;
                StatCov.(sprintf('%s_tcov',connmethod)).pval(ch1,ch2,:,:) = pval;
                StatCov.(sprintf('%s_tcov',connmethod)).ci(1,ch1,ch2,:,:) = cilo;
                StatCov.(sprintf('%s_tcov',connmethod)).ci(2,ch1,ch2,:,:) = cihi;
            case 'cov'
                Cor = cov(squeeze(C(ch1,ch2,:,:))'); % squeeze(C(ch1,ch2,:,:))*squeeze(C(ch1,ch2,:,:))';
                [dummy pval cilo cihi] = corrcoef(squeeze(C(ch1,ch2,:,:))','alpha',StatCov.alpha); % squeeze(C(ch1,ch2,:,:))*squeeze(C(ch1,ch2,:,:))';
                pval = pval.*~eye(size(pval));   % diagonals are always significant and we want to show variance
                ConnCov.(sprintf('%s_fcov',connmethod))(ch1,ch2,:,:) = Cor;
                StatCov.(sprintf('%s_fcov',connmethod)).pval(ch1,ch2,:,:) = pval;
%                 StatCov.(sprintf('%s_fcov',connmethod)).ci(1,ch1,ch2,:,:) = 0;
%                 StatCov.(sprintf('%s_fcov',connmethod)).ci(2,ch1,ch2,:,:) = 0;
                Cor = cov(squeeze(C(ch1,ch2,:,:))); %squeeze(C(ch1,ch2,:,:))'*squeeze(C(ch1,ch2,:,:));
                [dummy pval cilo cihi] = corrcoef(squeeze(C(ch1,ch2,:,:)),'alpha',StatCov.alpha); %squeeze(C(ch1,ch2,:,:))'*squeeze(C(ch1,ch2,:,:));
                pval = pval.*~eye(size(pval));   % diagonals are always significant
                ConnCov.(sprintf('%s_tcov',connmethod))(ch1,ch2,:,:) = Cor;
                StatCov.(sprintf('%s_tcov',connmethod)).pval(ch1,ch2,:,:) = pval;
%                 StatCov.(sprintf('%s_tcov',connmethod)).ci(1,ch1,ch2,:,:) = 0;
%                 StatCov.(sprintf('%s_tcov',connmethod)).ci(2,ch1,ch2,:,:) = 0;
            otherwise
                error('unknown method %s',method);
        end
        
    end
end


freqs = ConnCov.freqs;
times = ConnCov.erWinCenterTimes;

%%

Clim = 99.5;

figureHandles = [];

% plot time x time covariance matrices
ConnCov.erWinCenterTimes = times;
ConnCov.freqs = times;
ConnCov = hlp_filterConns(ConnCov,'connmethods',{sprintf('%s_tcov',connmethod)},'sigThresh',0);
[figureHandles(end+1) tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG,'Conn',hlp_absvalsq(ConnCov,{sprintf('%s_tcov',connmethod)},false,false,true),'Stats',StatCov,'EventMarkers',{{0-Offset,'r',':',2} {tM/1000-Offset, 'b', '-', 2}},'NodeLabels',anatnames,'SourceMarginPlot','none','MatrixLayout',{'Partial','UpperTriangle', sprintf('%s_tcov',connmethod), 'LowerTriangle',sprintf('%s_tcov',connmethod),'Diagonal','none'},'ColorLimits',Clim,'FrequencyScale','linear','Smooth2D',false,'Thresholding',{'Statistics','ThresholdingMethod','pval','PlotConfidenceIntervals',true,'AlphaSignificance',StatCov.alpha},'FrequencyMarkers',0,'FrequencyMarkerColor','r');


% plot freq x freq covariance matrices
ConnCov.freqs = freqs;
ConnCov.erWinCenterTimes = freqs;
[figureHandles(end+1) tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG,'Conn',hlp_absvalsq(ConnCov,{sprintf('%s_fcov',connmethod)},false,false,true),'Stats',StatCov,'NodeLabels',anatnames,'SourceMarginPlot','none','MatrixLayout',{'Partial','UpperTriangle', sprintf('%s_fcov',connmethod), 'LowerTriangle',sprintf('%s_fcov',connmethod),'Diagonal','none'},'ColorLimits',Clim,'FrequencyScale','linear','Smooth2D',false,'Thresholding',{'Statistics','ThresholdingMethod','pval','PlotConfidenceIntervals',true,'AlphaSignificance',StatCov.alpha});


%% also plot the freq-integrated and time-integrated spectra
[figureHandles(end+1) tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG,'Conn',EEG.CAT.Conn,'TimesToPlot',TimeRangeToCollapse,'FrequenciesToPlot',FreqRangeToCollapse(1):FreqRangeToCollapse(end),'EventMarkers',{{0-Offset,'r',':',2} {tM/1000-Offset, 'b', '-', 2}}, 'NodeLabels',anatnames,'SourceMarginPlot','none','MatrixLayout',{'Partial','UpperTriangle', connmethod, 'LowerTriangle',connmethod,'Diagonal','none'},'FrequencyScale','linear','Smooth2D',false);


% plot time x causality (integrate over freqs)
CollapsedConn = hlp_filterConns(EEG.CAT.Conn,'connmethods',{connmethod},'method',{'freq','net'},'frange',FreqRangeToCollapse,'freqdim',3,'timedim',4);
[figureHandles(end+1) tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG,'Conn',CollapsedConn,'TimesToPlot',TimeRangeToCollapse,'EventMarkers',{{0-Offset,'g',':',2} {tM/1000-Offset, 'b', '-', 2}}, 'NodeLabels',anatnames,'SourceMarginPlot','none','MatrixLayout',{'Partial','UpperTriangle', connmethod, 'LowerTriangle',connmethod,'Diagonal','none'},'FrequencyScale','linear','Smooth2D',false);

% plot freq x causality (integrate over time)
CollapsedConn = hlp_filterConns(EEG.CAT.Conn,'connmethods',{connmethod},'method',{'time','net'},'trange',TimeRangeToCollapse,'freqdim',3,'timedim',4);
[figureHandles(end+1) tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG,'Conn',CollapsedConn,'FrequenciesToPlot',FreqRangeToCollapse(1):FreqRangeToCollapse(2), 'NodeLabels',anatnames,'SourceMarginPlot','none','MatrixLayout',{'Partial','UpperTriangle', connmethod, 'LowerTriangle',connmethod,'Diagonal','none'},'FrequencyScale','linear','Smooth2D',false);


%% Dock the figures
 for i=figureHandles, set(i,'WindowStyle','docked'); end