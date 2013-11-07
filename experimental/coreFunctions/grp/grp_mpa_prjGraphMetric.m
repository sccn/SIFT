function varargout = grp_mpa_prjGraphMetric(varargin)
% [STUDY obj] = grp_mpa_prjGraphMetric(STUDY,ALLEEG,...)
%
% perform causal projection analysis using the Measure Projection Toolbox
% Output: 
%       STUDY:  modified STUDY with MPT dipoleAndMeasure object in STUDY.measureprojection.sift.object
%       obj:    (optional) dipoleAndMeasure object

% extract some stuff from inputs for arg defaults
ALLEEG = arg_extract(varargin,'ALLEEG',2);

if ~isempty(ALLEEG) && isempty(hlp_checkeegset(ALLEEG(1),{'conn'}));
    Conn            = ALLEEG(1).CAT.Conn;
    ConnNames       = hlp_getConnMethodNames(Conn);
    freqRangeDef    = [Conn.freqs(1) Conn.freqs(end)];
    timeRangeDef    = [Conn.erWinCenterTimes(1) Conn.erWinCenterTimes(end)];
    clear Conn;
else
    ConnNames = {''};
    [freqRangeDef, timeRangeDef] = deal([]);
end

g = arg_define(varargin, ...
        arg('STUDY',mandatory,[],'STUDY object'),  ...
        arg('ALLEEG',mandatory,[],'ALLEEG array'), ...
        arg({'connmethod','Estimator'},ConnNames,ConnNames,'Connectivity estimator to project'), ...
        arg({'timeRange','TimeRange'},timeRangeDef,[],'[Min Max] Time range to smooth (sec). Leave blank to use all time points','shape','row','type','denserealdouble'), ...
        arg({'freqRange','FrequencyRange'},freqRangeDef,[],'[Min Max] Frequency range to smooth (Hz). Leave blank to use all frequencies','type','expression','shape','row'), ...
        arg({'collapseTime','CollapseTime'},false,[],'Average across time before smoothing'), ...
        arg({'collapseFreqs','CollapseFreqs'},true,[],'Integrate across frequencies before smoothing'), ...
        arg({'rejOutsideBrain','RejectDipolesOutsideBrain'},true,[],'Reject dipoles outside the brain'), ...
        arg_sub({'graphMetric','GraphMetric'},@hlp_computeGraphMeasure,{},'Graph metric options','suppress',{'srcNodes','targNodes'}), ...
        arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );
    
% basic error checking
if g.collapseTime && g.collapseFreqs
    error('At most one dimension can be collapsed. Otherwise, there is nothing to smooth!');
end
if any(cellfun(@(x) isempty(hlp_checkeegset(x,{'conn'})), ALLEEG))
    error('You must first run SIFT analysis on all datasets in the STUDY');
end

waitbarstr = sprintf('Creating Measure Projection Object for %s',g.graphMetric.graphMetric);
if g.verb==2
    % reset the color list
    hlp_getNextUniqueColor(hsv(10),'reset');
    multiWaitbar(waitbarstr,'Reset','Color',hlp_getNextUniqueColor);
elseif g.verb==1
    fprintf([waitbarstr '...']);
end

% make sure that you have loaded STUDY in which all datasets have SIFT measures
dipoleAndMeasure = pr.dipoleAndMeasureOfStudySIFT(g.STUDY, g.ALLEEG);
% reject dipoles outside the brain
if g.rejOutsideBrain
    dipoleAndMeasure = dipoleAndMeasure.createSubsetInRelationToBrain();
end
    
for sidx = 1:length(g.ALLEEG)
    % get indices of dipoles for this dataset
    dipidx = dipoleAndMeasure.datasetIdAllConditions==sidx;
    % determine source and target node indices
    g.graphMetric.srcNodes  = dipoleAndMeasure.numberInDataset(dipidx);
    g.graphMetric.targNodes = g.graphMetric.srcNodes;
    % compute graph measure [nodes x freq x time]
    [gm nt nf] = ComputeGraphMeasure(g.ALLEEG(sidx).CAT.Conn, g);
    % subtract mean
    gm  = bsxfun(@minus, gm, mean(gm,1));
    % store vectorized measures
    dipoleAndMeasure.linearizedMeasure(:,dipidx) = gm(:,:)';
    % update progress
    if g.verb==2
        multiWaitbar(waitbarstr,sidx/length(g.ALLEEG));
    elseif g.verb==1, fprintf('.'); end
end

% update some fields
dipoleAndMeasure.time = nt;
dipoleAndMeasure.frequency = nf;
dipoleAndMeasure.measureLabel = g.graphMetric.graphMetric;
dipoleAndMeasure.numberOfMeasureDimensions = ndims(gm)-1;

% store results and prepare outputs
STUDY.measureprojection.sift.object =  dipoleAndMeasure;
varargout{1} = STUDY;
if nargout>1, varargout{2} = dipoleAndMeasure; end

if g.verb==2
    multiWaitbar(waitbarstr,'Close');
elseif g.verb==1, fprintf('\n'); end


% Helper functions
% --------------------------------------------------------------
function [gm newtimes newfreqs] = ComputeGraphMeasure(Conn,g)
% compute graph measures from a Connectivity object


% collapse connectivity matrices
if g.collapseFreqs
    % collapse across freqs
    Conn = hlp_filterConns(Conn,'connmethod',g.connmethod, ...
        'method',{'freq','net'},'frange',g.freqRange,'freqdim',3,'timedim',4,'verb',g.verb);
else
    % only select subset of freqs
    Conn = hlp_filterConns(Conn,'connmethod',g.connmethod, ...
        'method',{'freq','shrinkonly'}, ...
        'frange',g.freqRange, 'freqdim',3,'timedim',4,'verb',g.verb);
end
if g.collapseTime
    % collapse across time
    Conn = hlp_filterConns(Conn,'connmethod',g.connmethod, ...
        'method',{'time','mean'},'trange',g.timeRange,'freqdim',3,'timedim',4,'verb',g.verb);
else
    % only select subset of times
    Conn = hlp_filterConns(Conn,'connmethod',g.connmethod, ...
        'method',{'time','shrinkonly'}, ...
        'trange',g.timeRange,'freqdim',3,'timedim',4,'verb',g.verb);
end
newtimes = Conn.erWinCenterTimes;
newfreqs = Conn.freqs;
% compute graph measure
gm = hlp_computeGraphMeasure('cmatrix',Conn.(g.connmethod),g.graphMetric);









    
    
%  arg_sub({'mpaOpts','MPA_Options'}, ...
%             { ...
%             arg({'rejOutsideBrain','RejectDipolesOutsideBrain'},true,[],'Reject dipoles outside the brain'), ...
%             arg({'similarityThreshold','SimilarityThreshold'},0.9,[0 1],'Similarity Threshold'), ...
%             
%             similarityThreshold, significanceLevel, minNumberOfClusters, outlierThreshold
% 




% 
% % compute optimal Gaussian width
% [optimalGaussianWidth meanPredictionSimilarity gaussianStdValues] = pr.find_optimal_gaussian_width(dipoleAndMeasure);
% 
% % compute group-difference
% gmProj = pr.meanProjection(dipoleAndMeasure, dipoleAndMeasure.getPairwiseFishersZSimilarity, pr.headGrid, 'stdOfDipoleGaussian', optimalGaussianWidth);
% 
% % project measures
% gmProj = gmProj.createDomain(dipoleAndMeasure, g.mpaOpts, 0.05);
%     


% % build 'linearizedMeasure' under dipoleAndMeasure 
% for i=1:length(dipoleAndMeasure.datasetIdAllConditions)
%     outflow = ComputeGraphMeasure(ALLEEG(dipoleAndMeasure.datasetIdAllConditions(i)).CAT.Conn);
%     
% %     outflow = ALLEEG(dipoleAndMeasure.datasetIdAllConditions(i)).outflow(:,dipoleAndMeasure.numberInDataset(i));
%     % outflow = outflow(goodFreqIdx);
%     % outflow = log(outflow);
%     % outflow = outflow - groupMeanOutflow'; % subtract mean
%     outflow = bsxfun(@minus, outflow, mean(outflow)); % subtract mean
%     dipoleAndMeasure.linearizedMeasure(:,i) = outflow;
% end

