function varargout = grp_mpa_prjGraphMetric(varargin)
% [STUDY obj] = grp_mpa_prjGraphMetric(STUDY,ALLEEG,...)
%
% Perform causal projection analysis [1] using the Measure Projection Toolbox [2]
%
% Output: 
%       STUDY:  modified STUDY with MPT dipoleAndMeasure object in STUDY.measureprojection.sift.object
%       obj:    (optional) dipoleAndMeasure object
% 
% make sure that you have loaded STUDY in which all datasets have SIFT measures
%
% References: 
% 
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.9
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
% [2] Bigdely-Shamlo, N, Mullen, T, Kreutz-Delgado,  K, Makeig, S, 
%   Measure Projection Analysis: A Probabilistic  Approach  to EEG Source Comparison and Multi-Subject Inference,
%   NeuroImage (2013), Vol. 72, pp. 287-303, doi:   10.1016/j.neuroimage.2013.01.04
%
% Author: Tim Mullen 10-5-2013, SCCN/INC, UCSD. 
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


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

g = arg_define([0 2],varargin, ...
        arg_norep('STUDY',mandatory,[],'STUDY object'),  ...
        arg_norep('ALLEEG',mandatory,[],'ALLEEG array'), ...
        arg({'connmethod','Estimator'},ConnNames{1},ConnNames,'Connectivity estimator to project'), ...
        arg({'timeRange','TimeRange'},timeRangeDef,[],'[Min Max] Time range to smooth (sec). Leave blank to use all time points','shape','row','type','denserealdouble'), ...
        arg({'freqRange','FrequencyRange'},freqRangeDef,[],'[Min Max] Frequency range to smooth (Hz). Leave blank to use all frequencies','type','expression','shape','row'), ...
        arg({'collapseTime','CollapseTime'},false,[],'Average across time before smoothing'), ...
        arg({'collapseFreqs','CollapseFreqs'},true,[],'Integrate across frequencies before smoothing'), ...
        arg({'rejOutsideBrain','RejectDipolesOutsideBrain'},true,[],'Reject dipoles outside the brain'), ...
        arg_sub({'graphMetric','GraphMetric'},{},@hlp_computeGraphMeasure,'Graph metric options','suppress',{'srcNodes','targNodes'}), ...
        arg({'enableCache','EnableCaching'},true,[],'Enable caching of MPT object. This could use a lot of memory, but accelerate subsequent calls to this function'), ...
        arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );
    
% basic error checking
if g.collapseTime && g.collapseFreqs
    error('At most one dimension can be collapsed. Otherwise, there is nothing to smooth!');
end
if any(arrayfun(@(x) ~isempty(hlp_checkeegset(x,{'conn'})), ALLEEG))
    error('You must first run SIFT analysis on all datasets in the STUDY');
end

waitbarstr = ['Creating MPT Object for ',g.graphMetric.graphMetric];
if g.verb==2
    % reset the color list
    hlp_getNextUniqueColor(hsv(10),'reset');
    multiWaitbar(waitbarstr,'Reset','Color',hlp_getNextUniqueColor);
elseif g.verb==1
    fprintf([waitbarstr '...\n']);
end

% create MPT object
if g.enableCache
    dipoleAndMeasure = hlp_microcache('grp_mpa',@pr.dipoleAndMeasureOfStudySIFT,g.STUDY, g.ALLEEG);
else
    dipoleAndMeasure = pr.dipoleAndMeasureOfStudySIFT(g.STUDY, g.ALLEEG);
end
% reject dipoles outside the brain
if g.rejOutsideBrain
    if g.verb, fprintf('Removing dipoles outside the brain\n'); end
    dipoleAndMeasure = dipoleAndMeasure.createSubsetInRelationToBrain();
end

if g.verb, fprintf('Computing graph measures...\n'); end
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
dipoleAndMeasure.time = nt*1000;
dipoleAndMeasure.frequency = nf;
dipoleAndMeasure.measureLabel = g.graphMetric.graphMetric;
dipoleAndMeasure.numberOfMeasureDimensions = ndims(gm)-1;

% store results and prepare outputs
STUDY = g.STUDY;
STUDY.measureprojection.sift.object =  dipoleAndMeasure;
varargout{1} = STUDY;
if nargout>1, varargout{2} = dipoleAndMeasure; end

if g.verb==2
    multiWaitbar(waitbarstr,'Close');
elseif g.verb==1, fprintf('done\n'); end


% Helper functions
% --------------------------------------------------------------
function [gm newtimes newfreqs] = ComputeGraphMeasure(Conn,g)
% compute graph measures from a Connectivity object

if ~iscell(g.connmethod)
    g.connmethod = {g.connmethod};
end

% collapse connectivity matrices
if g.collapseFreqs
    % collapse across freqs
    Conn = hlp_filterConns(Conn,'connmethods',g.connmethod, ...
        'method',{'freq','net'},'frange',g.freqRange,'freqdim',3,'timedim',4,'verb',0);
else
    % only select subset of freqs
    Conn = hlp_filterConns(Conn,'connmethods',g.connmethod, ...
        'method',{'freq','shrinkonly'}, ...
        'frange',g.freqRange, 'freqdim',3,'timedim',4,'verb',0);
end
if g.collapseTime
    % collapse across time
    Conn = hlp_filterConns(Conn,'connmethods',g.connmethod, ...
        'method',{'time','mean'},'trange',g.timeRange,'freqdim',3,'timedim',4,'verb',0);
else
    % only select subset of times
    Conn = hlp_filterConns(Conn,'connmethods',g.connmethod, ...
        'method',{'time','shrinkonly'}, ...
        'trange',g.timeRange,'freqdim',3,'timedim',4,'verb',0);
end
newtimes = Conn.erWinCenterTimes;
newfreqs = Conn.freqs;
% compute graph measure
gm = hlp_computeGraphMeasure('cmatrix',Conn.(g.connmethod{1}),g.graphMetric);









    
    
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
% gmProj = gmProj.createDomain(gmProj, 0.9, 0.05);
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

