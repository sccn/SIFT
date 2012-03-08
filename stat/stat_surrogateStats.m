function [Stats ConnMean] = stat_surrogateStats(varargin)
%
% Return surrogate statistics based on a surrogate distribution (bootstrap, jacknife, etc).
% Each surrogate distribution should approximate the distribution of the
% estimator.
%
%
% Stats are performed across the 5th dimension (boostrap samples or subjects).
% This function calls the statcond() routine written by Arnaud Delorme.
%
% ===============================================
%    This function is under development and may
%    be unstable.
%    Please check sccn.uscd.edu/wiki/SIFT for
%    updated version
% ===============================================
%
% Input                         Information
% ----------------------------------------------------------------------------------------------------------------
% BootstrapConnectivity         Surrogate connectivity structure as returned
%                               by stat_surrogate()
%
% Optional                      Information
% ----------------------------------------------------------------------------------------------------------------
% NullDistribution              Null distribution structure as returned by
%                               stat_surrogate()
%
%
% ConnectivityMethods:          Connectivity estimator(s) to bootstrap
%                               All connectivity estimators must be named exactly as they appear in EEG.CAT.Conn
%                               Possible values: ''
%                               Default value  : 'n/a'
%                               Input Data Type: boolean
%
% MultipleComparisonCorrection: Correction for multiple comparisons.
%                               'numvars' does a bonferonni correction
%                               considering M^2 indep. degrees of freedom,
%                               where M is the dimension of the VAR model
%                               (i.e. number of channels)
%                               Possible values: 'none','fdr','bonferonni','numvars'
%                               Default value  : 'fdr'
%                               Input Data Type: string
%
% ComputeConfIntervals:         Compute confidence intervals
%                               Input Range  : Unrestricted
%                               Default value: 0
%                               Input Data Type: boolean
%
%     | Alpha:                  Significance level (two sided)
%                               A value of 0.95 will produce 95% confidence intervals.
%                               Input Range  : [0  1]
%                               Default value: 0.95
%                               Input Data Type: real number (double)
%
% statcondargs:                 List of paired arguments for statcond()
%                               Possible values: Unrestricted
%                               Default value  : 'mode','perm'
%                               Input Data Type: any evaluable Matlab expression.
%
% VerbosityLevel:               Verbosity level. 0 = no output, 1 = verbose output
%                               Possible values: 0,1
%                               Default value  : 1
%                               Input Data Type: real number (double)
%
% Output                        Information
% ----------------------------------------------------------------------------------------------------------------
% Stats                         Contains p-values, confidence interval, and
%                               other statistics for each connectivity
%                               estimator in BootstrapConnectivity. The
%                               dimensions for the stats matrices are generally
%                               the same (or +1) as the dimensions of the matrix
%                               for the original estimator.
%                               Output format is as follows:
%                               Stats.<estimator>.[pval, ci, thresh, ...]
% ConnMean                      Contains the mean of the bootstrap
%                               distribution for each connectivity
%                               estimator (if using bootstrap methods).
%                               This does not apply to PhaseRand methods
%
%
% See Also: stat_surrogate(), stat_analyticStats(), statcond()
%
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 5.
%   Available at: http://www.sccn.ucsd.edu/wiki/SIFT
%
% Author: Tim Mullen, 2011, SCCN/INC, UCSD.
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
PConn = arg_extract(varargin,'PConn',2);
if ~isempty(PConn(1))
    ConnNames   = hlp_getConnMethodNames(PConn(1));
    conndef     = ConnNames;
else
    ConnNames = {''};
    conndef = '';
end

if length(PConn)>2
    error('SIFT:stat_surrogateStats','A maximum of two datasets can be compared statistically');
elseif length(PConn)==2
    statTestDef = {'Hab'};
elseif length(PConn)==1
    if strcmpi(PConn.mode,'phaserand')
        statTestDef = {'Hnull'};
    else
        statTestDef = {'Hbase'};
    end
end


% statTestDef = {'Hnull','Hbase','Hab'};

g = arg_define([0 Inf],varargin, ...
    arg_norep({'PConn','BootstrapConnectivity'},mandatory,[],'Boostrap connectivity structure.'), ...
    arg_norep({'null','NullDistribution'},[],[],'Null distribution structure.'), ...
    arg_subswitch({'statTest','StatisticalTest'},statTestDef,...
    {'Hnull', ...
    {arg_norep({'dummy1'},[],[],'dummy')} ...
    'Hbase', ...
    { ...
    arg({'baseline','Baseline'},[],[],'Time range of baseline [Min Max] (sec). Will subtract baseline from each point.','shape','row','type','denserealdouble'), ...
    }, ...
    'Hab', ...
    {arg_norep({'dummy2'},[],[],'dummy')}, ...
    }, 'Statistical test to perform. Hnull: compare PConn to null distribution (must be provided in ''NullDistribution'' argument) (one-sided, unpaired). Hbase: Compare each sample in PConn to baseline mean (two-sided, unpaired). Hab: compute significance for two conditions (PConn(1)-PConn(2)) (two-sided, paired).'), ...
    arg({'connmethods','ConnectivityMethods'},conndef,ConnNames,'Connectivity estimator(s) to bootstrap. All connectivity estimators must be named exactly as they appear in EEG.CAT.Conn','type','logical'), ...
    arg({'mcorrection','MultipleComparisonCorrection'},'fdr',{'none','fdr','bonferonni','numvars'},'Correction for multiple comparisons. Note: ''numvars'' does a bonferonni correction considering M^2 indep. degrees of freedom, where M is the dimension of the VAR model (i.e. number of channels)'), ...
    arg({'computeci','ConfidenceIntervals'},true, [],'Compute empirical confidence intervals.'), ...
    arg({'alpha','Alpha'},0.05,[0 1],'Significance level. This is used for significance thresholds (one- or two-sided) and confidence intervals (two-sided). For example, a value of alpha=0.05 will produce p < alpha=0.05 thresholds and (1-alpha)*100 = 95% confidence intervals.'), ...
    arg('statcondargs',{'mode','perm'},{},'List of paired arguments for statcond()','type','expression','shape','row'), ...
    arg({'verb','VerbosityLevel'},1,{int32(0) int32(1)},'Verbosity level. 0 = no output, 1 = verbose output') ...
    );


clear PConn;

if ~ismember('mode',lower(g.statcondargs))
    g.statcondargs = [g.statcondargs 'mode', 'perm'];
end

if g.verb==0
    g.statcondargs = [g.statcondargs 'verbose', 'off'];
else
    g.statcondargs = [g.statcondargs 'verbose', 'on'];
end

if ~isempty(g.null)
    surrogateMode = g.null.mode;
else
    surrogateMode = g.PConn.mode;
end

for m=1:length(g.connmethods)
    
    switch g.statTest.arg_selection
        case 'Hnull'
            % Our null hypothesis is Conn(i,j)=0
            % perform one-sided statistical test against null hypothesis
            
            if length(g.PConn)>1
                error('SIFT:stat_surrogateStats','Please select a single condition for Hnull test.');
            end
            
            if strcmpi(surrogateMode,'phaserand')
                % We are testing with respect to a phase-randomized null
                % distribution. A p-value for rejection of the null hypothesis
                % can be obtained by computing the probability that the
                % observed connectivity is a random sample from the null distribution
                
                if isempty(g.null)
                    error('SIFT:stat_surrogateStats','You must provide a null distribution');
                end
                
                if g.verb
                    fprintf('\nTesting estimator %s\n',g.connmethods{m});
                    fprintf('Computing statistics against null hypothesis C(i,j)=0\n');
                    fprintf('Stats are based on %s distribution\n',g.PConn.mode);
                end
                
                if ndims(g.PConn.(g.connmethods{m}))==ndims(g.null.(g.connmethods{m}))
                    % a bootstrap distribution for the estimator was provided so
                    % we need to compute the mean of the bootstrap
                    % distribution for comparison to null distribution
                    g.PConn.(g.connmethods{m}) = stat_getDistribMean(g.PConn.(g.connmethods{m}));
                end
                
                [statval, df, Stats.(g.connmethods{m}).pval] = statcond( { }, 'mode','perm', 'surrog', g.null.(g.connmethods{m}), 'stats', g.PConn.(g.connmethods{m}), g.statcondargs{:},'tail','one');
                
                
                % compute the threshold based on 1-alpha percentile of null
                % distribution
                Stats.(g.connmethods{m}).thresh = prctile(g.null.(g.connmethods{m}),100-100*g.alpha,ndims(g.null.(g.connmethods{m})));
                
                % there are no confidence intervals to estimate here
                Stats.(g.connmethods{m}).ci = [];
            else
                % we are testing w.r.t. the approximate distribution of
                % the estimator itself. A p-value for rejection of the null
                % hypothesis can be obtained by computing the probability
                % that a sample from the estimator's distribution is less
                % than or equal to zero
                % univariate nonparametric test
                if g.verb
                    fprintf('\nTesting estimator %s\n',g.connmethods{m});
                    fprintf('Computing statistics against null hypothesis C(i,j)=0\n');
                    fprintf('Stats are based on %s distribution\n',g.PConn.mode);
                end
                
                sz = size(g.PConn.(g.connmethods{m}));
                [statval, df, Stats.(g.connmethods{m}).pval] = statcond( { }, 'mode','perm', 'surrog', g.PConn.(g.connmethods{m}), 'stats', zeros(sz(1:end-1)), g.statcondargs{:},'tail','one');
                Stats.(g.connmethods{m}).pval = 1-Stats.(g.connmethods{m}).pval;
                
                % compute two-sided confidence intervals
                % NOTE: need to check whether it's 2*g.alpha or just
                % g.alpha (should be twice that used for one-sided p-value
                % test above to ensure appropriate overlap of lower CI with
                % zero when non-signficant)
                Stats.(g.connmethods{m}).ci = stat_computeCI(g.PConn.(g.connmethods{m}),g.alpha,'both');
            end
            
        case 'Hab'
            % For conditions A and B, the null hypothesis is either
            % A(i,j)<=B(i,j), for a one-sided test, or
            % A(i,j)=B(i,j), for a two-sided test
            % A p-value for rejection of the null hypothesis can be
            % obtained by taking the difference of the distributions
            % computing the probability
            % that a sample from the difference distribution is non-zero
            % This is a paired nonparametric test
            if length(g.PConn)~=2
                error('SIFT:stat_surrogateStats','BootstrapConnectivity object must have two conditions for difference test');
            elseif size(g.PConn(1).(g.connmethods{m}))~=size(g.PConn(2).(g.connmethods{m}))
                error('SIFT:stat_surrogateStats','BootstrapConnectivity matrices must be of equal dimension for both conditions');
            end
            
            if g.verb
                fprintf('\nTesting estimator %s\n',g.connmethods{m});
                fprintf('Computing statistics against null hypothesis A(i,j)=B(i,j)\n');
                fprintf('This is a two-sided test for significant differences between conditions\n');
                fprintf('Stats are based on %s distributions\n',g.PConn(1).mode);
            end
            
            sz = size(g.PConn(1).(g.connmethods{m}));
            Pdiff = g.PConn(1).(g.connmethods{m})-g.PConn(2).(g.connmethods{m});
            [statval, df, Stats.(g.connmethods{m}).pval] = statcond( { }, 'mode','perm', 'surrog', Pdiff, 'stats', zeros(sz(1:end-1)), g.statcondargs{:},'tail','both');
            Stats.(g.connmethods{m}).pval = Stats.(g.connmethods{m}).pval;
            
            % compute two-sided confidence intervals
            % NOTE: need to check whether it's 2*g.alpha or just
            % g.alpha (should be twice that used for one-sided p-value
            % test above to ensure appropriate overlap of lower CI with
            % zero when non-signficant)
            Stats.(g.connmethods{m}).ci = stat_computeCI(Pdiff,g.alpha,'both');
            
        case 'Hbase'
            % For conditions A, the null hypothesis is
            % C(i,j)=baseline_mean(C). This is a two-sided test.
            % A p-value for rejection of the null hypothesis can be
            % obtained by obtaining the distribution of the difference from
            % baseline mean and computing the probability
            % that a sample from this distribution is non-zero
            % This is a paired nonparametric test
            
            if length(g.PConn)>1
                error('SIFT:stat_surrogateStats','Please select a single condition for Hbase test.');
            end
            
            if length(g.PConn.winCenterTimes)==1
                error('SIFT:stat_surrogateStats','Estimator must be time-varying to compute deviation from temporal baseline');
            end
            
            if g.verb
                fprintf('\nTesting estimator %s\n',g.connmethods{m});
                fprintf('Computing statistics against null hypothesis C(i,j)=baseline_mean(C)\n');
                fprintf('This is a two-sided test for significant deviation from baseline\n');
                fprintf('Stats are based on %s distributions\n',g.PConn.mode);
            end
            
            % compute the baseline deviation distribution
            % (Conn - baseline_mean(Conn))
            Pdiff = g.PConn.(g.connmethods{m}) - stat_getBaselineDistrib(g.PConn.(g.connmethods{m}),g.statTest.baseline,g.PConn.erWinCenterTimes);
            sz = size(Pdiff);
            
            [statval, df, Stats.(g.connmethods{m}).pval] = statcond( { }, 'mode','perm', 'surrog', Pdiff, 'stats', zeros(sz(1:end-1)), g.statcondargs{:},'tail','both');
            
            
            % compute two-sided confidence intervals
            % NOTE: need to check whether it's 2*g.alpha or just
            % g.alpha (should be twice that used for one-sided p-value
            % test above to ensure appropriate overlap of lower CI with
            % zero when non-signficant)
            Stats.(g.connmethods{m}).ci = stat_computeCI(Pdiff,g.alpha,'both');
            
            Stats.baseline = g.statTest.baseline;
    end
    
    % Correct for multiple comparisons
    switch g.mcorrection
        case 'fdr'
            Stats.(g.connmethods{m}).pval = fdr(Stats.(g.connmethods{m}).pval);
        case 'bonferonni'
            Stats.(g.connmethods{m}).pval = Stats.(g.connmethods{m}).pval * numel(Stats.(g.connmethods{m}).pval);
        case 'numvars'
            Stats.(g.connmethods{m}).pval = Stats.(g.connmethods{m}).pval * size(Stats.(g.connmethods{m}).pval,1)^2;
    end
    
    % check if a singleton dimension was squeezed out and, if so,
    % restore the singleton dim
    
    % check pval
    if isfield(Stats.(g.connmethods{m}),'pval')
        szp = size(g.PConn.(g.connmethods{m}));
        szs = size(Stats.(g.connmethods{m}).pval);
        [dummy dimidx] = setdiff(szp(1:end-1),szs);
        if ~isempty(dimidx)
            % a singleton dimension was squeezed out, restore it
            Stats.(g.connmethods{m}).pval = hlp_insertSingletonDim(Stats.(g.connmethods{m}).pval,dimidx+1);
        end
    end
    
    % check ci
    if isfield(Stats.(g.connmethods{m}),'ci')
        szp = size(g.PConn.(g.connmethods{m}));
        szs = size(Stats.(g.connmethods{m}).ci);
        [dummy dimidx] = setdiff(szp(1:end-1),szs(2:end));
        if ~isempty(dimidx)
            % a singleton dimension was squeezed out, restore it
            Stats.(g.connmethods{m}).ci = hlp_insertSingletonDim(Stats.(g.connmethods{m}).ci,dimidx+1);
        end
    end
    
    % check thresh
    if isfield(Stats.(g.connmethods{m}),'thresh')
        szp = size(g.PConn.(g.connmethods{m}));
        szs = size(Stats.(g.connmethods{m}).thresh);
        [dummy dimidx] = setdiff(szp(1:end-1),szs);
        if ~isempty(dimidx)
            % a singleton dimension was squeezed out, restore it
            Stats.(g.connmethods{m}).thresh = hlp_insertSingletonDim(Stats.(g.connmethods{m}).thresh,dimidx+1);
        end
    end
    
end


statcondargs = hlp_varargin2struct(g.statcondargs);
Stats.mode = statcondargs.mode;
Stats.correction = g.mcorrection;
Stats.alpha = g.alpha;

% return the mean of each distribution
if nargout>1
    if strcmpi(Stats.mode,'phaserand')
        % PConn should be the estimator itself
        ConnMean = [];
    else
        ConnMean = stat_getDistribMean(g.PConn);
    end
end


end





% ---------------------------
% | OLD SCRAPS OF CODE
% ---------------------------
%
%
%     if strcmpi('param',g.statcondargs(find(ismember(g.statcondargs,'mode'))+1))
%         % Perform parametric statistics
%
%         if ~isempty(g.null)
%             % compare connectivity to a null distribution
%             [statval, df, Stats.(g.connmethods{m}).pval] = statcond( {g.PConn.(g.connmethods{m}) g.null.(g.connmethods{m})} , g.statcondargs{:});
%         elseif length(g.PConn)==1
%             % univariate parametric test
%             [statval, df, Stats.(g.connmethods{m}).pval] = statcond( g.PConn.(g.connmethods{m}) , g.statcondargs{:});
%         else
%             % two-sample parametric difference test
%             [statval, df, Stats.(g.connmethods{m}).pval] = statcond( {g.PConn.(g.connmethods{m})} , g.statcondargs{:});
%         end
%     else
%         % Perform non-parametric statistics
%
%         if length(g.PConn)==2
%             % difference test
%             PConn = g.PConn(1).(g.connmethods{m})-g.PConn(2).(g.connmethods{m});
%             sz = size(PConn);
%         end
%
%         if ~isempty(g.null)
%             % compare connectivity to null distribution (e.g. mode=PhaseRand)
%             [statval, df, Stats.(g.connmethods{m}).pval] = statcond( { }, 'mode','perm', 'surrog', g.null.(g.connmethods{m}), 'stats', g.PConn.(g.connmethods{m}), g.statcondargs{:});
%         elseif length(g.PConn)==1
%             % univariate nonparametric test
%             sz = size(g.PConn.(g.connmethods{m}));
%             [statval, df, Stats.(g.connmethods{m}).pval] = statcond( { }, 'mode','perm', 'surrog', g.PConn.(g.connmethods{m}), 'stats', zeros(sz(1:end-1)), g.statcondargs{:});
%         else
%             % nonparametric difference test
%
%             [statval, df, Stats.(g.connmethods{m}).pval] = statcond( { }, 'mode','perm', 'surrog', PConn, 'stats', zeros(sz(1:end-1)), g.statcondargs{:});
%         end
%     end
