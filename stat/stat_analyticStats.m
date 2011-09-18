function Stats = stat_analyticStats(varargin)
%
% Compute analytic alpha-significance thresholds, p-values, and confidence
% intervals for select connectivity estimators (RPDC, nPDC, nDTF, DTF)
%
% Note that nPDC stats are valid only for the normalized magnitude-squared PDC
%
% This function requires the Matlab statistics toolbox
%
% ===============================================
%    This function is under development and may
%    be unstable.
%    Please check sccn.uscd.edu/wiki/SIFT for 
%    updated version
% ===============================================
%
%
% Input           Information
% --------------------------------------------------------------------------------------------------------------------
% ALLEEG:         EEG data structure containing connectivity structure in EEG.CAT.Conn
%
% Optional        Information
% --------------------------------------------------------------------------------------------------------------------
% Estimator:      Estimator for which to compute stats
%                 Possible values: ''
%                 Default value  : 'n/a'
%                 Input Data Type: boolean
%
% Statistic:      Statistical quantities to return
%                 Possible values: 'P-value','Threshold','ConfidenceInterval'
%                 Default value  : 'P-value','Threshold','ConfidenceInterval'
%                 Input Data Type: boolean
%
% Alpha:          Significance level
%                 This is used for significance thresholds (single-sided) and confidence intervals (two-sided). For
%                 example, a value of alpha=0.05 will produce p < alpha=0.05 p-value thresholds and (1-alpha)*100 =
%                 95% confidence intervals.
%                 Input Range  : [0  1]
%                 Default value: 0.05
%                 Input Data Type: real number (double)
%
% VerbosityLevel: Verbosity level
%                 I.E., 0 = no command-line output, 1 = verbose output
%                 Possible values: 0,1
%                 Default value  : 1
%                 Input Data Type: real number (double)
%
% Output          Information
% --------------------------------------------------------------------------------------------------------------------
% Stats           Contains p-values, confidence interval, and
%                 other statistics for each connectivity
%                 estimator in EEG.CAT.Conn. 
%                 The dimensions for the stats matrices are generally 
%                 the same (or +1) as the dimensions of the matrix 
%                 for the original estimator. 
%                 Output format is as follows:
%                 Stats.<estimator>.[pval, ci, thresh, ...]
% 
%                               
%
% See Also: stat_surrogate()
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
ALLEEG = arg_extract(varargin,'ALLEEG',1);

if length(ALLEEG)>1
    error('analytic stats currently available only for individual datasets');
end

if ~isempty(ALLEEG)
    if ~isempty(ALLEEG.CAT.Conn)
        Conn = ALLEEG.CAT.Conn(1);
        ConnNames   = hlp_getConnMethodNames(Conn);
        ConnNames = intersect(ConnNames,{'RPDC','nPDC'});   %,'PDC','nDTF','DTF'
        conndef     = ConnNames;
    else
        ConnNames = {''};
        conndef = '';
    end
else
    ConnNames = {''};
    conndef = '';
end
clear ALLEEG Conn;

% argument definition
g = arg_define([0 Inf],varargin, ...
    arg_norep('ALLEEG',mandatory,[],'EEG structure.'), ...
    arg({'estimator','Estimator'},conndef,ConnNames,'Estimator for which to compute stats.','type','logical'), ...
    arg({'statistic','Statistic'},{'P-value','Threshold','ConfidenceInterval'},{'P-value','Threshold','ConfidenceInterval'},'Statistical quantities to return.','type','logical'), ...
    arg({'alpha','Alpha'},0.05,[0 1],'Significance level. This is used for significance thresholds (single-sided) and confidence intervals (two-sided). For example, a value of alpha=0.05 will produce p < alpha=0.05 p-value thresholds and (1-alpha)*100 = 95% confidence intervals.'), ...
    arg({'verb','VerbosityLevel'},1,{int32(0) int32(1)},'Verbosity level. I.E., 0 = no command-line output, 1 = verbose output') ...
    );

if isempty(g.estimator), error('You must supply an estimator!'); end
try b=isempty(g.ALLEEG.CAT.Conn); catch, b=1; end
if b, error('ALLEEG must contain a CAT.Conn structure!'); end


for m=1:length(g.estimator)
    estimator = g.estimator{m};
    
    switch estimator
        case 'RPDC'
            
            if g.verb
                fprintf('Computing asymptotic statistics for RPDC. Please wait a moment...\n'); 
            end
            
            N = g.ALLEEG.trials*(g.ALLEEG.srate*g.ALLEEG.CAT.MODEL.winlen);
            df=fastif((g.ALLEEG.CAT.MODEL.morder == 1), 1 , 2);
            
            if ismember('P-value',g.statistic)
                Stats.RPDC.pval = 1-chi2cdf(g.ALLEEG.CAT.Conn.RPDC*N,df);
            end
            
            if ismember('Threshold',g.statistic)
                Stats.RPDC.thresh = chi2inv(1-g.alpha,df)/N;
                Stats.RPDC.thresh = Stats.RPDC.thresh(ones(size(g.ALLEEG.CAT.Conn.RPDC)));
            end
            
            if ismember('ConfidenceInterval',g.statistic)
                Stats.RPDC.ci(1,:,:,:,:) = ncx2inv(g.alpha/2,df,g.ALLEEG.CAT.Conn.RPDC*N)/N;    % lower bound
                Stats.RPDC.ci(2,:,:,:,:) = ncx2inv(1-g.alpha/2,df,g.ALLEEG.CAT.Conn.RPDC*N)/N;  % upper bound
            end
            
        case 'nPDC'
            if ~isreal(g.ALLEEG.CAT.Conn.nPDC)
                error('normalized PDC cannot be complex and must be magnitude-squared estimates. Please use hlp_absvalsq() first to obtain magnitude-squared estimates.');
            end
            
            if g.verb
                fprintf('Computing asymptotic statistics for nPDC. Please wait a moment...\n'); 
            end
            
            if g.ALLEEG.CAT.MODEL.morder==1
                df    = 2;
                scale = 0.5;
            else
                df    = 1;
                scale = 1;
            end
            
            N = g.ALLEEG.trials*(g.ALLEEG.srate*g.ALLEEG.CAT.MODEL.winlen);
            
            connmethods = [];
            if ~isfield(g.ALLEEG.CAT.Conn,'pdc_denom'),
                connmethods = {'pdc_denom'};
            end
            if ~isfield(g.ALLEEG.CAT.Conn,'Vpdc')
                connmethods = [connmethods 'Vpdc'];
            end
            
            if ~isempty(connmethods)
                % compute required estimates
                Conn = est_mvarConnectivity(g.ALLEEG,g.ALLEEG.CAT.MODEL,'connmethods',connmethods,'verb',g.verb,'freqs',g.ALLEEG.CAT.Conn.freqs);
                g.ALLEEG.CAT.Conn.Vpdc = Conn.Vpdc;
                g.ALLEEG.CAT.Conn.pdc_denom = Conn.pdc_denom;
                clear Conn;
            end
            
            
            if ismember('P-value',g.statistic)
                Stats.nPDC.pval = 1-chi2cdf((g.ALLEEG.CAT.Conn.nPDC.*g.ALLEEG.CAT.Conn.pdc_denom*N/scale)./g.ALLEEG.CAT.Conn.Vpdc , df);
            end
            
            if ismember('Threshold',g.statistic)
                Stats.nPDC.thresh = bsxfun(@rdivide,g.ALLEEG.CAT.Conn.Vpdc.*chi2inv(1-g.alpha,df)/scale,N*g.ALLEEG.CAT.Conn.pdc_denom);
            end
            
            if ismember('ConfidenceInterval',g.statistic)
                % TODO
                Stats.nPDC.ci = [];
%                 
%                 Stats.nPDC.ci(1,:,:,:,:) = ncx2inv(g.alpha/2,df,(g.ALLEEG.CAT.Conn.nPDC.*g.ALLEEG.CAT.Conn.pdc_denom*N/scale)./g.ALLEEG.CAT.Conn.Vpdc);    % lower bound
%                 Stats.nPDC.ci(2,:,:,:,:) = ncx2inv(1-g.alpha/2,df,(g.ALLEEG.CAT.Conn.nPDC.*g.ALLEEG.CAT.Conn.pdc_denom*N/scale)./g.ALLEEG.CAT.Conn.Vpdc);  % upper bound
            end
            
            
        case 'nDTF'
            
%             if g.verb
%                 fprintf('Computing asymptotic statistics for nDTF. Please wait a moment...\n'); 
%             end
            
            %             if ~isreal(g.ALLEEG.CAT.Conn.nDTF)
            %                 error('normalized PDC cannot be complex and must be magnitude-squared estimates. Please use hlp_absvalsq() first to obtain magnitude-squared estimates.');
            %             end
            %
            %
            %             N = g.ALLEEG.trials*(g.ALLEEG.srate*g.ALLEEG.CAT.MODEL.winlen);
            %
            %             connmethods = [];
            %             if ~isfield(g.ALLEEG.CAT.Conn,'pdc_denom'),
            %                 connmethods = {'pdc_denom'};
            %             end
            %             if ~isfield(g.ALLEEG.CAT.Conn,'Vpdc')
            %                 connmethods = [connmethods 'Vpdc'];
            %             end
            %
            %             if ~isempty(connmethods)
            %                 % compute required estimates
            %                 Conn = est_mvarConnectivity(g.ALLEEG,g.ALLEEG.CAT.MODEL,'connmethods',connmethods,'verb',g.verb,'freqs',g.ALLEEG.CAT.Conn.freqs);
            %                 g.ALLEEG.CAT.Conn.Vpdc = Conn.Vpdc;
            %                 g.ALLEEG.CAT.Conn.pdc_denom = Conn.pdc_denom;
            %                 clear Conn;
            %             end
            %
            %
            %             if ismember('P-value',g.statistic)
            %                 Stats.nDTF.pval = 1-chi2cdf((g.ALLEEG.CAT.Conn.nDTF.*g.ALLEEG.CAT.Conn.pdc_denom*N)./g.ALLEEG.CAT.Conn.Vpdc , 1);
            %             end
            %
            %             if ismember('Threshold',g.statistic)
            %                 Stats.nDTF.thresh = bsxfun(@rdivide,g.ALLEEG.CAT.Conn.Vpdc.*chi2inv(1-g.alpha,1),N*g.ALLEEG.CAT.Conn.pdc_denom);
            %             end
            %
            %             if ismember('ConfidenceInterval',g.statistic)
            %                 % TODO
            %                 Stats.nDTF.ci = [];
            %             end
            
            
        case 'DTF'
            % TODO
            
    end
    
    
    Stats.alpha = g.alpha;
end

