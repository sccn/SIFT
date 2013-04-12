function [Stats PConn g] = stat_analyticStats(varargin)
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
% EEG:            EEG data structure containing connectivity structure in EEG.CAT.Conn
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

% initialization
PConn = [];

% extract some stuff from inputs for arg defaults
EEG = arg_extract(varargin,{'EEG','ALLEEG'},1);

if length(EEG)>1
    error('analytic stats currently available only for single-condition datasets');
end

if ~isempty(EEG)
    if ~isempty(EEG.CAT.Conn)
        Conn = EEG.CAT.Conn(1);
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
clear EEG Conn;

% argument definition
g = arg_define([0 Inf],varargin, ...
    arg_norep({'EEG','ALLEEG'},mandatory,[],'EEG structure.'), ...
    arg({'estimator','Estimator'},conndef,ConnNames,'Estimator for which to compute stats.','type','logical'), ...
    arg({'statistic','Statistic'},{'P-value','Threshold','ConfidenceInterval'},{'P-value','Threshold','ConfidenceInterval'},'Statistical quantities to return.','type','logical','cat','Statistic'), ...
    arg({'alpha','Alpha'},0.05,[0 1],'Significance level. This is used for significance thresholds (single-sided) and confidence intervals (two-sided). For example, a value of alpha=0.05 will produce p < alpha=0.05 p-value thresholds and (1-alpha)*100 = 95% confidence intervals.','cat','Statistic'), ...
    arg_subtoggle({'genpdf','ComputePDF'},[],...
        {arg({'numSamples','NumSamples'},250,[10 10000],'Number of samples from distribution')} ...
    ,'Generate samples from analytic PDF. This can be treated as a surrogate distribution for statistical comparisons','cat','Statistic'), ...
    arg({'verb','VerbosityLevel'},1,{int32(0) int32(1)},'Verbosity level. I.E., 0 = no command-line output, 1 = verbose output') ...
    );

if isempty(g.estimator), error('You must supply an estimator!'); end
try b=isempty(g.EEG.CAT.Conn); catch, b=1; end
if b, error('EEG must contain a CAT.Conn structure!'); end

for m=1:length(g.estimator)
    estimator = g.estimator{m};
    
    switch estimator
        case 'RPDC'
            
            if g.verb
                fprintf('Computing asymptotic statistics for RPDC. Please wait a moment...\n'); 
            end
            
            N = g.EEG.CAT.trials*(g.EEG.srate*g.EEG.CAT.MODEL.winlen);
            
            % specify the degrees of freedom
            if g.EEG.CAT.MODEL.morder == 1
                df = 1;
            elseif g.EEG.CAT.Conn.freqs(1)==0 ...
                || g.EEG.CAT.Conn.freqs(end)==g.EEG.srate/2
                % adjust degrees of freedom at end-points
                df = 2*ones(size(g.EEG.CAT.Conn.RPDC),'single');
                df(:,:,[1 end],:) = 1;
            else
                df = 2;
            end
%             df=fastif((g.EEG.CAT.MODEL.morder == 1), 1 , 2);
            
            if ismember('P-value',g.statistic)
                Stats.RPDC.pval = 1-chi2cdf(g.EEG.CAT.Conn.RPDC*N,df);
            end
            
            if ismember('Threshold',g.statistic)
                Stats.RPDC.thresh = chi2inv(1-g.alpha,df)/N;
                Stats.RPDC.thresh = Stats.RPDC.thresh(ones(size(g.EEG.CAT.Conn.RPDC)));
            end
            
            if ismember('ConfidenceInterval',g.statistic)
                Stats.RPDC.ci(1,:,:,:,:) = ncx2inv(g.alpha/2,df,g.EEG.CAT.Conn.RPDC*N)/N;    % lower bound
                Stats.RPDC.ci(2,:,:,:,:) = ncx2inv(1-g.alpha/2,df,g.EEG.CAT.Conn.RPDC*N)/N;  % upper bound
            end
            
            if g.genpdf.arg_selection && nargout > 1
               
                % check if we are likely to exceed available memory and
                % notify user.
                bytesAvail = hlp_getAvailableMemory('bytes');
                bytesReq   = 4*(2*numel(g.EEG.CAT.Conn.RPDC)*g.genpdf.numSamples);
                
                if bytesReq > bytesAvail
                    res = questdlg2(sprintf('This operation will require at least %5.5g MiB. It appears you may not have sufficient memory to carry out this operation. Do you want to continue?',bytesReq/(1024^2)),'Memory check','Yes','No','No');
                    if strcmpi(res,'no')
                        return;
                    end
                end

                % generate data from analytic PDF
                nd = ndims(g.EEG.CAT.Conn.RPDC);
                PConn.RPDC = ncx2rnd(df,repmat(g.EEG.CAT.Conn.RPDC*N,[ones(1,nd) g.genpdf.numSamples]))/N;
            end
            
        case 'nPDC'
            if ~isreal(g.EEG.CAT.Conn.nPDC)
                error('normalized PDC cannot be complex and must be magnitude-squared estimates. Please use hlp_absvalsq() first to obtain magnitude-squared estimates.');
            end
            
            if g.verb
                fprintf('Computing asymptotic statistics for nPDC. Please wait a moment...\n'); 
            end
            
            if g.EEG.CAT.MODEL.morder==1
                df    = 2;
                scale = 0.5;
            else
                df    = 1;
                scale = 1;
            end
            
            N = g.EEG.CAT.trials*(g.EEG.srate*g.EEG.CAT.MODEL.winlen);
            
            connmethods = [];
            if ~isfield(g.EEG.CAT.Conn,'pdc_denom'),
                connmethods = {'pdc_denom'};
            end
            if ~isfield(g.EEG.CAT.Conn,'Vpdc')
                connmethods = [connmethods 'Vpdc'];
            end
            
            if ~isempty(connmethods)
                % compute required estimates
                Conn = est_mvarConnectivity(g.EEG,g.EEG.CAT.MODEL,'connmethods',connmethods,'verb',g.verb,'freqs',g.EEG.CAT.Conn.freqs);
                if isempty(Conn)
                    Stats = [];
                    return;
                end
                g.EEG.CAT.Conn.Vpdc = Conn.Vpdc;
                g.EEG.CAT.Conn.pdc_denom = Conn.pdc_denom;
                clear Conn;
            end
            
            
            if ismember('P-value',g.statistic)
                Stats.nPDC.pval = 1-chi2cdf((g.EEG.CAT.Conn.nPDC.*g.EEG.CAT.Conn.pdc_denom*N/scale)./g.EEG.CAT.Conn.Vpdc , df);
            end
            
            if ismember('Threshold',g.statistic)
                Stats.nPDC.thresh = bsxfun(@rdivide,g.EEG.CAT.Conn.Vpdc.*chi2inv(1-g.alpha,df)/scale,N*g.EEG.CAT.Conn.pdc_denom);
            end
            
            if ismember('ConfidenceInterval',g.statistic)
                % TODO
                Stats.nPDC.ci = [];
%                 
%                 Stats.nPDC.ci(1,:,:,:,:) = ncx2inv(g.alpha/2,df,(g.EEG.CAT.Conn.nPDC.*g.EEG.CAT.Conn.pdc_denom*N/scale)./g.EEG.CAT.Conn.Vpdc);    % lower bound
%                 Stats.nPDC.ci(2,:,:,:,:) = ncx2inv(1-g.alpha/2,df,(g.EEG.CAT.Conn.nPDC.*g.EEG.CAT.Conn.pdc_denom*N/scale)./g.EEG.CAT.Conn.Vpdc);  % upper bound
            end
            
            if g.genpdf.arg_selection && nargout > 1
               
                % check if we are likely to exceed available memory and
                % notify user.
                bytesAvail = hlp_getAvailableMemory('bytes');
                bytesReq   = 4*(2*numel(g.EEG.CAT.Conn.nPDC)*g.genpdf.numSamples);
                
                if bytesReq > bytesAvail
                    res = questdlg2(sprintf('This operation will require at least %5.5g MiB. It appears you may not have sufficient memory to carry out this operation. Do you want to continue?',bytesReq/(1024^2)),'Memory check','Yes','No','No');
                    if strcmpi(res,'no')
                        return;
                    end
                end

                % generate data from analytic PDF
                nd = ndims(g.EEG.CAT.Conn.nPDC);
                PConn.nPDC = ncx2rnd(df,repmat(g.EEG.CAT.Conn.nPDC*N,[ones(1,nd) g.genpdf.numSamples]))/N;
            end
            
        case 'nDTF'
            
%             if g.verb
%                 fprintf('Computing asymptotic statistics for nDTF. Please wait a moment...\n'); 
%             end
            
            %             if ~isreal(g.EEG.CAT.Conn.nDTF)
            %                 error('normalized PDC cannot be complex and must be magnitude-squared estimates. Please use hlp_absvalsq() first to obtain magnitude-squared estimates.');
            %             end
            %
            %
            %             N = g.EEG.CAT.trials*(g.EEG.srate*g.EEG.CAT.MODEL.winlen);
            %
            %             connmethods = [];
            %             if ~isfield(g.EEG.CAT.Conn,'pdc_denom'),
            %                 connmethods = {'pdc_denom'};
            %             end
            %             if ~isfield(g.EEG.CAT.Conn,'Vpdc')
            %                 connmethods = [connmethods 'Vpdc'];
            %             end
            %
            %             if ~isempty(connmethods)
            %                 % compute required estimates
            %                 Conn = est_mvarConnectivity(g.EEG,g.EEG.CAT.MODEL,'connmethods',connmethods,'verb',g.verb,'freqs',g.EEG.CAT.Conn.freqs);
            %                 g.EEG.CAT.Conn.Vpdc = Conn.Vpdc;
            %                 g.EEG.CAT.Conn.pdc_denom = Conn.pdc_denom;
            %                 clear Conn;
            %             end
            %
            %
            %             if ismember('P-value',g.statistic)
            %                 Stats.nDTF.pval = 1-chi2cdf((g.EEG.CAT.Conn.nDTF.*g.EEG.CAT.Conn.pdc_denom*N)./g.EEG.CAT.Conn.Vpdc , 1);
            %             end
            %
            %             if ismember('Threshold',g.statistic)
            %                 Stats.nDTF.thresh = bsxfun(@rdivide,g.EEG.CAT.Conn.Vpdc.*chi2inv(1-g.alpha,1),N*g.EEG.CAT.Conn.pdc_denom);
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

% construct PConn object
if ~isempty(PConn) || nargout > 2
    PConn.winCenterTimes    = g.EEG.CAT.Conn.winCenterTimes;
    PConn.erWinCenterTimes  = g.EEG.CAT.Conn.erWinCenterTimes;
    PConn.freqs             = g.EEG.CAT.Conn.freqs;
    PConn.mode              = 'analytic';
    PConn.resampleTrialIdx  = [];
    connfields = hlp_getConnMethodNames(PConn);
    for m=1:length(connfields)
%         % check dimensions
%         szp = size(PConn.(connfields{m}));
%         szs = size(g.EEG.CAT.Conn.(connfields{m}));
%         [dummy dimidx] = setdiff(szp(1:end-1),szs);
%         if ~isempty(dimidx)
%             % a singleton dimension was squeezed out, restore it
%             PConn.(connfields{m}) = hlp_insertSingletonDim(PConn.(connfields{m}),dimidx+1);
%         end

        % insert singleton dimensions if necessary
        if length(PConn.winCenterTimes)==1
            PConn.(connfields{m}) = hlp_insertSingletonDim(PConn.(connfields{m}),4);
        end
        if length(PConn.freqs)==1
            PConn.(connfields{m}) = hlp_insertSingletonDim(PConn.(connfields{m}),3);
        end
    end
end
