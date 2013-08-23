function [Conn PConn MCMC_LastState] = stat_bsplSmoother(varargin)
% PConn contains the distribution of the estimators stored in the last dimension
% This can be used with stat_surrogateStats to compute confidence
% intervals, p-values, etc in precisely the same manner as bootstrap distributions

% Author: Wes Thompson and Tim Mullen

% extract some stuff from inputs for arg defaults
Conn = arg_extract(varargin,'Conn',1);

if ~isempty(Conn)
    Conn            = Conn(1);
    ConnNames       = hlp_getConnMethodNames(Conn);
    freqRangeDef    = [Conn.freqs(1) Conn.freqs(end)];
    timeRangeDef    = [Conn.erWinCenterTimes(1) Conn.erWinCenterTimes(end)];
    clear Conn;
else
    ConnNames = {''};
    [freqRangeDef, timeRangeDef] = deal([]);
end

g = arg_define([0 Inf], varargin, ...
    arg_norep('Conn',mandatory), ...
    arg_nogui({'bsplmodel','BSplineModel'},struct([]),[],'Optional bivariate B-Spline model. This can be computed via stat_bsplsm_mkspl().'), ...
    arg_nogui({'MCMC_LastState'},struct([]),[],'Structure containing state of Gibbs sampler. This will be used to initialize the sampler'), ...
    arg({'connmethods','Estimator'},ConnNames,ConnNames,'Connectivity estimators to smooth','type','logical'), ...
    arg({'timeRange','TimeRange'},timeRangeDef,[],'[Min Max] Time range to smooth (sec). Leave blank to use all time points','shape','row','type','denserealdouble'), ...
    arg({'freqRange','FrequencyRange'},freqRangeDef,[],'[Min Max] Frequency range to smooth (Hz). Leave blank to use all frequencies','type','expression','shape','row'), ...
    arg({'collapseTime','CollapseTime'},false,[],'Average across time before smoothing'), ...
    arg({'collapseFreqs','CollapseFreqs'},true,[],'Integrate across frequencies before smoothing'), ...
    arg({'timeKnots','TimeKnots'},5,[],'Positions of spline knots along time dimension (sec). If a scalar, q, then q knots are evenly spaced from first to last timepoint.  A good heuristic is one knot every 5%','shape','row'), ...
    arg({'freqKnots','FreqKnots'},5,[],'Positions of spline knots along frequency dimension (Hz). If a scalar, q, then knots are evenly spaced from first to last frequency. A good heuristic is one knot every 5%','shape','row'), ...
    arg({'fpcaBasisDim','FPCABasisDim'},4,[0 Inf],'Number of FPCA basis functions. (univariate smoothing only)'), ...
    arg({'smoothingLayout','MatrixElementsToSmooth'},{'diagonals','off-diagonals'},{'off-diagonals'},'Which parts of the matrix to smooth. Diagonals (e.g. auto-connectivity) and off-diagonals (e.g. cross-connectivity) will be smoothed separately','type','logical'), ...
    arg({'nMCMCiters','NumMcmcIters','niter','niters'},1000, [1 Inf], 'Number of MCMC iterations for spline fitting'), ...
    arg({'burnInFraction','BurnInFractionForMCMC'},0.5,[0 0.99],'Fraction of initial MCMC samples to discard (burn in period). The number of MCMC samples is taken to be the smaller of MCMCitersOffDiag and MCMCitersDiag.'), ...
    arg({'normlog','NormalizeAndLogTranform'},true,[],'Transform data before smoothing. Normalize across last dim (time), add 1 and take logarithm. Inverse transform is applied after smoothing'), ...
    arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );

if isempty(g.smoothingLayout)
    error('SIFT:stat_bsplSmoother:badInput', ...
          'You must choose some matrix elements to smooth. See ''MatrixElementsToSmooth'' option'); 
end

if g.collapseTime && g.collapseFreqs
    error('At most one dimension can be collapsed');
end

numBurnInSamples = floor(g.burnInFraction*g.nMCMCiters);
g.niterToKeep    = g.nMCMCiters-numBurnInSamples;

if g.verb,
    fprintf(['I will discard %d burn-in samples.\n' ...
        'The distribution of the estimator will have %d samples\n'], ...
        numBurnInSamples,g.niterToKeep);
end

g.retMCMCLastState = (nargout > 1);

% collapse connectivity matrices
% -------------------------------------------------------------------------
Conn = g.Conn;
g = rmfield(g,'Conn');
bsplmodel = struct([]);

if g.collapseFreqs
    % collapse across freqs
    Conn = hlp_filterConns(Conn,'connmethods',g.connmethods, ...
        'method',{'freq','net'},'frange',g.freqRange,'freqdim',3,'timedim',4,'verb',g.verb);
else
    % only select subset of freqs
    Conn = hlp_filterConns(Conn,'connmethods',g.connmethods, ...
        'method',{'freq','shrinkonly'}, ...
        'frange',g.freqRange, 'freqdim',3,'timedim',4,'verb',g.verb);
end

if g.collapseTime
    % collapse across time
    Conn = hlp_filterConns(Conn,'connmethods',g.connmethods, ...
        'method',{'time','mean'},'trange',g.timeRange,'freqdim',3,'timedim',4,'verb',g.verb);
else
    % only select subset of times
    Conn = hlp_filterConns(Conn,'connmethods',g.connmethods, ...
        'method',{'time','shrinkonly'}, ...
        'trange',g.timeRange,'freqdim',3,'timedim',4,'verb',g.verb);
end


% smooth the time-series
% -------------------------------------------------------------------------
if isscalar(Conn.freqs) || isscalar(Conn.erWinCenterTimes)
    % univariate (1D) smoothing
    [Conn PConn] = univariate_smooth(Conn,g);
else
    % bivariate (2D) smoothing
    [Conn PConn] = bivariate_smooth(Conn,g);
    
end

PConn.mode      = 'bspline';
PConn.options   = g;
Conn.options    = g;


%% univariate_smooth()
%  Perform univariate smoothing (freq x causality or time x causality)
%  Smoothing is achieved via a Bayesian (Monte Carlo) B-spline smoother
%  with optional dimensionality reduction using Functional PCA
%  ------------------------------------------------------------------------

    function [Conn PConn] = univariate_smooth(Conn,g)
        
        smoothDiags     = ismember_bc('diagonals',g.smoothingLayout);
        smoothOffDiags  = ismember_bc('off-diagonals',g.smoothingLayout);
        
        ntimes = length(Conn.erWinCenterTimes);
        nfreqs = length(Conn.freqs);
        
        % obtain indices of knot positions
        if isscalar(Conn.freqs)
            % frequencies were collapsed, smooth across time
            
            if isscalar(g.timeKnots)
                % determine positions of a user-specified
                % number of knots with equal spacing
                order = 4;
                knots = quantile(1:ntimes,(0:1:(g.timeKnots-order+1)) ...
                                           / (g.timeKnots-order+1));
            else
                % user provided knot locations
                knots = getindex(Conn.erWinCenterTimes,g.timeKnots);
            end
            
        else % times were collapsed, smooth across frequency
            
            if isscalar(g.freqKnots)
                % determine positions of a user-specified
                % number of knots with equal spacing
                order = 4;
                knots = quantile(1:nfreqs,(0:1:(g.freqKnots-order+1)) ...
                                          / (g.freqKnots-order+1));
            else
                % user provided knot locations
                knots = getindex(Conn.freqs,g.freqKnots);
            end
            
        end
        
        % copy all supplementary fields into PConn
        PConn = hlp_splitstruct(Conn,setdiff_bc(fieldnames(Conn),hlp_getConnMethodNames(Conn)));
        
        % smooth each estimator in Conn
        % -----------------------------------------------------------------
        for m=1:length(g.connmethods)
            
            connmethod = g.connmethods{m};
            
            if g.verb
                waitbarstr = sprintf('Applying MCMC smoothing to %s',connmethod);
                multiWaitbar(waitbarstr,0,'Color',[0.8 0.0 0.1]);   
            end
            
            if g.normlog
                [Conn.(connmethod) mu sigma] = transform(double(Conn.(connmethod)));
            end
            
            % initalize smoothed conn matrix. last dimension contains 
            % independent realizations from MCMC iterations
            PConn.(connmethod) = zeros([size(Conn.(connmethod)) g.niterToKeep]);
            
            [fit_distrib MCMC_LastState] = stat_ubsplsm( ...
                                          {squeeze(Conn.(connmethod))}, ...
                                          knots,g.fpcaBasisDim,         ...
                                          g.nMCMCiters,g.niterToKeep,   ...
                                          smoothDiags,smoothOffDiags);
            
            if isscalar(PConn.freqs)
                % time was smoothed, insert singleton 'frequency' dimension
                PConn.(connmethod) = hlp_insertSingletonDim(fit_distrib{1},3);
            elseif isscalar(PConn.erWinCenterTimes)
                % frequency was smoothed
                PConn.(connmethod) = fit_distrib{1};
            end
            
            if g.normlog
                PConn.(connmethod) = untransform(PConn.(connmethod),mu,sigma);
            end
            
            % replace original connectivity estimates with mean smoothed 
            % connectivity estimates
            Conn.(connmethod) = mean(PConn.(connmethod),ndims(PConn.(connmethod)));
            
            if g.verb
                % cleanup waitbar
                multiWaitbar(waitbarstr, 'Close');
            end
            
        end
    end

%% subfunction bivariate_smooth()

    function [Conn PConn bsplmodel] = bivariate_smooth(Conn,g)
        % perform bivariate smoothing (time x freq x causality)
        
        skipDiags = ~ismember_bc('diagonals',g.smoothingLayout);
        
        ntimes = length(Conn.erWinCenterTimes);
        nfreqs = length(Conn.freqs);
        
        % obtain indices of knot positions
        if isscalar(g.timeKnots)
            % equally-spaced knots
            g.timeKnots = quantile(1:ntimes,(0:1:(g.timeKnots-1)) / (g.timeKnots-1));
            
        else
            g.timeKnots = getindex(Conn.erWinCenterTimes,g.timeKnots);
        end
        
        if isscalar(g.freqKnots)
            % equally-spaced knots
            g.freqKnots = quantile(1:nfreqs,(0:1:(g.freqKnots-1)) / (g.freqKnots-1));
        else
            g.freqKnots = getindex(Conn.freqs,g.freqKnots);
        end
        
        % * construct the b-spline basis functions
        
        if ~isempty(g.bsplmodel)
            % spline model already computed
            bsplmodel = g.bsplmodel;
        else
            % obtain the b-spline basis functions
            bsplmodel = stat_bsplsm_mkspl(ntimes,nfreqs,g.timeKnots,g.freqKnots,g.verb);
        end
        
        % * copy all supplementary fields into PConn
        PConn = hlp_splitstruct(Conn,setdiff_bc(fieldnames(Conn),hlp_getConnMethodNames(Conn)));
        
        
        % * smooth connectivity values
        for m=1:length(g.connmethods)
            
            connmethod = g.connmethods{m};
            
            if g.verb
                h=waitbar(0,sprintf('Applying MCMC smoothing to %s',connmethod));
            end
            
            C=Conn.(connmethod);
            
            nchs = size(C,1);
            
            % initalize smoothed conn matrix
            % last dimension contains independent
            % realizations from MCMC iterations
            PConn.(connmethod) = zeros([size(C) g.niterToKeep]);
            
            % smooth each connectivity pair independently
            % (note, this can be parallelized)
            curiter = 0;
            numPairs = nchs^2 - skipDiags*nchs;
            
            for i=1:nchs
                for j=1:nchs
                    
                    if skipDiags && i==j
                        continue;
                    end
                    
                    if g.verb
                        waitbar(curiter/numPairs,h,sprintf('Applying MCMC smoothing to %s (%d/%d)',connmethod,curiter,numPairs));
                        curiter = curiter+1;
                    end
                    
                    % smooth this C(i,j) connectivity pair
                    % NOTE: may need to multiply C by some large value (e.g. 1000)
                    PConn.(connmethod)(i,j,:,:,:) = stat_bsplsm_mcmc(double(squeeze(C(i,j,:,:))), bsplmodel,g.nMCMCiters);
                    
                end
            end
            
            % replace noisy connectivity estimates with
            % mean smoothed connectivity estimates
            Conn.(connmethod) = mean(PConn.(connmethod),ndims(PConn.(connmethod)));
            
            if g.verb, close(h); end
        end
        
    end

    function [X,mu,sigma] = transform(X)
        mu = 0; sigma = 1;
        
        % logit transform
        %         X = log(X./(1-X));
        
        %         [X mu sigma] = zscore(X,[],ndims(X));
        %         X = log(1000*X+1);
        X = log(X);
        %         X = 1000*log(X+10);
    end

    function X = untransform(X,mu,sigma)
        % inverse logit
        %         X = exp(X)./(1+exp(X));
        
        X = exp(X);
        
        %         X = (exp(X)-1)/1000;
        
        %         X = bsxfun(@times,X,sigma);
        %         X = bsxfun(@plus,X,mu);
    end

end



