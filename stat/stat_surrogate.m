function [PConn g] = stat_surrogate(varargin)
%
% Generate surrogate statistical distributions of SIFT connectivity and
% other estimators. This function can estimate bootstrap, jacknife, k-fold crossvalidation and
% phase-randomized (null) distributions.
%
%
% ===============================================
%    This function is under development and may
%    be unstable.
%    Please check sccn.uscd.edu/wiki/SIFT for 
%    updated version
% ===============================================
%
% === TODO ============================================
% - add memory check option
% - retest phaserand
% - add warning when attempting single-trial bootstrap
% =====================================================
%
%
% Input                    Information
% ------------------------------------------------------------------------------------------------------------------------------
% ALLEEG:                  EEG data structure containing connectivity structure in EEG.CAT.Conn
%
% configs:                 Config superstructure. Contains configuration structures as fields: prepcfg,modfitcfg,conncfg. 
%                          These are returned from the respective pre_prepData(), est_fitMVAR(), and est_mvarConnectivity() 
%                          functions
%
% Optional                 Information                                                                                           
% -------------------------------------------------------------------------------------------------------------------------------
% Mode:                    Resampling modes                                                                                      
%                          Bootstrap (Efron Bootstrap resampling with replacement), Jacknife (leave-one-out cross-validation),   
%                          Crossval (k-fold cross-validation), PhaseRand (Theiler phase randomization)                                                                    
%                          Possible values: 'Bootstrap','Jacknife','InverseJacknife','Crossval','PhaseRand'                      
%                          Default value  : 'Bootstrap'                                                                          
%                          Input Data Type: string                                                                               
%                                                                                                                                
%     | NumFolds:          Number of folds                                                                                       
%                          This performs k-fold resampling                                                                       
%                          Input Range  : Unrestricted                                                                           
%                          Default value: 10                                                                                     
%                          Input Data Type: real number (double)                                                                 
%                                                                                                                                
%     | NumPermutations:   Number of resamples                                                                                   
%                          This performs Theiler phase randomization                                                             
%                          Input Range  : Unrestricted                                                                           
%                          Default value: 200                                                                                    
%                          Input Data Type: real number (double)                                                                 
%                                                                                                                                
%     | NumPermutations:   Number of resamples                                                                                   
%                          This performs efron bootstrapping                                                                     
%                          Input Range  : Unrestricted                                                                           
%                          Default value: 200                                                                                    
%                          Input Data Type: real number (double)                                                                 
%                                                                                                                                
% AutoSave:                Autosave distributions                                                                                
%                          This will periodically save the computation to a mat file. In case of system crash, bootstrapping     
%                          can be resumed from this file                                                                         
%                          Input Range  : Unrestricted                                                                           
%                          Default value: 0                                                                                      
%                          Input Data Type: boolean                                                                              
%                                                                                                                                
%     | FileNamePrefix:    Prefix (including optional path) for autosave file                                                    
%                          Data is saved in <FileNamePrefix>.part                                                                
%                          Possible values: Unrestricted                                                                         
%                          Default value  : 'SIFT_boostrap'                                                                      
%                          Input Data Type: string                                                                               
%                                                                                                                                
%     | AutoSaveFrequency: Fractional increment between saves                                                                    
%                          For example, 0.25 = quarterly saves                                                                   
%                          Input Range  : [2.2204e-16           1]                                                               
%                          Default value: 0.25                                                                                   
%                          Input Data Type: real number (double)                                                                 
%                                                                                                                                
% ConnectivityMethods:     Connectivity estimator(s) to bootstrap                                                                
%                          All connectivity estimators must be named exactly as they appear in EEG.CAT.Conn                      
%                          Possible values: ''                                                                                   
%                          Default value  : 'n/a'                                                                                
%                          Input Data Type: boolean                                                                              
%                                                                                                                                
% VerbosityLevel:          Verbosity level. 0 = no output, 1 = text, 2 = graphical                                               
%                          Possible values: 0,1,2                                                                                
%                          Default value  : 2                                                                                    
%                          Input Data Type: real number (double) 
%
% Output                 Information                                                                                           
% -------------------------------------------------------------------------------------------------------------------------------
% PConn:                 Surrogate connectivity structure. Same format as EEG.CAT.Conn but with resamples stored in the last 
%                        dimension of each connectivity matrix (e.g. [M x M x Time x Freq x Resamples])   
%
% See Also: stat_analyticStats(), statcond(), stat_surrogateStats()
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
    error('surrogate stats currently available only for individual datasets');
end

if ~isempty(ALLEEG)
    if isfield(ALLEEG,'CAT') && isfield(ALLEEG.CAT,'Conn') && ~isempty(ALLEEG.CAT.Conn)
        Conn = ALLEEG.CAT.Conn(1);
        ConnNames   = hlp_getConnMethodNames(Conn);
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

g = arg_define([0 Inf],varargin, ...
    arg_norep('ALLEEG',mandatory,[],'EEG structure.'), ...
    arg_norep('configs',mandatory,[],'Config superstructure. Contains configuration structures as fields: prepcfg,modfitcfg,conncfg. These are returned from the respective pre_prepData(), est_fitMVAR(), and est_mvarConnectivity() functions'), ...
    arg_subswitch({'mode','Mode'},'Bootstrap', ...
    { ...
    'Bootstrap' { ...
    arg({'nperms','NumPermutations'},200,[],'Number of resamples. This performs efron bootstrapping'), ...
    }, ...
    'Jacknife', { ...
    arg_norep('dummy',[],[]), ...
    }, ...
    'InverseJacknife', { ...
    arg_norep('dummy2',[],[]), ...
    }, ...
    'Crossval' { ...
    arg({'nfolds','NumFolds'},10,[],'Number of folds. This performs k-fold resampling'), ...
    }, ...
    'PhaseRand', { ...
    arg({'nperms','NumPermutations'},200,[],'Number of resamples. This performs Theiler phase randomization'), ...
    }, ...
    }, 'Resampling modes. Bootstrap (Efron Bootstrap resampling with replacement), Jacknife (leave-one-out cross-validation), Crossval (k-fold cross-validation), PhaseRand (Theiler phase randomization)'), ...
    arg_subtoggle({'autosave','AutoSave'},'off', ...
    { ...
    arg({'savefname','FileNamePrefix'},'SIFT_boostrap','','Prefix (including optional path) for autosave file. Data is saved in <FileNamePrefix>.part'), ...
    arg({'savefrq','AutoSaveFrequency'},0.25,[eps 1-eps],'Fractional increment between saves. For example, 0.25 = quarterly saves'), ...
    }, 'Autosave distributions. This will periodically save the computation to a mat file. In case of system crash, bootstrapping can be resumed from this file'), ...
    arg({'connmethods','ConnectivityMethods'},conndef,ConnNames,'Connectivity estimator(s) to bootstrap. All connectivity estimators must be named exactly as they appear in EEG.CAT.Conn','type','logical'), ...
    arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );


arg_toworkspace(g);

if (ALLEEG.trials==1 || isempty(ALLEEG.trials)) ...
        && any(ismember(lower(mode.arg_selection),{'bootstrap','jacknife','inversejacknife','crossval'}))
    error(['Unable to compute bootstrap distributions for a single trial. ' char(10) ...
           'You must use either Analytic Statistics or Phase Randomization surrogate option']);
end
    
% Perform boostrapping
switch lower(mode.arg_selection)
    case 'bootstrap'
        nperms = mode.nperms;
        cv.training = @(x)(randsample(ALLEEG.trials, ALLEEG.trials, true));
    case 'jacknife'
        nperms = ALLEEG.trials;
        cv = cvpartition(ALLEEG.trials,'leaveout');
    case 'inversejacknife'
        nperms = ALLEEG.trials;
        cv = cvpartition(ALLEEG.trials,'leaveout');
    case 'crossval'
        nperms = mode.nfolds;
        cv = cvpartition(ALLEEG.trials,'kfold',nperms);
    case 'phaserand'
        nperms = mode.nperms;
end


% Initiate progress bar
if verb==1
    progress('init',sprintf('Resampling (%d/%d)...',0,nperms));
    tic
    t0 = toc;
elseif verb==2
    h=statusbar(sprintf('Resampling (x%d)...',nperms));
end

if autosave.arg_selection
    savefrq = round(autosave.savefrq*nperms);
end

for perm=1:nperms
    
    if strcmpi(mode.arg_selection,'inversejacknife') % (fit model to a random single-trial)
        % preprocess data
        EEG = pre_prepData('ALLEEG',ALLEEG,configs.prepData,'verb',0,'newtrials',cv.test(perm));
        
    elseif ~strcmpi(mode.arg_selection,'phaserand')
        % all bootstrapping/resampling modes
        
        % preprocess data
        EEG = pre_prepData('ALLEEG',ALLEEG,configs.prepData,'verb',1,'newtrials',cv.training(perm));
        
    else % (Phase randomization)
        
        % preprocess data (only once)
        if perm==1
            ALLEEG = pre_prepData('ALLEEG',ALLEEG,configs.prepData,'verb',0);
        end
        
        EEG = ALLEEG;
        
        % Create surrogate data by randomizing phases.
        % Multiply each fourier amplitude by e^{iw)
        % where w is a random phase chosen in [0 2pi]
        % (c.f. Theiler, et al 1997)
        % to generate hermitian phase distributions, I extract
        % the phase of a random matrix. This ensures the random
        % spectrum is conjugate symmetric
        EEG.CAT.srcdata = permute(EEG.CAT.srcdata,[2 1 3]);
        [npnts nchs ntr] = size(EEG.CAT.srcdata);
        
        % NFFT = 2^nextpow2(npnts); size(data,1);
        
        for tr=1:ntr
            EEG.CAT.srcdata(:,:,tr) = ...
                ifft(abs(fft(EEG.CAT.srcdata(:,:,tr))) ...
                .* exp(1i*angle(fft(rand(npnts,nchs)))), ...
                'symmetric');
        end
        
        EEG.CAT.srcdata = ipermute(EEG.CAT.srcdata,[2 1 3]);
        
    end
    
    
    % fit model
    EEG = pop_est_fitMVAR(EEG,'nogui',configs.fitMVAR,'verb',0);
    
    % calculate connectivity
    EEG = pop_est_mvarConnectivity(EEG,configs.mvarConnectivity,'verb',0);
    
    % append causal estimates to the last dimension of respective arrays in Conn
    if isempty(connmethods)
        connfields = hlp_getConnMethodNames(EEG.CAT.Conn);
    else
        connfields = connmethods;
    end
    
    if perm==1
        for m=1:length(connfields)
            if ~isempty(EEG.CAT.Conn.(connfields{m}))
                % initialize matrices
                % permutations will be temporarily stored in first
                % dimension
                PConn.(connfields{m}) = zeros([nperms size(EEG.CAT.Conn.(connfields{m}))],'single');
            end
        end
    end
    
    for m=1:length(connfields)
        if ~isempty(EEG.CAT.Conn.(connfields{m}))
            % store permutation
            PConn.(connfields{m})(perm,:,:,:,:,:,:,:) = single(EEG.CAT.Conn.(connfields{m}));
        end
    end
    
    % autosave checkpoint
    if autosave.arg_selection && ~isempty(autosave.savefname) && ~mod(perm,savefrq)
        save(sprintf('%s.part',autosave.savefname),'PConn','-v7.3');
    end
    
    % update progress bar
    if verb==1
        te = toc-t0;
        tt = ceil(te*(nperms-perm)/perm);
        progress(perm/nperms,sprintf('Resampling [%d/%d] (EL: %0.1fm / ETA: %0.1fm)',perm,nperms,ceil(te)/60,tt/60));
    elseif verb==2
        statusbar(perm/nperms,h);
    end
    
end

for m=1:length(connfields)
    sz = size(PConn.(connfields{m}));
    % put permutations in last dimension
    PConn.(connfields{m}) = permute(PConn.(connfields{m}),circshift((1:length(sz))',-1));
end
    
PConn.winCenterTimes = EEG.CAT.Conn.winCenterTimes;
PConn.erWinCenterTimes = EEG.CAT.Conn.erWinCenterTimes;
PConn.freqs = EEG.CAT.Conn.freqs;
PConn.mode = mode.arg_selection;

% clean up
if verb==1
    fprintf('\n')
elseif verb==2
    delete(h);
end

if autosave.arg_selection && ~isempty(autosave.savefname),
    if exist(sprintf('%s.part',autosave.savefname),'file')
        delete(sprintf('%s.part',autosave.savefname));
    end
end
