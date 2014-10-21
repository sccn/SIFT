function [FIT MCMC_LastState]=stat_ubsplsm(varargin)
% univariate b-spline smoothing with fPCA

% Y = cell array of [numChans x numChans x T] connectivities for each subject
% K = number of knots
% Q = number of fpca basis functions

% output:
% fit_distrib: cell array of [numChans x numChans x T x iterations]
%              posterior distribution for each subject

% Author: Tim Mullen and Wes Thompson, 2010-12, SCCN/INC, UCSD.
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

arg_define([0 1],varargin, ...
    arg_norep({'Y','TSData'},mandatory,[],sprintf(['Cell array of data to smooth.\n' ...
              'Generally, Y is a cell array where Y{i} is the [M x M x Q] matrix of M-channel connectivity for the ith subject.\n'])), ...
    arg({'smoothingLayout','MatrixElementsToSmooth'},{'diagonals','off-diagonals'},{'diagonals','off-diagonals'},'Which parts of the matrix to smooth. Diagonals (e.g. auto-connectivity) and off-diagonals (e.g. cross-connectivity) will be smoothed separately','type','logical'), ...
    arg_sub({'mcmc_opts','MCMC'},{}, ...
    {...
        arg_sub({'mcmc_init','InitMCMC'},{},@stat_ubsplsm_init,'Initialize MCMC','suppress','verb') ...
        arg_sub({'mcmc_run','RunMCMC'},{},@stat_ubsplsm_mcmc,'MCMC runtime options','suppress','verb') ...
    },'Perform Markov Chain Monte Carlo estimation'), ...
    arg_norep({'MCMC_InitState'},[],[],'Object containing initial state of Gibbs sampler. If supplied, this overrides mcmc_init'), ...
    arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );
    
% initialize some vars
init_mcmc = isempty(MCMC_InitState);    
nsubj=length(Y);
% determine the number of sources for each subject
M_i=cellfun(@(y) size(y,1),Y,'UniformOutput',true);

% determining smoothing layout
smoothDiags     = ismember_bc('diagonals',smoothingLayout);
smoothOffDiags  = ismember_bc('off-diagonals',smoothingLayout);
    
%% get smoothed time-varying connectivity coefficients

% Diagonals and off-diagonals often have different variances, so we smooth
% them separately

% get linear indices of diagonal and off-diagonal
%  elements of Y{si} conn tensor for each subject si
[dg_inds, od_inds] = deal({});
for si=1:nsubj
    mi = M_i(si);
    all_inds = 1:mi^2;
    dg_inds{si} = sub2ind([mi mi],1:mi,1:mi);
    od_inds{si} = setdiff_bc(all_inds,dg_inds{si});
end
% Note that we will "flatten" Y to cell array of form:
% CPairs = {s1(1,1) s1(1,2) s1(1,3) ... s1(N,1) s1(N,2) s1(N,3) ...
%           s2(1,1) s2(1,2) s2(1,3) ... }
% where sk(i,j) is the 1 x Q vector of connectivity values from channels 
% j to i, for the kth subject   

waitbarstr = 'Smoothing...';
% smooth diagonals
% -------------------------------------------------------------------------
if smoothDiags
    if verb==1
        fprintf(hlp_separator());
        fprintf('Smoothing diagonals...\n');
    end
    if verb==2
        multiWaitbar(waitbarstr,'Reset', ...
                     'Color',hlp_getNextUniqueColor('reset'));
    end

    % extract self-connectivity for each subject
    % in "flattened" form
    CPairs = flatten_tensor(Y,dg_inds);

    % initialize MCMC
    if init_mcmc
        MCMC_InitState_DG = hlp_microcache('ubsplsm',@stat_ubsplsm_init, ...
                                               mcmc_opts.mcmc_init, ...
                                               'T',length(CPairs{1}),'R',length(CPairs),...
                                               'verb',verb);
    else
        MCMC_InitState_DG = MCMC_InitState.diags;
    end
    
    % perform the smoothing for all subjects jointly
    [fit_dg, MCMC_LastState.diags]=stat_ubsplsm_mcmc('Y',CPairs,mcmc_opts.mcmc_run, ...
                                  'verb',verb,'MCMC_InitState',MCMC_InitState_DG);
else
    MCMC_LastState.diags = struct([]);
end


% smooth off-diagonals
% -------------------------------------------------------------------------
if smoothOffDiags
    if verb==1
        fprintf(hlp_separator());
        fprintf('Smoothing off-diagonals...\n');
    end
    if verb==2
        multiWaitbar(waitbarstr,1/3);
    end
    
    % extract cross-connectivity off-diagonal
    % elements for each subject in "flattened" form
    CPairs = flatten_tensor(Y,od_inds);

    % initialize MCMC
    if init_mcmc
        MCMC_InitState_OD = hlp_microcache('ubsplsm',@stat_ubsplsm_init, ...
                                            mcmc_opts.mcmc_init, ...
                                            'T',length(CPairs{1}),'R',length(CPairs), ...
                                            'verb',verb);
    else
        MCMC_InitState_OD = MCMC_InitState.offdiags;
    end
    
    % perform the smoothing
    [fit_od, MCMC_LastState.offdiags]=stat_ubsplsm_mcmc('Y',CPairs,mcmc_opts.mcmc_run, ...
                                     'verb',verb,'MCMC_InitState',MCMC_InitState_OD);
else
    MCMC_LastState.offdiags = struct([]);
end


% store data in [nchs x nchs x time x distribution]
if verb==1
    fprintf(hlp_separator());
    fprintf('Constructing final matrices\n');
end
if verb==2
    multiWaitbar(waitbarstr,2/3);
end
FIT = cell(1,nsubj);
ind=0;
ind_diag = 0;
for i=1:nsubj
    FIT{i} = zeros(M_i(i),M_i(i),size(fit_od,1),size(fit_od,3));
    for j1=1:M_i(i)
        for j2=1:M_i(i)
            if smoothOffDiags && (j1~=j2)
                ind=ind+1;
                FIT{i}(j1,j2,:,:)=squeeze(fit_od(:,ind,:));
            elseif smoothDiags && (j1==j2)
                ind_diag = ind_diag+1;
                FIT{i}(j1,j2,:,:)=squeeze(fit_dg(:,ind_diag,:));
            end
        end
    end
end

if verb==2
    multiWaitbar(waitbarstr,'Close');
    pause(0.1);
end

