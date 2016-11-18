function [ALLEEG cfg] = pop_est_eigenmode(ALLEEG,typeproc,varargin)
%
% Estimate eigenmodes (Principal Oscillation Patterns) of a VAR 
% process using the method described by T. Schneider and 
% A. Neumaier [1]. EEG.CAT.MODEL must be present (see est_fitMVAR). 
% Statistics available if arfit method used for model fitting. 
%
%
% Inputs:
%
%       ALLEEG:        EEG data structure with EEG.CAT.MODEL present
%                      estimated using ARFIT
%       typeproc:      Reserved for future use. Use 0
%
% Optional:         
%
%       <'Name',value> pairs as defined in est_eigenmode()       
%
% Outputs:
%
%   ALLEEG.EIGMODE
%       .modes:         Columns contain the estimated eigenmodes of the VAR model.  
%       .modeconf:      Margins of error for the components of the estimated 
%                       eigenmodes, such that (S +/- Serr) are approximate 95% 
%                       confidence intervals for the individual components of 
%                       the eigenmodes.
%       .period:        The first row contains the estimated oscillation period
%                       for each eigenmode. The second row contains margins of
%                       error such that ( period(1,k) +/- period(2,k) ) are
%                       approximate 95% confidence intervals for the period of
%                       eigenmode .modes(:,k).  For a purely relaxatory eigenmode, 
%                       the period is infinite (Inf). For an oscillatory eigenmode, 
%                       the periods are finite.
%       .dampingTime:   The first row contains the estimated oscillation
%                       damping time for each eigenmode. The second row contains 
%                       margins of error such that ( dampingTime(1,k) +/- dampingTime(2,k) )
%                       are approximate 95% confidence intervals for the damping time of
%                       eigenmode .modes(:,k)
%       .exctn:         The excitation of an eigenmode measures its dynamical importance
%                       and is returned as a fraction exctn that is normalized such that
%                       the sum of the excitations of all eigenmodes equals one.
%       .lambda:        The columns contain the eigenvalues of the
%                       eigenmodes
%   cfg:            Argument specification structure.
%
% See Also:  est_eigenmode(), est_fitMVAR()
%
% References:
%
% [1] Schneider T, Neumaier A (2001) Algorithm 808: ARfit---a matlab package
% for the estimation of parameters and eigenmodes of multivariate 
% autoregressive models. ACM Transactions on Mathematical Software 27:58-65
% http://www.gps.caltech.edu/~tapio/arfit/
% 
% Author: Tim Mullen 2010, SCCN/INC, UCSD. 
% This function uses a mod of armode() from the ARfit package by T. Schneider and 
% A. Neumaier [1]
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


if nargin<2
    typeproc = 0;
end

fcnName     = strrep(mfilename,'pop_','');
fcnHandle   = str2func(fcnName);

% check the dataset
res = hlp_checkeegset(ALLEEG,{'model'});
if ~isempty(res)
    error(['SIFT:' fcnName],res{1});
end

if isfield(ALLEEG(1).CAT.configs,fcnName)
    % get default configuration (from prior use) and merge with varargin
    varargin = [hlp_struct2varargin(ALLEEG(1).CAT.configs.(fcnName)) varargin];
end

if strcmpi(typeproc,'nogui')
    % get the config from function
    cfg = arg_tovals(arg_report('rich',fcnHandle,[{'EEG',ALLEEG(1)},varargin]),false);
else
    % render the GUI
    [PGh figh] = feval(['gui_' fcnName],ALLEEG(1),varargin{:});
    
    if isempty(PGh)
        % user chose to cancel
        cfg = [];
        return;
    end
    
    % get the specification of the PropertyGrid
    ps = PGh.GetPropertySpecification;
    cfg = arg_tovals(ps,false);
end

drawnow;

if strcmpi(typeproc,'cfg_only')
    return;
end


% initialize progress bar
if cfg.verb==2 && length(ALLEEG)>1
    waitbarTitle = 'Estimating eigenmodes';
    
    multiWaitbar(waitbarTitle,'Reset');
    multiWaitbar(waitbarTitle,'ResetCancel',true);
    multiWaitbar(waitbarTitle,...
                 'Color', [0.8 0.0 0.1],  ...
                 'CanCancel','on',        ...
                 'CancelFcn',@(a,b)disp('[Cancel requested. Please wait...]'));
end

% execute the low-level function
for cnd=1:length(ALLEEG)
    [ALLEEG(cnd).CAT.EIGMODE] = feval(fcnHandle,'EEG',ALLEEG(cnd),cfg);
    
    if ~isempty(cfg)
        % store the configuration structure
        ALLEEG(cnd).CAT.configs.(fcnName) = cfg;
    end
    
    if cfg.verb==2 && length(ALLEEG)>1
        % update waitbar
        drawnow;
        cancel = multiWaitbar(waitbarTitle,cnd/length(ALLEEG));
        if cancel && hlp_confirmWaitbarCancel(waitbarTitle)
            break;
        end
    end
    
end

% cleanup progress bar
if cfg.verb==2 && length(ALLEEG)>1
    multiWaitbar(waitbarTitle,'Close');
end
