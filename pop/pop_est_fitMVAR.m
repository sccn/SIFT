function [ALLEEG cfg] = pop_est_fitMVAR(ALLEEG,typeproc,varargin)
%
% Fit an Adaptive Autoregressive Model to EEG data. See [1] for
% mathematical details on model fitting.
%
%
% Input:
%
%   ALLEEG:         Array of EEGLAB datasets to model.
%   typeproc:       Reserved for future use. Use 0
%
% Optional:         
%
%   <'Name',value> pairs as defined in est_fitMVAR()
%   
% Output:
%
%   ALLEEG:         EEG structure(s) containing model in ALLEEG.CAT.MODEL
%   cfg:            Argument specification structure.
%
%
% See Also: est_fitMVAR()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Section 6.5.1 
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
% 
% Author: Tim Mullen 2010, SCCN/INC, UCSD. 
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
res = hlp_checkeegset(ALLEEG,{'cat'});
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
    waitbarTitle = 'Fitting VAR Models';
    
    multiWaitbar(waitbarTitle,'Reset');
    multiWaitbar(waitbarTitle,'ResetCancel',true);
    multiWaitbar(waitbarTitle,...
                 'Color', [0.8 0.0 0.1],  ...
                 'CanCancel','on',        ...
                 'CancelFcn',@(a,b)disp('[Cancel requested. Please wait...]'));
end

% execute the low-level function
for cnd=1:length(ALLEEG)
    [ALLEEG(cnd).CAT.MODEL] = feval(fcnHandle,'EEG',ALLEEG(cnd),cfg);
    
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
