function [ALLEEG_out cfg] = pop_sim_varmodel(typeproc,varargin)
%
% Preprocess EEG dataset(s) for connectivity analysis. See [1] for
% mathematical details on preprocessing steps.
%
%
% Input:
%
%   typeproc:       if 'nogui' don't generate GUI
%
% Optional:
%
%   <'Name',value> pairs as defined in sim_varmodel()
%
% Output:
%
%   EEG:            Simulated EEG structure(s)
%   cfg:            Argument specification structure.
%
%
% See Also: sim_varmodel()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Section 6.5.1
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% Author: Tim Mullen 2013, SCCN/INC, UCSD.
% Email:  tim@sccn.ucsd.edu
%

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

% set default output
ALLEEG_out = [];
cfg = [];
        
fcnName     = strrep(mfilename,'pop_','');
fcnHandle   = str2func(fcnName);

if strcmpi(typeproc,'nogui')
    % get the default config from function and overload supplied args
    cfg = arg_tovals(arg_report('rich',fcnHandle,varargin),false);
else
    % render the GUI
    [PGh figh] = feval(['gui_' fcnName],varargin{:});
    
    if isempty(PGh)
        % user chose to cancel
        return;
    end
    
    % get the specification of the PropertyGrid
    ps = PGh.GetPropertySpecification;
    cfg = arg_tovals(ps,false);
end

drawnow;

if ~cfg.makeEEGset.arg_selection
    error('SIFT:sim_varmodel',['If using pop_' fcnName '(), you must enable the BuildEEGLABStructure option.\n' ...
                               'Use ' fcnName '() from the command-line to return a raw dataset']);
end

if strcmpi(typeproc,'cfg_only')
    return;
end

% execute the low-level function
[ALLEEG_out] = feval(fcnHandle,cfg);
    
if ~isempty(cfg)
    % store the configuration structure
    ALLEEG_out.CAT.configs.(fcnName) = cfg;
end


