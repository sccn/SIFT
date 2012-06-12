function [ALLEEG cfg] = pop_est_fitMVAR(ALLEEG,typeproc,varargin)
%
% Preprocess EEG dataset(s) for connectivity analysis. See [1] for
% mathematical details on preprocessing steps.
%
%
% Input:
%
%   ALLEEG:         Array of EEGLAB datasets to preprocess.
%   typeproc:       Reserved for future use. Use 0
%
% Optional:         
%
%   <'Name',value> pairs as defined in pre_prepData()
%   
% Output:
%
%   ALLEEG:         Prepocessed EEG structure(s)
%   cfg:            Argument specification structure.
%
%
% See Also: pre_prepData()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Section 6.5.1 
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
% 
% Author: Tim Mullen 2009, SCCN/INC, UCSD. 
% Email:  tim@sccn.ucsd.edu
% 
% Revised Jan 2010.

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

% check the dataset
res = hlp_checkeegset(ALLEEG,{'cat'});
if ~isempty(res)
    error('SIFT:est_fitMVAR',res{1});
end

if isfield(ALLEEG(1).CAT.configs,'fitMVAR')
    % get default configuration (from prior use) and merge with varargin
    varargin = [hlp_struct2varargin(ALLEEG(1).CAT.configs.fitMVAR) varargin];
end

if strcmpi(typeproc,'nogui')
    % get the config from function
    cfg = arg_tovals(arg_report('rich',@est_fitMVAR,[{'EEG',ALLEEG(1)},varargin]));
else
    % render the GUI
    [PGh figh] = gui_est_fitMVAR(ALLEEG,varargin{:});

    if isempty(PGh)
        % user chose to cancel
        cfg = [];
        return;
    end

    % get the specification of the PropertyGrid
    ps = PGh.GetPropertySpecification;
    cfg = arg_tovals(ps,false);
end

% save([fnpath SLASH '@configs' SLASH 'preprep.cfg'],'cfg');

% fit the MVAR model
for cond=1:length(ALLEEG)
    [ALLEEG(cond).CAT.MODEL] = est_fitMVAR('EEG',ALLEEG(cond),cfg);
end

% save the configuration file
for cond=1:length(ALLEEG)
    ALLEEG(cond).CAT.configs.fitMVAR = cfg;
end

