function [ALLEEG cfg] = pop_pre_prepData(ALLEEG,typeproc,varargin)
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



eeg_options;
if ~option_computeica
    fprintf('Enabling "computeica" option in EEGLAB prefs\n');
    pop_editoptions( 'option_computeica', 1);
%     error('Please enable the "computeica" option in eeg_options');
end

for cond = 1:length(ALLEEG)
    if ~isfield(ALLEEG(cond),'CAT') || ~isfield(ALLEEG(cond).CAT,'curComps')
        ALLEEG(cond).CAT.curComps = 1:size(ALLEEG(cond).icaweights,1);
    end
    if  ~isfield(ALLEEG(cond),'CAT') || ~isfield(ALLEEG(cond).CAT,'MODEL')
        ALLEEG(cond).CAT.MODEL = [];
    end
end

splashscreen;

if isunix
    SLASH = '/';
else
    SLASH = '\';
end

% [fnpath fnname] = fileparts(which('pop_pre_prepData'));
% if isempty(varargin)
%     if exist('preprep.cfg','file')
%         load('preprep.cfg','-mat');
%     else
%         cfg = [];
%     end
%     varargin = {cfg};
% end

% render the GUI
[PGh figh] = gui_prepData(ALLEEG,varargin{:});

if isempty(PGh)
    % user chose to cancel
    cfg = [];
    return;
end

% get the specification of the PropertyGrid
ps = PGh.GetPropertySpecification;
cfg = arg_tovals(ps,false);

% save([fnpath SLASH '@configs' SLASH 'preprep.cfg'],'cfg');

% execute the low-level function
[ALLEEG cfg] = pre_prepData('ALLEEG',ALLEEG,cfg);
