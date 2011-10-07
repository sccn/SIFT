function [ALLEEG cfg] = pop_stat_surrogate(ALLEEG,typeproc,varargin)
%
% Return surrogate statistics based on a surrogate distribution (bootstrap, jacknife, etc).
% Each surrogate distribution should approximate the distribution of the
% estimator. 
%
%
% Input:
%
%   ALLEEG:         EEGLAB dataset to preprocess.
%   typeproc:       Reserved for future use. Use 0
%
% Optional:         
%
%   <'Name',value> pairs as defined in stat_surrogate()
%   
% Output:
%
%   ALLEEG:         EEG structure(s) with Stats object stored in
%                   ALLEEG.CAT.Stats
%   cfg:            Argument specification structure.
%
%
% See Also: stat_surrogate()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
% 
% Author: Tim Mullen 2010-2011, SCCN/INC, UCSD. 
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


if isunix
    SLASH = '/';
else
    SLASH = '\';
end

if (~isfield(ALLEEG(1).CAT,'Conn')), 
    errordlg2('Please compute connectivity first!','Surrogate Statistics');
end
    
% [fnpath fnname] = fileparts(which('pop_stat_surrogate'));
% if isempty(varargin)
%     if exist('preprep.cfg','file')
%         load('preprep.cfg','-mat');
%     else
%         cfg = [];
%     end
%     varargin = {cfg};
% end


% render the GUI
[PGh figh] = gui_surrogate(ALLEEG(1),varargin{:});

if isempty(PGh)
    % user chose to cancel
    cfg = [];

    return;
end

% get the specification of the PropertyGrid
ps = PGh.GetPropertySpecification;
cfg = arg_tovals(ps,false);

drawnow
    
% save([fnpath SLASH '@configs' SLASH 'preprep.cfg'],'cfg');

for cnd=1:length(ALLEEG)
    res=questdlg2(['In the following dialog, please load the original (not processed by SIFT) dataset corresponding to Dataset '...
             fastif(isempty(ALLEEG(cnd).filename),num2str(cnd),ALLEEG(cnd).filename) '.' char(10) ...
             'Press Cancel to use the currently selected EEG datasets.'], ...
             'Loading data for Surrogate Statistics...','OK','Cancel','OK');
         switch lower(res)
             case 'ok'
                tmp = pop_loadset;
                
                if isempty(tmp)
                    return;
                else
                    ALLEEGorig(cnd) = tmp;
                end
             case 'cancel'
                 ALLEEGorig(cnd) = ALLEEG(cnd);
         end
end

% execute the low-level function
for cnd=1:length(ALLEEG)
    [ALLEEG(cnd).CAT.PConn] = stat_surrogate('ALLEEG',ALLEEGorig(cnd),'configs',ALLEEG(cnd).CAT.configs,cfg);
end

