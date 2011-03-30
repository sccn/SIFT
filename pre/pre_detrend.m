function [EEG g] = pre_detrend(varargin)
% 
% Detrend each channel and trial using a least-squares linear fit.
% Alternately, one can center the data by removing the mean from each
% trial. See [1] for more details.
%
% Inputs:
%
%   EEG:        EEG data structure
%
% Optional:     <'Name',Value> pairs
%
%     VerbosityLevel:   Verbosity level. 0 = no output, 1 = text, 2 = graphical                                           
%                       Possible values: 0,1,2                                                                            
%                       Default value  : 1                                                                                
%                       Input Data Type: real number (double)                                                             
% 
%     DetrendingMethod: Detrending options                                                                                
%                       Linear: removes the least-squares fit of a straight line from each trial. 
%                       Constant: removes the mean from each trial (centering)                                                                  
%                       Possible values: 'linear','constant'                                                              
%                       Default value  : 'linear'                                                                         
%                       Input Data Type: boolean  
% Outputs:
%
%   EEG:        processed EEG structure
%   g:          argument specification structure
%
% See Also: pop_pre_prepData(), pre_prepData(), detrend()
%
% Refences:
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

g = arg_define([0 1], varargin, ...
        arg_norep('EEG',mandatory), ...
        arg_nogui({'verb','VerbosityLevel'},1,{int32(0) int32(1) int32(2)}, ...
                   'Verbosity level. 0 = no output, 1 = text, 2 = graphical'), ...
        arg({'method','DetrendingMethod'},{'linear'},{'linear','constant'}, ...
             'Detrending options. Linear: removes the least-squares fit of a straight line from each trial. Constant: removes the mean from each trial (centering)','type','logical')); 

% commit ALLEEG variable to workspace
[data g] = hlp_splitstruct(g,{'EEG'});
arg_toworkspace(data);
clear data;

for i=1:length(g.method)
    m = g.method{i};
    if g.verb && strcmpi(m,'mean'), disp('Centering data...'); end
    if g.verb && strcmpi(m,'linear'), disp('Detrending data...'); end
    
    % apply the detrending/centering
    for ch=1:size(EEG.icaact,1)
        EEG.data(ch,:,:) = detrend(squeeze(EEG.data(ch,:,:)),m);
    end
end

if g.verb, disp('Done!'); end


