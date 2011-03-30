
function res = hlp_checkeegset(EEG,checks)
%
% Check whether EEG dataset contains one or more substructures.
%
% Inputs:
% 
%   EEG         EEG dataset
%   checks      Cell vector containing one or more checks to perform
%               Possible checks: {'cat','conn','model'}.
% Outputs:
%
%   res         cell array containing results of checks where res{i} is the 
%               result of checking for checks{i}
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift/
% 
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
        checks = ['cat','conn','model'];
    end
    
    res = {};
    
    for cnd=1:length(EEG)
        for i=1:length(checks)
            switch lower(checks{i})
                case 'cat'
                    try
                        EEG(cnd).CAT;
                    catch
                        res = [res, 'error:hlp_checkEEGset. EEG must contain CAT structure'];
                    end
                case 'conn'
                    try
                        EEG(cnd).CAT.Conn;
                    catch
                        res = [res, 'error:hlp_checkEEGset. EEG.CAT must contain Conn structure'];
                    end
                case 'model'
                    try
                        EEG(cnd).CAT.MODEL;
                    catch
                        res = [res, 'error:hlp_checkEEGset. EEG.CAT must contain MODEL structure'];
                    end
            end
        end
    end
    
    