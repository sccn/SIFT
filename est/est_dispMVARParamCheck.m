function checkcode = est_dispMVARParamCheck(ALLEEG,g)
% display results of sanity checks on MVAR parameters
%
% Input:
%
%   ALLEEG             EEG data structure.
%   g                  struct of MVAR parameters as returned from
%                      est_fitMVAR()
%
%
% Output:
%
%   checkcode:         'ok'      - all checks passed
%                      'error'   - parameters result in critical error
%                      'warning' - possible problem with chosen params
%
% See Also: est_checkMVARParams(), est_fitMVAR()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
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


checkcode = 'ok';

% check that parameters are OK
try
    [infostring warnstring errstring] = est_checkMVARParams(ALLEEG,g);
catch err
    %         errordlg2(err.message,'Error in MVAR configuration!');
    checkcode = err;
    return;
end


if g.verb>0
    % display summary on command line
    for str=1:length(infostring)
        if ~isempty(errstring{str})
            fprintf('ERROR\t');
        elseif ~isempty(warnstring{str})
            fprintf('WARNING\t');
        else
            fprintf('OK\t');
        end
        fprintf(infostring{str});
        fprintf(warnstring{str});
        fprintf(errstring{str});
        fprintf('\n');
    end
end

if ~all(ismember(errstring,''))
    checkcode = 'error';
elseif ~all(ismember(warnstring,''))
    checkcode = 'warning';
end