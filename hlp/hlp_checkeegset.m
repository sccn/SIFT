
function res = hlp_checkeegset(EEG,checks)
%
% Check whether EEG dataset(s) contains one or more substructures.
%
% Inputs:
%
%   EEG         EEG dataset or array of datasets
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

res = {};

if nargin<2
    checks = ['cat','conn','model'];
end

if ~iscell(checks)
    checks = {checks};
end

for cnd=1:length(EEG)
    for i=1:length(checks)
        switch lower(checks{i})
            case 'cat'
                if ~isfield(EEG(cnd),'CAT') || isempty(EEG(cnd).CAT)
                    res = [res, sprintf(['SIFT:hlp_checkEEGset:\nEEG must contain CAT structure\n' ...
                        'You probably need to complete the pre-processing step first'])];
                end
            case 'conn'
                if ~isfield(EEG(cnd),'CAT') ...
                    || ~isfield(EEG(cnd).CAT,'Conn') ...
                    || isempty(EEG(cnd).CAT.Conn)
                    
                    res = [res, sprintf(['SIFT:hlp_checkEEGset:\nEEG.CAT must contain Conn structure\n' ...
                        'You need to estimate connectivity first'])];
                end
            case 'model'
                if ~isfield(EEG(cnd),'CAT') ...
                    || ~isfield(EEG(cnd).CAT,'MODEL') ...
                    || isempty(EEG(cnd).CAT.MODEL)
                    res = [res, sprintf(['SIFT:hlp_checkEEGset:\nEEG.CAT must contain MODEL structure\n' ...
                        'You need to fit a model first'])];
                end
            case 'pconn'
                if ~isfield(EEG(cnd),'CAT') ...
                    || ~isfield(EEG(cnd).CAT,'PConn') ...
                    || isempty(EEG(cnd).CAT.PConn)
                    
                    res = [res, sprintf(['SIFT:hlp_checkEEGset. EEG.CAT must contain PConn structure\n' ...
                        'You need to estimate connectivity distributions first.\n' ...
                        'See stat_* functions or type ''help stat''.'])];
                end
            case 'pnull'
                if ~isfield(EEG(cnd),'CAT') ...
                    || ~isfield(EEG(cnd).CAT,'Pnull') ...
                    || isempty(EEG(cnd).CAT.Pnull)
                    
                    res = [res, sprintf(['SIFT:hlp_checkEEGset. EEG.CAT must contain Pnull structure\n' ...
                        'You need to estimate a null distribution first.\n' ...
                        'See stat_* functions or type ''help stat''.'])];
                end
            case 'stats'
                if ~isfield(EEG(cnd),'CAT') ...
                    || ~isfield(EEG(cnd).CAT,'Stats') ...
                    || isempty(EEG(cnd).CAT.Stats)
                    
                    res = [res, sprintf(['SIFT:hlp_checkEEGset. EEG.CAT must contain Stats structure\n' ...
                        'You need to estimate statistics first'])];
                end
            case 'configs'
                if ~isfield(EEG(cnd),'CAT') ...
                    || ~isfield(EEG.CAT,'configs') ...
                    || isempty(EEG.CAT.configs)
                    res = [res, sprintf('SIFT:hlp_checkEEGset. EEG.CAT does not contain a configs structure\n')];
                end
        end
    end
end

