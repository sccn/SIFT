% std_sift() - Compute model for SIFT
%
% Usage:    
% >> CAT = std_sift(EEG, options);
%
% Required inputs:
%   EEG          - EEG dataset
%
% Optional inputs:
%  options
% 
% Outputs:
%   CAT - CAT structure of SIFT
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, 2013-

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2013, arno@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
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

function CAT = std_sift(EEG, varargin)
    
    if nargin < 2
        help std_sift;
        return;
    end;
    
    EEG = sift_xxxxxx(EEG, varargin{:});
    CAT = EEG.CAT;
    