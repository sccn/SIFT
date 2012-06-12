
function varargout = pop_est_mvarConnectivity(ALLEEG,typeproc,varargin)
%
% Computes connectivity estimates from a precomputed MODEL 
% and stores the result in ALLEEG.CAT.Conn
% 
% Input:
%
%   ALLEEG          EEG structure with fields CAT.MODEL
%
% Optional:
%
%   'connmethods':  cell array of strings denoting connectivity methods to
%                   compute (parenthesized acronym from list below)
% 
%                     DIRECTED TRANSFER FUNCTION MEASURES:
%                         Directed Tranfer Function (DTF)
%                         Normalized DTF (nDTF)
%                         Direct DTF (dDTF)
%                         Direct DTF (with full causal normalization) (dDTF08)
%                         Full-frequency DTF (ffDTF)
%                      PARTIAL DIRECTED COHERENCE MEASURES
%                         Partial Directed Coherence (PDC)
%                         Normalized PDC (nPDC)
%                         Generalized Partial Directed Coherence (GPDC)
%                         Partial Directed Coherence Factor (PDCF)
%                         Renormalized Partial Directed Coherence (RPDC)
%                      GRANGER-GEWEKE CAUSALITY MEASURES
%                         Granger-Geweke Causality (GGC)
%                      SPECTRAL COHERENCE MEASURES
%                         Complex Coherence (Coh)
%                         Imaginary Coherence (iCoh)
%                         Partial Coherence (pCoh)
%                         Multiple Coherence (mCoh)
%                      SPECTRAL DENSITY MEASURES
%                         Complex Spectral Density (S)
%   'absvalsq':     Boolean (true,false) determining whether to return the
%                   square of the absolute value of complex measures (def:
%                   true)
%   'spectraldecibels': Boolean (true,false) determining whether to return
%                       the spectral power in units of decibels
%                       (10*log10(Power)) (def: false)
%
% Output:
%
%   ALLEEG          EEG structure with results in EEG.CAT.Conn.
%                   Conn.(connmethod) is a [num_chans x num_chans x
%                   num_freqs x num_time] connectivity matrix
%   params          The options used in the connectivity estimation
%   
% See Also: est_mvarConnectivity(), pop_est_fitMVAR()
%
% References: 
% 
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
% 
% Author: Tim Mullen, 2010, SCCN/INC, UCSD. 
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

res = hlp_checkeegset(ALLEEG,{'cat','model'});
if ~isempty(res)
    error(res{end});
end

if isempty(ALLEEG(1).CAT.MODEL)
    error('SIFT:pop_est_mvarConnectivity','You must fit an MVAR model first');
end

if nargin<2
    typeproc = 0;
end

if isfield(ALLEEG(1).CAT.configs,'mvarConnectivity')
    % get default configuration (from prior use) and merge with varargin
    varargin = [hlp_struct2varargin(ALLEEG(1).CAT.configs.mvarConnectivity) varargin];
end

if strcmpi(typeproc,'nogui')
    % get the config from function
    cfg = arg_tovals(arg_report('rich',@est_mvarConnectivity,[{'EEG',ALLEEG(1) 'MODEL', ALLEEG(1).CAT.MODEL},varargin]));
else
    % render the GUI
    [PGh figh] = gui_est_mvarConnectivity(ALLEEG,varargin{:});

    if isempty(PGh)
        % user chose to cancel
        return;
    end

    % get the specification of the PropertyGrid
    ps = PGh.GetPropertySpecification;
    cfg = arg_tovals(ps,false);
end


% now calculate connectivity
for cond=1:length(ALLEEG)
    ALLEEG(cond).CAT.Conn = est_mvarConnectivity(ALLEEG(cond),ALLEEG(cond).CAT.MODEL,cfg);

    % clear any existing visualization GUI config files
    try, ALLEEG(cond).CAT.configs.TimeFreqGrid = []; catch, end;
    try, ALLEEG(cond).CAT.configs.BrainMovie3D = []; catch, end;
    
    ALLEEG(cond).CAT.configs.mvarConnectivity = cfg;
end

varargout{1} = ALLEEG;
varargout{2} = cfg;


