
function varargout = pop_est_mvarConnectivity(ALLEEG,varargin)
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


varargout{1} = ALLEEG;

if nargin < 2
   popup = 1; 
else
    var = hlp_mergeVarargin(varargin{:});
    myargs = {'connmethods'  'cell'    {'DTF','dDTF','dDTF08','ffDTF','nDTF','GGC', 'iCoh','Coh','S','pCoh','pCoh2','mCoh','RPDC','GPDC','nPDC','PDCF','PDC'}     {'DTF'}; ...
              'absvalsq'     'boolean' []   true; ...
              'spectraldecibels'    'boolean'   []      false; ...
              'freqs',              'real'      []     (1 : fix(ALLEEG(1).srate/2)-1); ...
              'verb',               'boolean'   [0 2]  2};
    g = finputcheck(var, [myargs; hlp_getDefaultArglist('est')], 'pop_est_mvarConnectivity','ignore','quiet');
    if ischar(g), error(g); end
    g.connmethods = unique(g.connmethods);
    if nargout > 1, params = g; end
    popup=0;
end

res = hlp_checkeegset(ALLEEG,{'cat','model'});
if ~isempty(res)
    error(res{end});
end

g.winlen = ALLEEG(1).CAT.MODEL.winlen;
g.winstep = ALLEEG(1).CAT.MODEL.winstep;
connmethods = {'','DTF','nDTF','dDTF','dDTF08','ffDTF','','PDC','nPDC','GPDC','PDCF','RPDC','','GGC','','Coh','iCoh','pCoh','mCoh','','S'};
connmethodsFullNames = {'+ DIRECTED TRANSFER FUNCTION MEASURES', ...
                        '     Directed Tranfer Function (DTF)',...
                        '     Normalized DTF (nDTF)',...
                        '     Direct DTF (dDTF)',...
                        '     Direct DTF (with full causal normalization)',...
                        '     Full-frequency DTF (ffDTF)',...
                        '+ PARTIAL DIRECTED COHERENCE MEASURES', ...
                        '     Partial Directed Coherence (PDC)',...
                        '     Normalized PDC (nPDC)', ...
                        '     Generalized Partial Directed Coherence (GPDC)', ...
                        '     Partial Directed Coherence Factor (PDCF)',...
                        '     Renormalized Partial Directed Coherence (RPDC)', ...
                        '+ GRANGER-GEWEKE CAUSALITY MEASURES', ...
                        '     Granger-Geweke Causality (GGC)', ...
                        '+ SPECTRAL COHERENCE MEASURES',...
                        '     Complex Coherence (Coh)', ...
                        '     Imaginary Coherence (iCoh)', ...
                        '     Partial Coherence (pCoh)', ...
                        '     Multiple Coherence (mCoh)', ...
                        '+ SPECTRAL DENSITY MEASURES', ...
                        '     Complex Spectral Density', ...
                        };

if popup
    
    geomhoriz = {1 1 1 1 1 [1 2]};
    uilist = { ...
               { 'Style', 'text', 'string', 'Select connectivity measures to calculate' }...
               { 'Style', 'text', 'string', '(hold Ctrl to select multiple)' }...
               { 'Style', 'listbox', 'string', connmethodsFullNames, 'tag', 'lstConnmethods','Value',2,'Min',1,'Max',20} ...
               { 'Style', 'checkbox', 'string', 'return squared amplitude of complex measures', 'value', true, 'tag','chkSquareModulus'}...
               { 'Style', 'checkbox', 'string', 'convert spectral density to decibels', 'value', false, 'tag','chkSpectralDecibels'}...
               { 'Style', 'text', 'string', 'Frequencies (Hz)'} ...
               { 'Style', 'edit', 'string', ['1: ' num2str(fix(ALLEEG(1).srate/2)-1)], 'tag','edtFreqs'}...
			 };

	[ tmp1 tmp2 strhalt result ] = inputgui( 'geometry', geomhoriz, 'geomvert',[1 1 10 1 1 1], 'uilist',uilist, 'helpcom','pophelp(''pop_est_mvarConnectivity'');', ...
					   'title','Calculate Connectivity Measures');
	if isempty( tmp1 ), return; end;
    
    if isempty(result.edtFreqs)
        errordlg2('You must specify a list of frequencies to calculate','Connectivity Estimation');
        return;
    end
    
    g.freqs = str2num(result.edtFreqs);
    g.connmethods = strtrim(connmethods(result.lstConnmethods));
    g.connmethods(cellfun(@(x)(isempty(x)),g.connmethods))=[];
    g.absvalsq = result.chkSquareModulus;
    g.spectraldecibels = result.chkSpectralDecibels;
    g.verb = 2;
end


% now calculate connectivity
for cond=1:length(ALLEEG)
    ALLEEG(cond).CAT.Conn = est_mvarConnectivity(ALLEEG(cond),ALLEEG(cond).CAT.MODEL,g);
    
    if g.absvalsq
        ALLEEG(cond).CAT.Conn = hlp_absvalsq(ALLEEG(cond).CAT.Conn,g.connmethods,false,g.verb);
    end
    
    if g.spectraldecibels && isfield(ALLEEG(cond).CAT.Conn,'S')
        ALLEEG(cond).CAT.Conn.S = 10*log10(ALLEEG(cond).CAT.Conn.S);
    end
    
    % clear any existing visualization GUI config files
    try, ALLEEG(cond).CAT.configs.TimeFreqGrid = []; catch, end;
    try, ALLEEG(cond).CAT.configs.BrainMovie3D = []; catch, end;
    
    ALLEEG(cond).CAT.configs.mvarConnectivity = g;
end

varargout{1} = ALLEEG;
varargout{2} = g;


