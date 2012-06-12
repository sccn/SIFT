function [AA] = hlp_aamp(varargin)

% Obtain the instantaneous analytic amplitude within a given band of a collection of
% processes in EEG.CAT.srcdata. This uses the hilbert transform and EEGLAB's eegfilt().
%
% Input:
%
%   EEG                 Preprocessed EEG structure.
%   AmplitudePassBand   The [lo hi] pass-band for which to estimate instantaneous amplitude
%
% Outputs
%   AA                  analytic amplitude of same dimensions as
%                       EEG.CAT.srcdata
%
% See Also: eegfilt(), hilbert(), est_PMGC()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/SIFT/
% 
% 
% Author: Tim Mullen 2012, SCCN/INC, UCSD. 
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
        arg_norep({'EEG'},mandatory),...
        arg({'ampband','AmplitudePassBand'},[80 150],[],'The [lo hi] pass-band to use for amplitude (Hz)','shape','row'), ...
        arg({'verb','Verbosity'},true,[],'Verbose output') ...
    );


% obtain the analytic amplitude of each time-series
[nchs npnts ntr] = size(g.EEG.CAT.srcdata);

if g.verb,
    h = waitbar(0,'Computing analytic amplitude using Hilbert transform...');
end


AA = zeros(nchs, npnts, ntr);
filtd = eegfilt(g.EEG.CAT.srcdata(:,:),g.EEG.srate,g.ampband(1),[],npnts);
filtd = eegfilt(filtd,g.EEG.srate,[],g.ampband(2),npnts);
filtd = reshape(filtd,[nchs, npnts, ntr]);

if nchs>ntr
    % hilbert transform each trial separately (all channels at once)
    for tr=1:ntr
        if g.verb,
            waitbar(tr/ntr,h,sprintf('Computing analytic amplitude using Hilbert transform (%d/%d)...',tr,ntr));
        end
        AA(:,:,tr)=abs(hilbert(filtd(:,:,tr)'))'; 
    end
else
    % hilbert transform each channel separately (all trials at once)
    for ch=1:nchs
        if g.verb,
            waitbar(ch/nchs,h,sprintf('Computing analytic amplitude using Hilbert transform (%d/%d)...',ch,nchs));
        end
        AA(ch,:,:) = abs(hilbert(filtd(ch,:,:))); 
    end
end

if g.verb
    close(h);
end
