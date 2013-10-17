
function Conn = est_cfc_pacmi(varargin)
% Calculate the Modulation Index measure for event-related phase-amplitude 
% coupling. Based on code by Ryan Canolty [1].
%
% References: 
%
% [1] Canolty RT, Edwards E, Dalal SS, Soltani nchs, Nagarajan SS, Kirsch HE, 
%     Berger MS, Barbaro NM, Knight RT. High gamma power is phase-locked to 
%     theta oscillations in human neocortex. Supporting Online Material. 
%     Science. 313(5793):1626-8.
%     www.sciencemag.org/cgi/content/full/313/5793/1626/DC1
%
% Author: Tim Mullen, SCCN/INC/UCSD.

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

EEG = arg_extract(varargin,'EEG',1);
MyComponentNames = [];
if ~isempty(EEG)
    EEG = EEG(1);
    if isfield(EEG,'CAT') && isfield(EEG.CAT,'curComponentNames') && ~isempty(EEG.CAT.curComponentNames)
        MyComponentNames = EEG.CAT.curComponentNames;
    elseif ~isfield(EEG,'CAT')
        MyComponentNames = 1:size(EEG.icaact,1);
        MyComponentNames = strtrim(cellstr(num2str(MyComponentNames'))');
    else
        MyComponentNames = EEG.CAT.curComps;
        MyComponentNames = strtrim(cellstr(num2str(MyComponentNames'))');
    end
    
    clear EEG
end

g = arg_define([0 1], varargin, ...
        arg_norep({'EEG'},mandatory),...
        arg({'ampband','AmplitudePassBand'},[80 150],[],'The [lo hi] pass-band to use for amplitude (Hz)','shape','row'), ...
        arg({'phaseband','PhasePassBand'},[3 7],[],'The [lo hi] pass-band to use for phase (Hz)','shape','row'), ...
        arg({'numsurrogate','NumberOfSurrogateValues'},200,[],'Number of surrogate values for computing statistics') ...
    );


minskip=g.EEG.srate;                   %% time lag must be at least this big
maxskip=g.EEG.CAT.pnts-g.EEG.srate;    %% time lag must be smaller than this 

data = g.EEG.CAT.srcdata;
[nchs npnts ntr] = size(data);

if ntr > 1
    error('SIFT:est_cfc_pacmi:trials_not_supported','Multi-trial (epoched) data not currently supported for PAC');
end

skip=ceil(npnts.*rand(g.numsurrogate*2,1)); 
skip(skip>maxskip)=[]; 
skip(skip<minskip)=[]; 
skip=skip(1:g.numsurrogate,1); 
surrogate_m=zeros(g.numsurrogate,nchs);

% HG analytic amplitude time series
amplitude=abs(hilbert(eegfilt(data,g.EEG.srate,g.ampband(1),g.ampband(2))')); 
% theta analytic phase time series
phase=angle(hilbert(eegfilt(data,g.EEG.srate,g.phaseband(1),g.phaseband(2))'));
z=amplitude.*exp(1i*phase);
% mean of z over time, prenormalized value 
m_raw=mean(z);

% compute surrogate values 
for s=1:g.numsurrogate
    surrogate_amplitude=[amplitude(skip(s):end,:) ; amplitude(1:skip(s)-1,:)]; 
    surrogate_m(s,:)=abs(mean(surrogate_amplitude.*exp(1i*phase))); 
    disp(g.numsurrogate-s)
end

% fit gaussian to surrogate data, uses normfit.m from MATLAB Statistics toolbox 
[surrogate_mean,surrogate_std]=normfit(surrogate_m); 

% normalize length using surrogate data (z-score) 
m_norm_length   = (abs(m_raw)-surrogate_mean)./surrogate_std; 
m_norm_phase    = angle(m_raw);
m_norm          = m_norm_length.*exp(1i*m_norm_phase);


