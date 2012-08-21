


function Conn = est_tfpac(EEG,varargin)
% Calculate the event-related phase-amplitude coupling for a collection of 
% component activations.


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


% 
% g = finputcheck(varargin, ...
%                 { 'alpha'         'real'     [0 0.2]                  [];
%                   'baseboot'      'float'    []                       0;
%                   'boottype'      'string'   {'times','trials','timestrials'}  'timestrials';
%                   'detrend'       'string'   {'on','off'}              'off';
%                   'freqs'         'real'     [0 Inf]                  [0 srate/2];
%                   'freqs2'        'real'     [0 Inf]                  [];
%                   'freqscale'     'string'   { 'linear','log' }       'linear';
%                   'itctype'       'string'   {'phasecoher','phasecoher2','coher'}  'phasecoher';
%                   'nfreqs'        'integer'  [0 Inf]                  [];
%                   'lowmem'        'string'   {'on','off'}              'off';
%                   'method'        'string'   { 'mod','corrsin','corrcos','latphase' }         'mod';
%                   'naccu'         'integer'  [1 Inf]                   250;
%                   'newfig'        'string'   {'on','off'}              'on';
%                   'padratio'      'integer'  [1 Inf]                   2;
%                   'rmerp'         'string'   {'on','off'}              'off';
%                   'rboot'         'real'     []                        [];
%                   'subitc'        'string'   {'on','off'}              'off';
%                   'subwin'        'real'     []                        []; ...
%                   'gammapowerlim' 'real'     []                        []; ...
%                   'powerlim'      'real'     []                        []; ...
%                   'powerlat'      'real'     []                        []; ...
%                   'gammabase'     'real'     []                        []; ...
%                   'timesout'      'real'     []                        []; ...
%                   'ntimesout'     'integer'  []                        200; ...
%                   'tlimits'       'real'     []                        [0 frame/srate];
%                   'title'         'string'   []                        '';
%                   'vert'          { 'real','cell' }  []                [];
%                   'wavelet'       'real'     [0 Inf]                   0;
%                   'wavelet2'      'real'     [0 Inf]                   [];
%                   'winsize'       'integer'  [0 Inf]                   max(pow2(nextpow2(frame)-3),4) }, 'pac');
%               
%               
%               

% EEG.pnts=size(EEG.data,2);     %% number of sample points in raw signal
minskip=g.EEG.srate;             %% time lag must be at least this big
maxskip=g.EEG.pnts-g.EEG.srate;    %% time lag must be smaller than this 

g.EEG.icaact = g.EEG.icaact(g.components,:);
[M N] = size(g.EEG.icaact);

skip=ceil(g.EEG.pnts.*rand(g.numsurrogate*2,1)); 
skip(skip>maxskip)=[]; 
skip(skip<minskip)=[]; 
skip=skip(1:g.numsurrogate,1); 
surrogate_m=zeros(g.numsurrogate,M);

% HG analytic amplitude time series, uses eegfilt.m from EEGLAB toolbox 
% (http://www.sccn.ucsd.edu/eeglab/) 
amplitude=abs(hilbert(eegfilt(g.EEG.icaact,g.EEG.srate,g.ampband(1),g.ampband(2))')); 
% theta analytic phase time series, uses EEGLAB toolbox 
phase=angle(hilbert(eegfilt(g.EEG.icaact,g.EEG.srate,g.phaseband(1),g.phaseband(2))')); %% complex-valued composite signal 
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


