  
function [C peaks] = hlp_filterConns(EEG, Conn, varargin)
%
% Apply a set of selection rules to a connectivity matrix
%
% Inputs:
%
%   EEG         EEG data structure
%   ````1`
%   Conn        [M x M x num_freqs x num_times] connectivity matrix
%
% Options:
%
%       'connThresh'    -     [real]. Connectivity threshold 
%                             (only keep C<connThresh)
%       'prcThresh'     -     [real] The upper percent [1-100] of filtered 
%                             connections to keep {def: 100}
%       'method'              cell array of dimensionality reduction methods to apply in a specified order. e.g.,
%                             {'freq','net','time','peak'} will first integrate over freq, then find peak over time. 
%                             If method is a string, then this is taken to be the compression global dimension reduction 
%                             method applied to all dims in the order {'time','freq',...}.
%                             allowed methods:  'mean', 'net', 'peak', 'shrinkonly'  
%        ...
%
% Outputs:
%
%   C       filtered connectivity matrix
%   peaks   frequency peaks (if method = 'peak')
%
% See Also: est_mvarConnectivity()
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


% ----------------------------------------------------------
% TODO:  
% need to provide ability to specify order in which dim-compression
% methods are applied. e.g., {3,'net',4,'peak'} would first 
% integrate over the third dimension (frequency) and then find the peaks
% over the 4th dim (time)
% ONE OPTION:  allow 'method' = a cell array as above
% ----------------------------------------------------------
% Key point, need to have EEG.CAT.dims = {'chans','chans','freq','times'}
% or some such...
%


if nargin<1
    help hlp_filterConns;
    return;
end


if length(size(Conn))>4
    error('Connectivity Matrix cannot be greater than 4-D!');
end

% parse inputs
g = finputcheck(varargin,...
    {'connThresh'   ''          []          0; ...  % absolute thresholding to apply after significance thresholding. can be a scalar or a matrix of same size as EEG.CAT.C. Can be logical C(C~=thresh) = 0 or real-valued C(C<thresh)=0
     'prcThresh'    'real'      [eps 100]   100;...
     'frange'       'real'      []          [];...  % same units as EEG.CAT.freqs
     'trange'       'real'      []          [];...  % same units as EEG.times
     'dtrange'      'real'      []          [];...
     'sigThresh'    ''          []          [];...  % can be a scalar or a matrix of same size as EEG.CAT.C. Can be logical C(C~=thresh) = 0 or real-valued C(C<thresh)=0
     'badchans'     'integer'   []          [];...
     'method'       ''          ''          'mean';...        % cell array of dimensionality reduction methods to apply in a specified order. e.g.,
                                                              % {'freq','net','time','peak'} will first integrate over freq, then find peak over time. 
                                                              % If method is a string, then this is taken to be the compression global dimension reduction 
                                                              % method applied to all dims in the order {'time','freq',...}
     'peakWidth'    'integer'   []          2;...
     'chanlocs'     ''          []          [];...
     'distRad'      'real'      []          [];...
     'strictRad'    'boolean'   []          1;...  
     'metric',      'string'    []          'manhattan';...
     'diags',       'string'    {'off','on'} 'off'; ...
     'freqdim',     'real'      []          []; ...
     'timedim',     'real'      []          []; ...
     },'hlp_filterConns','error','quiet');

if ischar(g)
    error(g);
end

nmethods = 3;

if ~isempty(g.timedim)
    timedim = g.timedim;
else
    timedim=find(ismember(EEG.CAT.dims,'time'));
end

if ~isempty(g.freqdim)
   freqdim = g.freqdim; 
else
    freqdim=find(ismember(EEG.CAT.dims,'freq'));
end

if isempty(freqdim)
    disp('warning: no ''freq'' entry found in EEG.CAT.dims -- assuming dim 3');
    freqdim = 3;
end
if isempty(timedim)
    disp('warning: no ''time'' entry found in EEG.CAT.dims -- assuming dim 4');
    timedim = 4;
end

if ischar(g.method), 
    ONEMETHOD = true;
    g.method = {'time',g.method,'freq',g.method,'delay',g.method};
else
    ONEMETHOD = false;
end

% TODO: implement dtrange
if ~isempty(g.dtrange)
    warning('time delays not implemented yet!');
    g.dtrange = [];
end

% if ~isfield(EEG.CAT,'sigdist') || isempty(EEG.CAT.sigdist)
%     g.sigThresh = [];
% end

if ~isfield(EEG.CAT,'freqs') || isempty(EEG.CAT.freqs)
    g.frange = [];  % no freq. dim
elseif isempty(g.frange)
    g.frange = [EEG.CAT.freqs(1) EEG.CAT.freqs(end)]; 
end

if ~isfield(EEG.CAT,'times') || isempty(EEG.CAT.times)
    g.trange = [];  % no time dim
elseif isempty(g.trange)
        g.trange = [EEG.CAT.times(1) EEG.CAT.times(end)];
end

C = EEG.CAT.C;  % make copy of connectivity matrix


% apply significance thresholding
if ~isempty(g.sigThresh)
    if islogical(g.sigThresh) && isequal(size(C),size(g.sigThresh))
        C(~g.sigThresh)=0;
    else
        C(C<g.sigThresh) = 0;
    end
%     warning('significance thresholding not implemented yet!');
end

if ~strcmpi(g.method{2},'thresh')
    if ~isempty(g.trange) && ~isempty(g.frange) && ONEMETHOD && isequal(g.method{2},'peak2d')
        % apply 2D peak identification
        fprintf('applying method=%s to dim=%s and %s\n',g.method{2},g.method{1},g.method{3});
        C = collapseDim(C,'dim',4,'range',getindex(EEG.CAT.times,g.trange), ...
            'method','shrinkonly');     % extract time range 
        C = collapseDim(C,'dim',3,'range',getindex(EEG.CAT.freqs,g.frange), ...
            'method','shrinkonly');     % extract freq range
        [C peaks] = collapseDim(C,'dim',4,'range',[], 'method','peak2d', ...
            'peakWidth',g.peakWidth,'minpeakheight',g.connThresh); % find peaks
    else

        for mm=1:2:length(g.method)

            % collapse time dim
            if ~isempty(g.trange) && (isequal(g.method{mm},'time'))
                fprintf('applying method=%s to dim=%s\n',g.method{mm+1},g.method{mm});
                [C peaks.times] = collapseDim(C,'dim',timedim,...
                    'range',getindex(EEG.CAT.times,g.trange),'method',g.method{mm+1},...
                    'dx',EEG.CAT.times(2)-EEG.CAT.times(1),'peakWidth',g.peakWidth, ...
                    'minpeakheight',g.connThresh);
            end

            % collapse freq dim
            if ~isempty(g.frange) && (isequal(g.method{mm},'freq'))
                fprintf('applying method=%s to dim=%s\n',g.method{mm+1},g.method{mm});
                [C peaks.freqs] = collapseDim(C,'dim',freqdim,...
                    'range',getindex(EEG.CAT.freqs,g.frange),'method',g.method{mm+1},...
                    'dx',EEG.CAT.freqs(2)-EEG.CAT.freqs(1),'peakWidth',g.peakWidth, ...
                    'minpeakheight',g.connThresh);
            end

            % filter delays
            if ~isempty(g.dtrange) && (isequal(g.methods{mm},'delay')) 
                % TODO: implement time delay thresholding
            end

        end

    end
end

sz=size(C);
if all(sz(3:end)==1)
    C = squeeze(C);
end


% ensure C is now 2-D and square
if length(size(C))>2 || size(C,1)~=size(C,2)
    warning('Connectivity matrix was not successfully collapsed to 2-D square mat')
else

    % remove diagonals
    if strcmpi(g.diags,'off')
        C=C.*~eye(size(C));
    end
    
    % apply connectivity threshold
    if g.connThresh
        %         g.connThresh(isspace(g.connThresh))=[]; % rem whitespace
        %         if strcmp(g.connThresh(end),'%')        % rem trailing symbol
        %             g.connThresh(end) = [];
        %         end
        %         g.connThresh = str2num(g.connThresh);
        
        C(C<g.connThresh) = 0;
        
    end
    
    % apply distance radius
    if g.distRad
        [C] = disthresh(C, 'R',g.distRad,'chanlocs',g.chanlocs,...
                                         'strictRad',g.strictRad,...
                                         'metric',g.metric);
    end


    % apply percentile threshold
    if ~isempty(g.prcThresh) && g.prcThresh<100
            % get specified percentile
            prc = prctile(C(:),100-g.prcThresh);
            C(C<prc)=0;
    end
    
end



 
         