
function [Conn peaks] = hlp_collapseConn(varargin)
% Apply a set of selection rules to a connectivity matrix
%
% Inputs
%
% Name:                   Information                                                                                     
% ------------------------------------------------------------------------------------------------------------------------
% Conn:                   SIFT Connectivity object. Can also be a 'PConn' (bootstrap) structure.
%                         Input Data Type: struct 
%
% ConnectivityMethods:    Connectivity method names                                                                       
%                         Cell array of connectivity method names.                                                        
%                         Input Data Type: boolean                                                                        
%                                                                                                                         
% DimensionToCollapse:    Dimensions to collapse                                                                          
%                         Input Data Type: string                                                                         
% --------------------                                                                                                    
%                                                                                                                         
%     Frequency:          Collapse across frequency dimension                                                             
%                         Input Data Type: boolean                                                                        
%     ----------                                                                                                          
%                                                                                                                         
%         Range:          Value range [min max]                                                                           
%                         Can also an [N x 2] matrix where each row contains a [min max] range to select.                 
%                         Input Data Type: real number (double)                                                           
%                                                                                                                         
%         CollapseMethod: Collapse method                                                                                 
%                         Possible values: {'sum','net','mean','median','peak','peak2d','getrange','maxmag','max','min'}  
%                         Default value  : 'net'                                                                          
%                         Input Data Type: string                                                                         
%                                                                                                                         
%         Dimension:      Measure dimension                                                                               
%                         This determines the dimension of Conn.(methodname) to collapse. If empty, we will try to        
%                         automatically determine dimension from Conn.dims                                                
%                         Input Range  : [1  Inf]                                                                         
%                         Default value: 3                                                                                
%                         Input Data Type: real number (double)                                                           
%                                                                                                                         
%         Order:          Order to apply transformation                                                                   
%                         Possible values: {'',''}                                                                      
%                         Default value  : ''                                                                            
%                         Input Data Type: real number (double)                                                           
%                                                                                                                         
%     Time:               Collapse across time dimension                                                                  
%                         Input Data Type: boolean                                                                        
%     -----                                                                                                               
%                                                                                                                         
%         Range:          Value range [min max]                                                                           
%                         Can also an [N x 2] matrix where each row contains a [min max] range to select.                 
%                         Input Data Type: real number (double)                                                           
%                                                                                                                         
%         CollapseMethod: Collapse method                                                                                 
%                         Possible values: {'sum','net','mean','median','peak','peak2d','getrange','maxmag','max','min'}  
%                         Default value  : 'net'                                                                          
%                         Input Data Type: string                                                                         
%                                                                                                                         
%         Dimension:      Measure dimension                                                                               
%                         This determines the dimension of Conn.(methodname) to collapse. If empty, we will try to        
%                         automatically determine dimension from Conn.dims                                                
%                         Input Range  : [1  Inf]                                                                         
%                         Default value: 4                                                                                
%                         Input Data Type: real number (double)                                                           
%                                                                                                                         
%         Order:          Order to apply transformation                                                                   
%                         Possible values: {'',''}                                                                      
%                         Default value  : ''                                                                            
%                         Input Data Type: real number (double)                                                           
%                                                                                                                         
% Verbosity:              Verbose output                                                                                  
%                         Input Data Type: boolean 
%
% Outputs
%
% Name:                   Information                                                                                     
% ------------------------------------------------------------------------------------------------------------------------
% Conn:                   Collapsed SIFT Connectivity object.
% 
% peaks:                  Structure containing locations of peaks (if applicable)
%                         peaks.time and peaks.freqs are cell arrays
%                         where the kth element is the locations of peaks
%                         for the kth time or frequency bin
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
% Author: Tim Mullen 2008, 2010, SCCN/INC, UCSD.
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


Conn = arg_extract(varargin,'Conn',[],[]);
cnames = hlp_microcache('hlp_collapseConn',@hlp_getConnMethodNames,Conn);
if ~isstruct(Conn) || isempty(Conn)
    Conn = struct([]);
end

g = arg_define([0 Inf], varargin, ...
    arg_norep({'Conn','Connectivity'},mandatory,[],'Connectivity structure. Can also be a PConn structure.'), ...
    arg({'connmethods','ConnectivityMethods'},cnames,cnames,'Connectivity method names. Cell array of connectivity method names.'), ...
    arg_sub({'coldim','DimensionToCollapse'},{}, ...
    { ...
        arg_subtoggle({'freq','Frequency','Freq'},'off', ...
            { ...
            arg({'range','Range'},[],[],'Value range [min max]. Can also an [N x 2] matrix where each row contains a [min max] range to select.','shape','matrix'), ...
            arg({'method','CollapseMethod'},'net',{'sum','net','mean','median','peak','peak2d','getrange','maxmag','max','min'},'Collapse method'), ...
            arg({'dim','Dimension'},3,[1 Inf],'Measure dimension. This determines the dimension of Conn.(methodname) to collapse. If empty, we will try to automatically determine dimension from Conn.dims') ...
            arg({'order','Order'},1,{1 2},'Order to apply transformation'), ...
            }, 'Collapse across frequency dimension'), ...
         arg_subtoggle({'time','Time'},'off', ...
            { ...
            arg({'range','Range'},[],[],'Value range [min max]. Can also an [N x 2] matrix where each row contains a [min max] range to select.','shape','matrix') ...
            arg({'method','CollapseMethod'},'net',{'sum','net','mean','median','peak','peak2d','getrange','maxmag','max','min'},'Collapse method'), ...
            arg({'dim','Dimension'},4,[1 Inf],'Measure dimension. This determines the dimension of Conn.(methodname) to collapse. If empty, we will try to automatically determine dimension from Conn.dims') ...
            arg({'order','Order'},2,{1 2},'Order to apply transformation'), ...
            }, 'Collapse across time dimension'), ...
        },'Dimensions to collapse.'), ...
    arg({'verb','Verbosity'},false,[],'Verbose output') ...
    );
    
g.peakWidth = 2;  %FIXME: change to option
Conn = g.Conn;

collapsetime = g.coldim.time.arg_selection;
collapsefreq = g.coldim.freq.arg_selection;

if ~collapsetime && ~collapsefreq
    peaks = [];
    return;
end

% Determine freq/time dimensions
if collapsetime && ~isempty(g.coldim.time.dim)
    dim = g.coldim.time.dim;
elseif isfield(Conn,'dims')
    dim=find(ismember_bc(Conn.dims,'time'));
else
    error('hlp_collapseConn could not automatically determine the dimension for ''time''. Perhaps it is the 4th dimension?');
end
if collapsefreq && ~isempty(g.coldim.freq.dim)
    dim = g.coldim.freq.dim;
elseif isfield(Conn,'dims')
    dim=find(ismember_bc(Conn.dims,'freq'));
else
    error('hlp_collapseConn could not automatically determine the dimension for ''freqs''. Perhaps it is the 3rd dimension?');
end
% Determine default freq and time range
if ~collapsefreq || isempty(g.coldim.freq.range)
    g.coldim.freq.range = [Conn.freqs(1) Conn.freqs(end)];
end
if ~collapsetime || isempty(g.coldim.time.range)
    g.coldim.time.range = [Conn.erWinCenterTimes(1) Conn.erWinCenterTimes(end)];
end
if length(Conn.freqs)==1
    Conn.freqs = repmat(Conn.freqs,1,2);
end
if length(Conn.erWinCenterTimes)==1
    Conn.erWinCenterTimes = repmat(Conn.erWinCenterTimes,1,2);
    Conn.winCenterTimes   = repmat(Conn.winCenterTimes,1,2);
end
if ~collapsefreq || any(isnan(g.coldim.freq.range(:)))
    g.coldim.freq.range = [];
end
if ~collapsetime || any(isnan(g.coldim.time.range(:)))
    g.coldim.time.range = [];
end
% get the indices for the time/frequency range
frangeidx = zeros(size(g.coldim.freq.range));
for k=1:size(g.coldim.freq.range,1)
    frangeidx(k,:) = getindex(Conn.freqs,g.coldim.freq.range(k,:));
end
trangeidx = zeros(size(g.coldim.time.range));
for k=1:size(g.coldim.time.range,1)
    trangeidx(k,:) = getindex(Conn.erWinCenterTimes,g.coldim.time.range(k,:));
end

% determine sequence of collapse operations
colseq  = setdiff_bc(fieldnames(g.coldim),'arg_direct');
colseq  = colseq(cellfun(@(fn) g.coldim.(fn).arg_selection,colseq));
% reorder operations by preferred order
seqorder = zeros(1,length(colseq));
for k=1:length(colseq), seqorder(k) = g.coldim.(colseq{k}).order; end
if length(unique_bc(seqorder))<length(seqorder)
    error('SIFT:hlp_collapseConn', ...
          'Please specify a unique value of ''Order'' for each dimension'); 
end
colseq = hlp_vec(colseq(seqorder));

% dimorder = g.coldim.dimorder;
% dimorder = regexprep(dimorder,{'[fF]requency','[Ff]req'},'freq');
% dimorder = regexprep(dimorder,{'[tT]ime'},'time');

if isempty(g.connmethods)
    g.connmethods = hlp_getConnMethodNames(Conn); end
if ~iscell(g.connmethods)
    g.connmethods = {g.connmethods}; end

% Collapse each conn method individually
for m=1:length(g.connmethods)
    connmethod = g.connmethods{m};
   
    C = Conn.(connmethod);
    peaks.time = [];
    peaks.freqs = [];
    
    for dn = 1:length(colseq)
        dim_name = colseq{dn};
        dim_idx = g.coldim.(dim_name).dim;
        col_method = g.coldim.(dim_name).method;
        
        % collapse time dim
        if collapsetime && (isequal(dim_name,'time')) && ~isempty(trangeidx)
            if g.verb
                fprintf('applying method=%s to dim=%s\n',dim_name,dim_idx);
            end
            % iterate over sets of time range indices (bands)
            [Ctmp, ptmp] = deal(cell(1,size(trangeidx,1)));
            for k=1:size(trangeidx,1)
                [Ctmp{k}, peaks.time{k}] = collapseDim(C,'dim',dim_idx,...
                    'range',trangeidx(k,:),'method',col_method,...
                    'dx',Conn.erWinCenterTimes(2)-Conn.erWinCenterTimes(1),'peakWidth',g.peakWidth, ...
                    'minpeakheight',-Inf);
            end
            % replace original data matrix 
            % (FIXME: can optimize with cell2mat(Ctmp) when collapsing over
            % freq)
            C = [];
            for k=1:length(Ctmp)
                C = cat(dim_idx,C,Ctmp{k});
            end
        end

        % collapse freq dim
        if collapsefreq && ~isempty(frangeidx) && (isequal(dim_name,'freq'))
            if g.verb
                fprintf('applying method=%s to dim=%s\n',dim_name,dim_idx);
            end
            [Ctmp, ptmp] = deal(cell(1,size(frangeidx,1)));
            for k=1:size(frangeidx,1)
                [Ctmp{k}, peaks.freqs{k}] = collapseDim(C,'dim',dim_idx,...
                    'range',frangeidx(k,:),'method',col_method,...
                    'dx',Conn.freqs(2)-Conn.freqs(1),'peakWidth',g.peakWidth, ...
                    'minpeakheight',-Inf);
            end
            C = [];
            for k=1:length(Ctmp)
                C = cat(dim_idx,C,Ctmp{k});
            end
        end
    end
            
    % combine results
    
    sz=size(C);
    if all(sz(3:end)==1)
        C = squeeze(C);
    end
    
    Conn.(connmethod) = C;
    
end

if ~isempty(frangeidx) && ~strcmpi(g.coldim.freq.method,'getrange')
    Conn.collapsedFreqs = [];
    if size(frangeidx,1) > 1
        for k=1:size(frangeidx,1)
            Conn.collapsedFreqs{k} = Conn.freqs(frangeidx(k,1):frangeidx(k,end));
        end
        Conn.freqs = cellfun(@(x) median(x),Conn.collapsedFreqs);
    else
        Conn.collapsedFreqs = Conn.freqs(frangeidx(1):frangeidx(end));
        Conn.freqs = median(Conn.collapsedFreqs);
    end
    % frequency dimension has been collapsed
    % replace with median frequency within each collapsed interval
    % (FIXME: THIS ASSUMES EQUAL SPACING BETWEEN FREQS)
    % Conn.freqs = median(Conn.freqs(frangeidx),2)';
end

if ~isempty(trangeidx) && ~strcmpi(g.coldim.time.method,'getrange')
    Conn.collapsedTimes = [];
    Conn.collapsedWinCenterTimes = [];
    if size(trangeidx,1) > 1
        for k=1:size(trangeidx,1)
            Conn.collapsedTimes{k} = Conn.erWinCenterTimes(trangeidx(k,1):trangeidx(k,end));
            Conn.collapsedWinCenterTimes{k} = Conn.winCenterTimes(trangeidx(k,1):trangeidx(k,end));
        end
        Conn.erWinCenterTimes   = cellfun(@(x) median(x),Conn.collapsedTimes);
        Conn.winCenterTimes     = cellfun(@(x) median(x),Conn.collapsedWinCenterTimes);
    else
        Conn.collapsedTimes = Conn.erWinCenterTimes(trangeidx(1):trangeidx(end));
        Conn.erWinCenterTimes   = median(Conn.collapsedTimes);
        Conn.winCenterTimes     = median(Conn.collapsedWinCenterTimes);
    end
    % time dimension has been collapsed
    % replace with median time point within each collapsed interval 
    % (FIXME: THIS ASSUMES EQUAL SPACING BETWEEN TIMES)
    % Conn.erWinCenterTimes   = median(Conn.erWinCenterTimes(trangeidx),2)';
    % Conn.winCenterTimes     = median(Conn.winCenterTimes(trangeidx),2)';
end
