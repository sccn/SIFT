function Conn = hlp_sift_emptyconn(varargin)
% create a new (empty) Conn datastructure with default fields initialized
% varargin is an optional set of <name,value> pairs indicating fields and 
% associated values to store into the set

if nargin > 1 && mod(length(varargin),2)
    error('SIFT:hlp_sift_emptyconn', ...
        ['Additional arguments must be a list of <name,value> pairs. ' ...
         'At least one name is missing its corresponding value.']);
end

Conn.winCenterTimes     = [];
Conn.erWinCenterTimes   = [];
Conn.freqs              = [];
Conn.dims               = {'var_to','var_from','freq','time'};

% append additional arguments
if nargin > 1
    Conn = hlp_varargin2struct(varargin,Conn);    
end