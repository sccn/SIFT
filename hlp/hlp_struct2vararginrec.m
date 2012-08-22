function var = hlp_struct2vararginrec(cfg)
% recursively replace all structs with varargin
% see also: EEGLAB:vararg2str()
% Author: Tim Mullen, 2012, SCCN/INC/UCSD

var = hlp_struct2varargin(cfg);

for i=1:length(var),
    if isstruct(var{i}),
        var{i} = hlp_struct2vararginrec(var{i});
    end
end
