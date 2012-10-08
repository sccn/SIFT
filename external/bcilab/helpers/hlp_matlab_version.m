function v = hlp_matlab_version()
% Get the MATLAB version in a numeric format that can be compared with <, >, etc.
persistent vers;
if isempty(vers)
    vers = explode(version,'.'); vers = str2num(vers{1})*100 + str2num(vers{2}); end
v = vers;