function v = hlp_gcc_version()
% Get the GCC version in a numeric format that can be compared with <, >, etc.
[a,b]=system(['gcc --version']); v= str2num(strrep(b(11:15),'.','')); %#ok<ASGLU>