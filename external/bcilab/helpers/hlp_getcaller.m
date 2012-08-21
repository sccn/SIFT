function [name,file] = hlp_getcaller()
% Find the name & file of the calling function.
% [Name,File] = hlp_getcaller()
%
% Out: 
%   Name: MATLAB name of the calling function, if any
%   File: file name of the calling function, if any
%   
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-04-15

try
    throw; %#ok<LTARG> % fast way to get an exception
catch
    e = lasterror; %#ok<LERR>
    if length(e.stack) > 2
        name = e.stack(3).name;
        file = e.stack(3).file;    
    else
        name = [];
        file = [];
    end
end
