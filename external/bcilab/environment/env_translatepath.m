function filename = env_translatepath(filename)
% Translates platform-independent directories into a system-specific directories.
% System-Path = hlp_translatepath(Independent-Path)
%
% BCILAB supports platform-independent paths for all its scripts and IO functions, which allows for script portability from
% system to system. It is especially important when data paths are mounted at different locations, depending on access mode and operating system.
% A side effect of portable path names is that there is one unique expression which computes any given dataset, such as 
% "flt_iir(flt_reref((flt_resample(io_loadset('data:/Projects/Test/test1.vhdr'),[],200))),[],[4 6 25 30])", and this in turn 
% allows to share the same dataset caches (which are indexed by expression) across machines and operating systems, minimizing 
% redundant computations.
%
% In:
%   Independent-Path : platform-independent path; may contain forward and/or backward slashes
%                      (forward slashes generally preferred), and may refer to locations such as
%                      store:/ (the store path) or data:/ (one of the data paths);
%                      can also be a relative path
%
% Out:
%   System-Path : system-specific path (with slashes corrected and locations resolved)
%                 if multiple data paths are present, the one where the minimum number
%                 of directories (and files) would have to be created to write to the 
%                 given file is selected.
%   
% Examples:
%   env_translatepath('data:/projects/test.mat');
%    -> can translate into, e.g. 'C:\Programme\MATLAB\work\mydata\projects\test.mat'
%
% See also:
%   env_startup()
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-06-29

global bcilab;

% turn the path into a system-dependent one
filename = strrep(strrep(filename,'\',filesep),'/',filesep);

% resolve location references
if strmatch('store:',filename)
    filename = [bcilab.paths.store_path filename(1+length('store:'):end)]; 
elseif strmatch('resources:',filename)
    filename = [bcilab.paths.resource_path filename(1+length('resources:'):end)]; 
elseif strmatch('temp:',filename)
    filename = [bcilab.paths.temp_path filename(1+length('temp:'):end)]; 
elseif strmatch('data:',filename)
    rest = filename(1+length('data:'):end);
    bestpath = 1; bestlen = -1;
    if length(bcilab.paths.data_paths) > 1
        fpieces = explode(rest,filesep);
        % find the data path that contains the longest prefix of the filename
        for pidx=1:length(bcilab.paths.data_paths)
            p = bcilab.paths.data_paths{pidx};
            % for each prefix of the filename (starting with the longest one)
            for k=length(fpieces):-1:0
                % check if the data path plus the first k pieces of the filename exists
                if exist([p sprintf([filesep '%s'],fpieces{1:k})],'file')
                    % found a match - check if it is a new length record among all our data paths...
                    if k>bestlen
                        bestlen = k;
                        bestpath = pidx;
                    end
                    break;
                end
            end
        end
    end
    % resolve the reference using that data path which matches most of the filename,
    % where, if multiple data paths are equally well suited, the first one of them is taken
    filename = [bcilab.paths.data_paths{bestpath} rest];
end
