function hash = hlp_cryptohash(data,fromfile)
% Compute an MD5 hash of a file, string or generic data structure.
% Hash = hlp_cryptohash(Data,FromFile)
%
% In:
%   Data : data to be hashed; can be a filename, a string, or any other MATLAB data type.
%
%   FromFile : if true, data is interpreted as a file name (default: false)
%
% Out:
%   Hash : MD5 hash (decimal, lowercase) of the Data
%
%						Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                       2010-10-10


if exist('fromfile','var') && fromfile
    % take data as a file name
    if ~ischar(data)
        error('To represent a file name, Data should be a string.'); end
    f = fopen(data,'r');
    try
        data = fread(f,Inf);
        fclose(f);
    catch
        fclose(f);
        rethrow(lasterror); %#ok<LERR>
    end
else
    % take data literally
    if ~ischar(data)
        data = hlp_serialize(data); end
end

% use Java to hash the data (idea from Michael Kleder)
hasher = java.security.MessageDigest.getInstance('MD5');
hasher.update(uint8(data));
hash = dec2hex(typecast(hasher.digest,'uint8'),2)';
hash = lower(hash(:)');