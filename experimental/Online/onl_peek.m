function buf = onl_peek(streamname,len)
% Peek into an online stream (generates an EEG-set like view into it)
%
% In:
%   Streamname : Name of the online stream in the workspace
%
%   ViewLength : length of the view that should be generated; should not be longer than the
%                buffer capacity, in seconds (default: 10)
%
% Out:
%   Signal : An EEGLAB data set that represents an view into the stream's most recent 
%            ViewLength seconds
%
% Example:
%   % get the last 5 seconds of the stream
%   EEG = onl_peek('mystream',5)
%
% See also:
%   onl_newpredictor, onl_newstream, onl_append
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-03

if ~exist('streamname','var')
    error('You need to pass at least the name of a previously created stream.'); end
if ~exist('len','var')
    len = 10; end

try
    % get the stream from the base workspace
    stream = evalin('base',streamname);
catch e
    if strcmp(e.identifier,'MATLAB:badsubscript')
        error('BCILAB:onl_predict:improper_resolve','The raw data required by the predictor does not list the name of the needed source stream; this is likely a problem in onl_newpredictor');
    else
        error('BCILAB:onl_predict:stream_not_found',['The stream named ' streamname ' was not found in the base workspace.']);
    end
end


% by extracting only the new samples, rawdata behaves like a stateful pipeline stage
samples_to_get = min(stream.buffer_len, round(stream.srate*len));
% copy buffer and read out .data field from circular .buffer field
buf = stream;
buf.data = stream.buffer(:, 1+mod(stream.smax-samples_to_get:stream.smax-1,stream.buffer_len));
[buf.nbchan,buf.pnts,buf.trials] = size(buf.data);
buf.xmax = buf.smax/buf.srate;
buf.xmin = buf.xmax - (buf.pnts-1)/buf.srate;
