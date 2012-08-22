% epoch2continuous (alpha version)
%
% epoch2continuous.m converts an epoched dataset into a continuous one.
% Segments of data will be concatenated using a 'boundary' event.
%
% USAGE
%
% EEG = epoch2continuous(EEG);
%
% Author: Javier Lopez-Calderon
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% July 7, 2011
%
% Feedback would be appreciated and credited.
% NOTE: No ICA fields have been tested until this version.
%

function EEG = epoch2continuous(EEG)
if nargin<1
      help epoch2continuous
      return
end
if isempty(EEG.epoch)
      error('epoch2continuous() only works for epoched data!')
end
% tic;
[xlat, indx]  = unique([EEG.event.latency], 'first');
neegevent  = length(indx);
% new type
typearray  = {EEG.event(indx).type};
% new duration
if isfield(EEG.event, 'duration')
      durarray = {EEG.event(indx).duration};
else
      durarray = num2cell(ones(1,neegevent));
end
% new urevent
urarray  = {EEG.event(indx).urevent};
nepoch   = EEG.trials;
nepopnts = EEG.pnts;
latsamarray = zeros(1, neegevent);
% new continuous latencies
for i=1:neegevent
      ep  = EEG.event(indx(i)).epoch;
      lat = EEG.epoch(ep).eventlatency{EEG.epoch(ep).event == indx(i)};
      latsam =  getindex(EEG.times, lat);
      latsamarray(i) = latsam + (ep-1)*nepopnts;
end
latsamarray = num2cell(latsamarray); % cell
% new boundaries
latbound    = num2cell(nepopnts:nepopnts:nepopnts*nepoch);
nbound      = length(latbound);
boundarray  = repmat({'boundary'}, 1, nbound);
% concatenates events and boundaries info
typearray   = [typearray boundarray];
latsamarray = [latsamarray latbound];
urarray     = [urarray repmat({0},1,nbound)];
durarray    = [durarray repmat({0},1,nbound)];
neegevent   = length(typearray); % new length
% Builts new EEG
EEG.trials  = 1;
EEG.xmin    = 0;
EEG.data    = reshape(EEG.data , EEG.nbchan,nepoch*nepopnts);
EEG.event   = [];
% Events
[EEG.event(1:neegevent).type    ] = typearray{:};
[EEG.event(1:neegevent).latency ] = latsamarray{:};
[EEG.event(1:neegevent).urevent ] = urarray{:};
[EEG.event(1:neegevent).duration] = durarray{:};
EEG.epoch  = [];
EEG.times  = [];
EEG.epochdescription = {};
EEG.reject = [];
EEG.pnts   = length(EEG.data);
EEG.xmax   = length(EEG.data)/EEG.srate;
% check everything
% tproce = toc;
EEG    = eeg_checkset( EEG, 'eventconsistency' );
% fprintf('Complete in %.1f seconds!\n', tproce);