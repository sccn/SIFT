function [newMin newMax] = hlp_scaleLimits(varargin)
% this function implements an exponential window moving average to
% adapt a set of [min max] limits based on incoming data

arg_define([0 4],varargin, ...
    arg_norep({'values'},mandatory,[],'data values'), ...
    arg_norep({'lastMin','LastMin'},mandatory,[],'last mininum value'), ...
    arg_norep({'lastMax','LastMax'},mandatory,[],'last maximum value'), ...
    arg_norep({'numberOfRunsSoFar'},mandatory,[],'number of previous calls to this function'), ...
    arg({'adaptationHL'},10,[0 Inf],'half-life of exponential win moving average. In frames'), ...
    arg({'updateInterval'},1,[],'num frames between updates'), ...
    arg({'bufferTime'},10,[],'num frames to wait before adapting') ...
    );

% memory factor for exponential window moving average
MEMFACTOR = 2/((adaptationHL * 2.8854)+1);

newMin = lastMin;
newMax = lastMax;

if (numberOfRunsSoFar < bufferTime || mod(numberOfRunsSoFar,updateInterval) == 0) && numberOfRunsSoFar < Inf
    
    if numberOfRunsSoFar < bufferTime
        % min/max are based on absolute min/max for this window
        newMin = min(values(:));
        newMax = max(values(:));
    else
        % apply exponential-window moving average to calculate new min/max
        newMin  = MEMFACTOR * min(values(:)) + (1-MEMFACTOR) * lastMin;
        newMax  = MEMFACTOR * max(values(:)) + (1-MEMFACTOR) * lastMax;
    end
    
end