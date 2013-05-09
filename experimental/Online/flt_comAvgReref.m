function EEG = flt_comAvgReref(EEG,goodchannels)

    % do common average rereferencing (excluding bad channels)
    EEG.data = bsxfun(@minus,EEG.data,mean(EEG.data(goodchannels,:),1));