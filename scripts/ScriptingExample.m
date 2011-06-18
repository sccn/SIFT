

%% If using channels:

pop_editoptions;
% Uncheck the option labeled "If set, scale ICA components to RMS in microvolt"
% Click "OK"

% Now we will replace the ICA soution with a "fake" solution that always
% copies the channel data into the ICA component (icaact) field
% WARNING, THIS WILL *DELETE* YOUR ICA SOLUTION
EEG.etc.icaweights_beforeIdentity = EEG.icaweights;
EEG.etc.icasphere_beforeIdentity = EEG.icasphere;
EEG.etc.icawinv_beforeIdentity = EEG.icawinv;

[EEG.icaweights EEG.icasphere EEG.icawinv] = deal(eye(EEG.nbchan));

EEG.icaact = [];

% recompute the ICA activations
EEG = eeg_checkset(EEG,'ica');
ALLEEG(CURRENTSET) = EEG;



%% SETUP PARAMETERS

Components = [3 5 6 7 14 18 19 28]; %[6 7 10 11 12 14 17 22 25];
WindowLengthSec = 0.3;
WindowStepSizeSec = 0.02;



%% Preprocess the data
ComponentsToKeep = strtrim(cellstr(num2str(Components')));

[EEG config] = pre_prepData('ALLEEG',ALLEEG(1),'VerbosityLevel',2,'NewSamplingRate',256,'NormalizeData',{'Method',{'time','ensemble'}},'SelectComponents',{'ComponentsToKeep',ComponentsToKeep});

% gui version
% EEG2 = pop_pre_prepData(ALLEEG(1),0,'VerbosityLevel',2,'NormalizeData',{'Method',{'time','ensemble'}},'SelectComponents',{'ComponentsToKeep',{'10','13'}});

%% identify optimal model orders
IC = pop_est_selModelOrder(EEG,0,'icselector',{'aic','sbc','hq'},'algorithm','vieira-morf','morder',[1 30],'winlen',WindowLengthSec,'winstep',WindowStepSizeSec,'prctWinToSample',10,'verb',1,'plot',0);

% % plot the results
handles = vis_plotOrderCriteria(IC);
% for i=1:length(handles)
%     saveas(handles(i),sprintf('orderResults%d.fig',i));
% end
% close(handles);


% pick an optimal model order
ModelOrder = ceil(mean(IC{1}.sbc.popt));

%% fit the MVAR model
[EEG cfg] = pop_est_fitMVAR(EEG,0,'algorithm','vieira-morf','morder',ModelOrder,'winlen',WindowLengthSec,'winstep',WindowStepSizeSec,'verb',1);

% ALLEEG.CAT.MODEL now contains the model parameters

%% (Optional) Validate the fitted model
[whitestats PC stability] = pop_est_validateMVAR(EEG,0,'checkWhiteness',true,'whitenessCriteria',{'Ljung-Box','ACF','Box-Pierce','Li-McLeod'},...
                                                'checkConsistency',true,'checkStability',true,'alpha',0.05,'prctWinToSample',20,'verb',2,'plot',0);

% % plot the results
handles = vis_plotModelValidation(whitestats,PC,stability);
% for i=1:length(handles)
%     saveas(handles(i),sprintf('validationResults%d.fig',i));
% end
% close(handles);

%% Compute connectivity

% use the dDTF08 (has good freq resolution and produces sparse connectivity
% graph and normalizes across all time,freq,and pairs to allow better
% comparison between causal values at different time points, freqs, pairs)
EEG = pop_est_mvarConnectivity(EEG,'verb',true,'freqs',[3 : 45],'connmethods',{'DTF','dDTF','dDTF08','ffDTF','nDTF','GGC', 'iCoh','Coh','S','pCoh','mCoh','GPDC','nPDC','PDCF','PDC'},'absvalsq',true,'spectraldecibels',true);

%% Now we can visualize whatever we want
% for cond=1:length(ALLEEG)
%     pop_vis_TimeFreqGrid(ALLEEG(cond));
% end
cfg = pop_vis_TimeFreqGrid(EEG,cfg);





