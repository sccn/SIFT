%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SCRIPTING EXAMPLE FOR THE SOURCE INFORMATION FLOW TOOLBOX (SIFT)    %%%
%%% SIFT Version: 0.5-alpha                                             %%%
%%%                                                                     %%%
%%% This example demonstrates how to use SIFT from the command-line or  %%%
%%% in a script. This example applies to SIFT 0.5-alpha.                %%%
%%% For additional information on the below steps, please consult the   %%%
%%% SIFT manual located at http://sccn.ucsd.edu/wiki/SIFT               %%%
%%% Author: Tim Mullen (C) 2011, SCCN, INC, UCSD                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% STEP 1: Load Data

% We will begin by loading up the 'RespWrong.set' dataset located in the 
% /Data/ folder within the Sample Data package
% (you can download this package from the SIFT website or at 
% ftp://sccn.ucsd.edu/pub/tim/SIFT/SIFT_SampleData.zip)

EEG = pop_loadset;

%% OPTIONAL STEP: Analyzing Channels

% If you want to analyze channel data rather than components you can use
% this "hack". Namely, we will replace the ICA soution with a "fake" 
% solution that always copies the channel data into the ICA component
% (icaact) field. This ensures that SIFT will always work on the channel
% data rather than component data. NOTE: this functionality will be
% replaced by a more flexible interface in the Beta version.

h=msgbox('Next the EEGLAB options GUI will pop up. Uncheck the option labeled "If set, scale ICA components to RMS in microvolt" and click OK');
waitfor(h);
pop_editoptions;


% Now we will replace the ICA soution with a "fake" solution that always
% copies the channel data into the ICA component (icaact) field
% WARNING, THIS WILL *DELETE* YOUR ICA SOLUTION, uncomment the following
% lines to create a backup of the current ICA solution

% EEG.etc.icaweights_beforeIdentity = EEG.icaweights;
% EEG.etc.icasphere_beforeIdentity = EEG.icasphere;
% EEG.etc.icawinv_beforeIdentity = EEG.icawinv;

[EEG.icaweights EEG.icasphere EEG.icawinv] = deal(eye(EEG.nbchan));
EEG.icaact = [];

% copy the data into ICA 'activations' and verify dataset
EEG = eeg_checkset(EEG,'ica');


%% STEP 2: Define key Processing Parameters

Components = [8 11]; %13 19 20 23 38 39];   % these are the components/channels to which we'll fit our multivariate model
WindowLengthSec = 0.5;                   % sliding window length in seconds
WindowStepSizeSec = 0.25;                % sliding window step size in seconds
NewSamplingRate = [];                    % new sampling rate (if downsampling)


%% STEP 3: Pre-process the data

ComponentsToKeep = strtrim(cellstr(num2str(Components')));  % convert list of components to cell array of strings

[EEG prepcfg] = pre_prepData('ALLEEG',EEG,'VerbosityLevel',2,'NewSamplingRate',NewSamplingRate,'NormalizeData',{'Method',{'ensemble'}},'SelectComponents',{'ComponentsToKeep',ComponentsToKeep});

%% STEP 4: Identify the optimal model order

% Here we compute various model order selection criteria for varying model
% orders (e.g. 1 to 30) and visualize the results

% compute model order selection criteria...
IC = pop_est_selModelOrder(EEG,0,'icselector',{'aic','sbc','hq','ris'},'algorithm','vieira-morf','morder',[1 30],'winlen',WindowLengthSec,'winstep',WindowStepSizeSec,'prctWinToSample',100,'verb',1,'plot',0);

% ... and plot the results
handles = vis_plotOrderCriteria(IC);

% If you want to save this figure you can uncomment the following lines:
%
% for i=1:length(handles)
%     saveas(handles(i),sprintf('orderResults%d.fig',i));
% end
% close(handles);

% Finally, we can automatically select the model order which minimizes one
% of the criteria (or you can set this manually based on above figure)
ModelOrder = ceil(mean(IC{1}.hq.popt));

% As an alternative to using the minimum of the selection criteria over 
% model order, you can find the "elbow" in the plot of model order versus
% selection criterion value. This is useful in cases where the selection
% criterion does not have a clear minimum. For example, the lines below
% plot and select the elbow location (averaged across windows) for the AIC 
% criterion
%
% vis_plotOrderCriteria(IC,{},{},'elbow');
% ModelOrder = ceil(mean(IC{1}.aic.pelbow));



%% STEP 5: Fit the VAR model

% Once we have identified our optimal model order, we can fit our VAR model.

% Fit a model using a sliding-window approach with the vieira-morf
% lattice filter algorithm
[EEG modfitcfg] = pop_est_fitMVAR(EEG,0,'algorithm','vieira-morf','morder',ModelOrder,'winlen',WindowLengthSec,'winstep',WindowStepSizeSec,'verb',1);


% Note that EEG.CAT.MODEL now contains the model structure with
% coefficients (in MODEL.AR), prediction errors (MODEL.PE) and other
% self-evident information

% Alternately, we can fit the VAR parameters using a Kalman filter (see
% doc est_fitMVARKalman for more info on arguments)
%
% EEG.CAT.MODEL = est_fitMVARKalman(EEG,0,'updatecoeff',0.0005,'updatemode',2,'morder',ModelOrder,'verb',2,'downsampleFactor',50);


%% STEP 6: Validate the fitted model

% Here we assess the quality of the fit of our model w.r.t. the data. This
% step can be slow.

% We can obtain statistics for residual whiteness, percent consistency, and
% model stability ...
[whitestats PC stability] = pop_est_validateMVAR(EEG,0,'checkWhiteness',true,'whitenessCriteria',{'Ljung-Box','ACF','Box-Pierce','Li-McLeod'},...
                                                'checkConsistency',true,'checkStability',true,'alpha',0.05,'prctWinToSample',20,'verb',2,'plot',0);

% ... and then plot the results
handles = vis_plotModelValidation(whitestats,PC,stability);

% If you want to save this figure you can uncomment the following lines:
%
% for i=1:length(handles)
%     saveas(handles(i),sprintf('validationResults%d.fig',i));
% end
% close(handles);


% To automatically determine whether our model accurately fits the data you
% can write a few lines as follows (replace 'acf' with desired statistic):
%
% if ~all(whitestats{1}.acf.w)
%     msgbox('Residuals are not completely white!');
% end


%% STEP 7: Compute Connectivity

% Next we will compute various dynamical quantities, including connectivity,
% from the fitted VAR model. We can compute these for a range of
% frequencies (here 3-45 Hz). See 'doc est_mvarConnectivity' for a complete
% list of available connectivity and spectral estimators.

[EEG conncfg] = pop_est_mvarConnectivity(EEG,'verb',true,'freqs',(3 : 45),'connmethods',{'dDTF08','nDTF','Coh','S','pCoh','nPDC'},'absvalsq',true,'spectraldecibels',true);


%% STEP 8: Compute Statistics

% reload the datasets
EEGfresh = pop_loadset;

% first we obtain the bootstrap distributions for each condition
[PConn(1)] = stat_bootstrap(EEGfresh(1), 100, struct('conncfg',conncfg,'modfitcfg',modfitcfg,'prepcfg',prepcfg(1)));
[PConn(2)] = stat_bootstrap(EEGfresh(2), 100, struct('conncfg',conncfg,'modfitcfg',modfitcfg,'prepcfg',prepcfg(2)));

% next we compute the between-condition pvalues
Stats = stat_bootSigTest(PConn,'fdr');


%% STEP 9: Visualize the Connectivity estimates in a Time-Frequency Grid

[figureHandles tfgridcfg] = vis_TimeFreqGrid('ALLEEG',EEG,'Conn',EEG.CAT.Conn,'Stats',Stats,'MatrixLayout',{'Partial','UpperTriangle', 'dDTF08', 'LowerTriangle','dDTF08','Diagonal','S'},'ColorLimits',99.9,'Baseline',[-1.75 -0.5],'Smooth2D',true);

% You can also partially populate the GUI via a call to the pop_ function:
%
%[figureHandles tfgridcfg] = pop_vis_TimeFreqGrid(EEG,'MatrixLayout',{'Partial','UpperTriangle', 'dDTF08', 'LowerTriangle','dDTF08','Diagonal','S'},'ColorLimits',99.9,'Baseline',[-1.75 -0.5],'Smooth2D',true);


%% STEP 10: Visualize the Connectivity estimates in a 3D Brain-Movie

cfg=pop_vis_causalBrainMovie3D(EEG,'ConnectivityMethod','dDTF08','FreqCollapseMethod','integrate','FrequenciesToCollapse',(3:8),'NodeColorMapping','CausalFlow','FooterPanelDisplaySpec',{'ICA_ERPenvelope',{'icaenvelopevars', '8'},{'backprojectedchans' 'B2'}},'BrainMovieOptions',{'RenderCorticalSurface',{'VolumeMeshFile' 'standard_BEM_vol.mat', 'Transparency' 0.7}});    %


