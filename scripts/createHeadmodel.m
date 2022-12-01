% create head model

% try to use the same matrix for the Leadfield
eeglabp = fileparts(which('eeglab'));

sourceModelPath = fullfile(eeglabp, 'functions', 'supportfiles','head_modelColin27_5003_Standard-10-5-Cap339.mat'));
tmpSource   = load('-mat', sourceModelPath);
a = rand(60, 1000);
EEG = pop_importdata('dataformat','array','nbchan',0,'data','a','srate',100,'pnts',0,'xmin',0);
chanlList =  { 'AF7' 'FP1' 'Afz' 'FP2' 'AF8' 'F7' 'F3' 'F1' 'Fz' 'F2' 'F4' 'F8' 'FT7' 'FC5' 'FC3' 'FC1' 'FCz' 'FC2' 'FC4' 'FC6' 'FT8' 'T7' 'C5' 'C3' 'C1' 'Cz' 'C2' 'C4' 'C6' 'T8' 'TP7' 'CP5' 'CP3' 'CP1' 'CPz' 'CP2' 'CP4' 'CP6' 'TP8' 'P7' 'P5' 'P3' 'P1' 'Pz' 'P2' 'P4' 'P6' 'P8' 'P9' 'PO7' 'PO5' 'PO3' 'Poz' 'PO4' 'PO6' 'PO8' 'P10' 'PO9' 'PO9h' 'O1' 'Oz' 'O2' 'PO10h' 'PO10' };
[chans, ind1, ind2 ] =intersect(tmpSource.labels, chanlList);
EEG.chanlocs = struct('labels', chans);
EEG=pop_chanedit(EEG, 'lookup','/System/Volumes/Data/data/matlab/eeglab/plugins/dipfit/standard_BEM/elec/standard_1005.elc');
EEG = pop_dipfit_settings( EEG, 'hdmfile','/System/Volumes/Data/data/matlab/eeglab/plugins/dipfit/standard_BEM/standard_vol.mat','coordformat','MNI',...
    'mrifile','/System/Volumes/Data/data/matlab/eeglab/plugins/dipfit/standard_BEM/standard_mri.mat',...
    'chanfile','/System/Volumes/Data/data/matlab/eeglab/plugins/dipfit/standard_BEM/elec/standard_1005.elc','coord_transform',[0 0 0 0 0 -1.5708 1 1 1] ,'chansel',[1:64] );
%EEG = pop_dipfit_settings( EEG, 'hdmfile','/System/Volumes/Data/data/matlab/eeglab/plugins/dipfit/standard_BEM/standard_vol.mat','coordformat','MNI','mrifile','/System/Volumes/Data/data/matlab/eeglab/plugins/dipfit/standard_BEM/standard_mri.mat','chanfile','/System/Volumes/Data/data/matlab/eeglab/plugins/dipfit/standard_BEM/elec/standard_1005.elc','coord_transform',[0.90458 -16.6716 2.972 0.093457 0.0027509 -1.5727 0.99632 0.9087 0.97075] ,'chansel',[1:64] );
EEG = pop_leadfield(EEG, 'sourcemodel','/System/Volumes/Data/data/matlab/eeglab/functions/supportfiles/head_modelColin27_5003_Standard-10-5-Cap339.mat','sourcemodel2mni',[0 -24 -45 0 0 -1.5708 1000 1000 1000] ,'downsample',1);

tmpModel    = load('-mat', '/System/Volumes/Data/data/matlab/eeglab/plugins/dipfit/standard_BEM/standard_vol.mat');
surfFile      = '/System/Volumes/Data/data/matlab/eeglab/plugins/dipfit/standard_BEM/standard_vol_surfdata.mat';
hdmFile       = '/System/Volumes/Data/data/matlab/eeglab/plugins/dipfit/standard_BEM/standard_vol_head_model.mat';
leadFieldFile = '/System/Volumes/Data/data/matlab/eeglab/plugins/dipfit/standard_BEM/standard_vol_leadfield.mat';
%coordChans = [[EEG.chanlocs.X]; [EEG.chanlocs.Y]; [EEG.chanlocs.Z]; ones(1,64)];
%coordChans = traditionaldipfit(EEG.dipfit.coord_transform)*coordChans;

surfData(3).vertices = tmpSource.cortex.vertices;
surfData(3).faces    = tmpSource.cortex.faces;
save('-mat', surfFile, 'surfData');
sourceRois = {'temporalpole L' 'temporalpole R' }; %'frontalpole R'  'frontalpole L' }; %,'superiorparietal R' 'Cingulum_Mid_L' 'Parietal_Sup_L','Frontal_Sup_Medial_R','Precentral_R'};

leadField.K = tmpSource.K(ind1,:); %[ EEG.dipfit.sourcemodel.leadfield{:} ];
leadField.L = tmpSource.L; %
save('-mat', leadFieldFile, '-struct', 'leadField')
clear metadata;
metadata.surfaces         = surfData;
metadata.surfacesFilename = surfFile;
metadata.fiducials    = [];
metadata.atlas.label  = tmpSource.atlas.label;
metadata.atlas.color  = tmpSource.atlas.colorTable;
metadata.channelSpace = tmpSource.channelSpace(ind1,:)*1000; %  coordChans(1:3,:)';
metadata.label        = { EEG.chanlocs.labels }; % tmpSource.labels(ind1)
metadata.leadField    = tmpSource.K(ind1,:); % not necessary
metadata.L            = tmpSource.L; % not necessary
metadata.leadFieldFile = leadFieldFile; % not necessary
save('-mat', hdmFile, 'metadata');

hlp_microcache('clear');
hmObj = hlp_validateHeadModelObject(hdmFile);
