
% retrieve specific study dataset(s) from ALLEEG
MYEEG1 = ALLEEG(STUDY.design(STUDY.currentdesign).cell(1).dataset);  % cond 1
MYEEG2 = ALLEEG(STUDY.design(STUDY.currentdesign).cell(2).dataset);  % cond 2

% retreive only trials matching study design
ALLEEG(STUDY.design(STUDY.currentdesign).cell(2).dataset).icaact(:,:, STUDY.design(STUDY.currentdesign).cell(2).trials{1})

% extract data with selected trials from a study
eeg_getdatact(ALLEEG(STUDY.design(STUDY.currentdesign).cell(2).dataset), 'component', 3, 'trialindices', STUDY.design(STUDY.currentdesign).cell(2).trials{1});
eeg_getdatact(ALLEEG(STUDY.design(STUDY.currentdesign).cell(2).dataset), 'component', [ 3 5 6], 'trialindices', STUDY.design(STUDY.currentdesign).cell(2).trials{1});
dat = eeg_getdatact(ALLEEG(STUDY.design(STUDY.currentdesign).cell(2).dataset), 'component', [ 3 5 6], 'trialindices', STUDY.design(STUDY.currentdesign).cell(2).trials{1});
