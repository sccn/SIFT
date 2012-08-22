
cfg.timeRange = [-1 -0.5];
vis_causalBrainMovie3D(EEG,EEG.CAT.Conn,');
figh = gcf;

cfg.timeRange = [1 1];
cfg.BMopts.figurehandle = figh;
for i=1:length(EEG.CAT.Conn.erWinCenterTimes)
    Conntmp = struct('dDTF',EEG.CAT.Conn.dDTF(:,:,:,i), ...
              'winCenterTimes',1, ...
              'erWinCenterTimes',1, ...
              'freqs', EEG.CAT.Conn.freqs);
    vis_causalBrainMovie3D(EEG,Conntmp,cfg);
end



