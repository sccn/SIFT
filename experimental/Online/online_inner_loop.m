% grab the last 500 ms of data from the stream
    EEG = onl_peek('mystream',WindowLengthSec);
    
    
    % insert missing fields (replace with initialization in run_readbiosemi)
    f = fieldnames(s);
    for i = 1:length(f)
        EEG.(f{i}) = s.(f{i});
    end
    
    % pre-process in SIFT
    %     EEG = pre_prepData('ALLEEG',EEG,'VerbosityLevel',0,'NewSamplingRate',NewSamplingRate,'NormalizeData',{'Method',{'time'}},'SelectComponents',{'ComponentsToKeep',ComponentsToKeep});
    %     EEG.times = linspace(EEG.xmin*1000,EEG.xmax*1000,EEG.pnts);
    
    EEG = onl_init_SIFT(EEG,Components,ComponentsToKeep);
    
    tic
    
    % fit MVAR model
    EEG = pop_est_fitMVAR(EEG,0,'normalize',{'time'},'algorithm','gladmm','morder',ModelOrder,'winlen',WindowLengthSec,'winstep',WindowStepSizeSec,'epochTimeLims',[EEG.xmin EEG.xmax],'verb',0,'gladmm',struct('lambda',0.0005,'verb',0));
    
    % get connectivity
    EEG.CAT.Conn = est_mvarConnectivity(EEG,EEG.CAT.MODEL,'connmethods',connmethods,'freqs',1:40,'verb',0);
   
    tfit = toc;
    
    
    if CHECK_WHITENESS
        [whitestats PC stability] = pop_est_validateMVAR(EEG,0,'checkWhiteness',true,'whitenessCriteria',{'Ljung-Box','ACF'},...
                                                'checkConsistency',false,'checkStability',false,'alpha',0.05,'prctWinToSample',100,'verb',0,'plot',0);
        fprintf('\n\n------------------------------\n');
        fprintf('whiteness test:\n');
        fprintf('ACF: p=%0.10g (%s)\n',whitestats{1}.acf.pval,fastif(whitestats{1}.acf.w,'white','not white'));
        fprintf('LBP: p=%0.10g (%s)\n',whitestats{1}.ljungbox.pval,fastif(whitestats{1}.ljungbox.w,'white','not white'));
        fprintf('------------------------------\n\n');
    end
    
    if init
        % initialize brainmovie
        hbmcfg = gui_causalBrainMovie3D_online(EEG,EEG.CAT.Conn,struct('arg_direct',0));
        waitfor(hbmcfg,'UserData','init');
        %         BMCFG.BMopts.speedy = true;
        bmoptschanged = true;
        init = false;
        figh = [];
        bmvars = [];
        
        % cross-fade feature
%         fadetimer = timer;


        
%         eeglines=plot(dataax,EEG.times,squeeze(EEG.CAT.srcdata(1,:,:)));
%         
%         eegplot(EEG.CAT.srcdata,'srate',EEG.srate);
%         eeglines = get(findobj(gcf,'tag','eegaxis'),'children');
        
%         OldConn = EEG.CAT.Conn;
    else
%         figh = BMCFG.BMopts.figurehandle;
        
%         norm(OldConn.dDTF08(:)-EEG.CAT.Conn.dDTF08(:));
%         OldConn = EEG.CAT.Conn;
    end
    
    if bmoptschanged
        mode = 'init_and_render';
        bmoptschanged = false;
        bmvars = [];
%         BMCFG
        
    end
    
    %     fprintf('data norm: %0.5g\n', norm(EEG.data(:)));
    
    %     norm(EEG.CAT.Conn.dDTF08(:))
    
%     plot(ax,EEG.CAT.Conn.freqs, squeeze(EEG.CAT.Conn.dDTF08(1,3,:,1)));
    
    
    % render brainmovie
     
    tic
    
    opts = hlp_mergeVarargin(BMCFG.BMopts,'figurehandle',figh,'mode',mode,'speedy',true,'vars',bmvars);
    [tmp handles] = vis_causalBrainMovie3D(EEG,EEG.CAT.Conn,BMCFG,'timeRange',[],'BMopts',opts);
    
    tdraw = toc;
    
    % testing cross-fade feature
    
%     if ~isempty(figh)
% %         surfaces = findobj(get(findobj(figh,'tag','brain1'),'children'),'type','surface'); % all surfaces in brain
% 
%         surfaces = findobj(figh,'tag','tmpmov','type','surface'); % all edges/nodes
%         
%         stop(fadetimer);
%         set(fadetimer,'Period',0.01, ...  % round(1000*((tdraw+tfit)/10))/1000
%             'TimerFcn',{@(obj,evnt,handles,factor) hlp_scaleAlpha(handles,factor),surfaces,0.8});
%         start(fadetimer);
%         
% %         hlp_scaleAlpha(surfaces,0.5);
%     end
        

%     set(eeglines,'ydata',EEG.CAT.srcdata(1,:,:));
    
%     EEG.data = EEG.CAT.srcdata;
%     vis_hist(EEG);

%     eegplot(EEG.CAT.srcdata,'srate',EEG.srate);
    
%     for i=1:length(eeglines)
%         set(eeglines(i),'ydata',EEG.CAT.srcdata(i,:,:));
%     end
    
    figh = handles.figurehandle;
        
    bmvars = tmp.BMopts.vars;
    
    mode = 'render';
    
    
    if BENCHMARK
        fprintf('model fitting: %0.5g\n',tfit);
        fprintf('rendering    : %0.5g\n',tdraw);
        fprintf('total        : %0.5g\n',tfit+tdraw);
    end