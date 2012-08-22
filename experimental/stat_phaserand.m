function [PConn] = stat_phaserand(ALLEEG, nperms, configs, verb,savefname,savefrq,connmethods)

% NOTE: this function is deprecated

% Check inputs
if nargin < 3
    help stat_phaserand;
    return;
end
if nargin < 4
	verb = 2;
end
if nargin < 5
    savefname = '';
end
if nargin < 6
    savefrq = round(nperms/4);  % number of perms between saves
end
if nargin < 7
    connmethods = {};
end


            
if verb==1
   progress('init',sprintf('Randomizing phases (%d/%d)...',0,nperms)); 
   tic
   t0 = toc;
elseif verb==2
    h=statusbar(sprintf('Randomizing phases (x%d)...',nperms));
end

% Perform phase randomization to generate surrogate null distribution
    
for perm=1:nperms

    EEG = ALLEEG;
    
    % preprocess data
    EEG = pre_prepData('ALLEEG',EEG,configs.prepcfg,'verb',0);

    % Create surrogate data by randomizing phases.
    % Multiply each fourier amplitude by e^{iw)
    % where w is a random phase chosen in [0 2pi]
    % (c.f. Theiler, et al 1997)
    % to generate hermitian phase distributions, I extract
    % the phase of a random matrix. This ensures the random
    % spectrum is conjugate symmetric
    EEG.CAT.srcdata = permute(EEG.CAT.srcdata,[2 1 3]);
    [npnts nchs ntr] = size(EEG.CAT.srcdata);

    % NFFT = 2^nextpow2(npnts); size(data,1);

    for tr=1:ntr
        EEG.CAT.srcdata(:,:,tr) = ...
            ifft(abs(fft(EEG.CAT.srcdata(:,:,tr))) ...
            .* exp(1i*angle(fft(rand(npnts,nchs)))), ...
            'symmetric');
    end

    EEG.CAT.srcdata = ipermute(EEG.CAT.srcdata,[2 1 3]);


    % fit model
    EEG = pop_est_fitMVAR(EEG,0,configs.modfitcfg,'verb',0);

    % calculate connectivity
    EEG = pop_est_mvarConnectivity(EEG,configs.conncfg,'verb',0);

    % append causal estimates to the last dimension of respective arrays in Conn
    if isempty(connmethods)
        connfields = hlp_getConnMethodNames(EEG.CAT.Conn);
    else
        connfields = connmethods;
    end

    if perm==1
        for m=1:length(connfields)
            if ~isempty(EEG.CAT.Conn.(connfields{m}))
                PConn.(connfields{m}) = single([]);
            end
        end
    end

    for m=1:length(connfields)
        if ~isempty(EEG.CAT.Conn.(connfields{m}))
            dim = ndims(EEG.CAT.Conn.(connfields{m}));
            PConn.(connfields{m}) = cat(dim+1,PConn.(connfields{m}),single(EEG.CAT.Conn.(connfields{m})));
        end
    end

    % save checkpoint
    if ~isempty(savefname) && ~mod(perm,savefrq)
        save(sprintf('%s.part',savefname),'PConn','-v7.3');
    end

    if verb==1
        te = toc-t0;
        tt = ceil(te*(nperms-perm)/perm);
        progress(perm/nperms,sprintf('Randomizing phases [%d/%d] (EL: %0.1fm / ETA: %0.1fm)',perm,nperms,ceil(te)/60,tt/60));
    elseif verb==2
        statusbar(perm/nperms,h);
    end

end

PConn.winCenterTimes = EEG.CAT.Conn.winCenterTimes;
PConn.erWinCenterTimes = EEG.CAT.Conn.erWinCenterTimes;
PConn.freqs = EEG.CAT.Conn.freqs;


% clean up
if verb==1
    fprintf('\n')
elseif verb==2
    delete(h);
end
   
if ~isempty(savefname), 
    if exist(sprintf('%s.part',savefname),'file')
            delete(sprintf('%s.part',savefname));
    end
end
