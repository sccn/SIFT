% ---------------
% Plotting curves
% ---------------
% R is the freq x 1 matrix for chA --> chB causality
% Rstat is the freq x [low high] bootstrap significance threshold mat
% times is vector of times
% freqs is vector of freqs
% mbase is the vector (1 x freqs) of the mean of the baseline
% g is struct containing parameters and options
%
% TODO :   re-write this function to plot time curves

function plotallcurves(R, Rboot, times, freqs, mbase, g)

	% compute angles
	% --------------
	Rangle = angle(R);
	pos = get(gca,'position'); % plot relative to current axes
	q     = [pos(1) pos(2) 0 0];
	s = [pos(3) pos(4) pos(3) pos(4)];
	if ~isreal(R)
        R = abs(R);
        Rraw =R; % raw coherence values
        if ~isnan(g.baseline)
         	R = R - repmat(mbase',[1 g.timesout]); % remove baseline mean
        end;
	else
        Rraw = R;
        setylim = 0;
	end;

    % time unit
    % ---------
    if times(end) > 10000
        times = times/1000;
        timeunit = 's';
    else
        timeunit = 'ms';
    end;

    % legend
    % ------
    alllegend = {};
    if strcmpi(g.plotmean, 'on') && freqs(1) ~= freqs(end)
      alllegend = { [ num2str(freqs(1)) '-' num2str(freqs(end)) 'Hz' ] };
    else
      for index = 1:length(freqs)
         alllegend{index} = [ num2str(freqs(index)) 'Hz' ];
      end;
    end;
    
	fprintf('\nNow plotting...\n');
	if strcmpi(g.plotamp, 'on')
      %
      % Plot coherence amplitude in top panel
      %
      if strcmpi(g.plotphase, 'on'), subplot(2,1,1); end; 
      if isempty(g.maxamp), g.maxamp = 0; end;
      plotcurve(times, R, 'maskarray', Rboot, 'title', 'Coherence amplitude', ...
                'xlabel', [ 'Time (' timeunit ')' ], 'ylabel', '0-1', 'ylim', g.maxamp, ...
                'vert', g.vert, 'marktimes', g.marktimes, 'legend', alllegend, ...
                'linewidth', g.linewidth, 'highlightmode', g.highlightmode, 'plotmean', g.plotmean);
	end;
	
	if strcmpi(g.plotphase, 'on')
      %
      % Plot coherence phase lags in bottom panel
      %
      if strcmpi(g.plotamp, 'on'), subplot(2,1,2); end; 
      plotcurve(times, Rangle/pi*180, 'maskarray', Rboot, 'val2mask', R, 'title', 'Coherence phase', ...
                'xlabel', [ 'Time (' timeunit ')' ], 'ylabel', 'Angle (deg.)', 'ylim', [-180 180], ...
                'vert', g.vert, 'marktimes', g.marktimes, 'legend', alllegend, ...
                'linewidth', g.linewidth, 'highlightmode', g.highlightmode, 'plotmean', g.plotmean);
	end
	
	if strcmpi(g.plotamp, 'on') | strcmpi(g.plotphase, 'on')
       try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;
       if (~isempty(g.title)) % plot title
          h(13) = textsc(g.title, 'title');
       end
	
       %
       %%%%%%%%%%%%%%% plot topoplot() %%%%%%%%%%%%%%%%%%%%%%%
       %
       if (~isempty(g.topovec))
          h(15) = subplot('Position',[-.1 .43 .2 .14].*s+q);
          if size(g.topovec,2) <= 2
             topoplot(g.topovec(1),g.elocs,'electrodes','off', ...
                'style', 'blank', 'emarkersize1chan', 10);
          else
             topoplot(g.topovec(1,:),g.elocs,'electrodes','off');
          end;
          axis('square')
          
          h(16) = subplot('Position',[.9 .43 .2 .14].*s+q);
          if size(g.topovec,2) <= 2
             topoplot(g.topovec(2),g.elocs,'electrodes','off', ...
                'style', 'blank', 'emarkersize1chan', 10);
          else
             topoplot(g.topovec(2,:),g.elocs,'electrodes','off');
          end;
          axis('square')
       end
       
       try, axcopy(gcf); catch, end;
	end;