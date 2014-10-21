function hDipPlot = test_plotGroupDipoles(HeadModelFileName,GroupEEG,ALLEEG,ClustOrder,Coloring,ClustCentrTraj,ClustCentrTrue,CovScaleFact)
    % CentrTraj are the estimated cluster centroids (x,y,z) over MCMC
    % iterations, stored in GroupEEG.CAT.MCMC_LastState.(methodname).S_BAR
    
    if nargin<3 || isempty(ALLEEG)
        PlotAllDipoles = false;
    else
        PlotAllDipoles = true;
    end
    if nargin<4
        ClustOrder = [];
    end
    if nargin<5
        Coloring = 'clusters'; % 'subjects'
    end
    if nargin<6
        ClustCentrTraj = [];
    end
    if nargin<7
        ClustCentrTrue = [];
    end
    if nargin<8
        CovScaleFact = [];
    end
    
    % get cluster centroids as [M x 3]
    S_BAR   = reshape([GroupEEG.dipfit.model.posxyz],3,[])';
    % get cluster covariance matrices as [3 x 3 x M]
    SIGMA_S = reshape([GroupEEG.dipfit.model.covmat],3,3,[]);
    % reorder sources according to specified ordering
    if ~isempty(ClustOrder)
        S_BAR      = S_BAR(ClustOrder,:);
        SIGMA_S    = SIGMA_S(:,:,ClustOrder);
    end
    
    if PlotAllDipoles
        % get all subject-level dipoles and store in [N x 1] cell array
        SubjDipoles  = arrayfun(@(EEG_i) reshape([EEG_i.dipfit.model.posxyz],3,[])', ...
                                 ALLEEG,'UniformOutput',false);
        if strcmpi(Coloring,'clusters')
            ClusterIDs  = arrayfun(@(EEG_i) [EEG_i.dipfit.model.clustid], ...
                                 ALLEEG,'UniformOutput',false);
        end
        N = length(SubjDipoles);  % number of subjects
    end
    
    M = size(S_BAR,1);   % number of clusters
    
    % centroid dipole size
    centrDipSize = 100;
    % subject dipole sizes
    subjDipSize  = centrDipSize/2;
    % marker sizes for (I),(G),(E)
    MarkerSize = 15;
    % bgcolors = distinguishable_colors(length(hmObj.atlas.label),[1 1 1]);
    bgcolors = [1 1 1; 0 0 0; 0.5 0.5 0.8];
        
    if ~isempty(HeadModelFileName) && ischar(HeadModelFileName)
        % load the head model
        hmObj = headModel.loadFromFile(HeadModelFileName);

        % plot head model
        hDipPlot = hmObj.plotHeadModel;
        hold on;

        % sensors and labels off
        set(hDipPlot.hSensors,'Visible','off');
        set(hDipPlot.hLabels,'Visible','off');
        % scalp off
        set(hDipPlot.hScalp,'FaceAlpha',0);
        % atlas off
        set(hDipPlot.hSkull,'FaceAlpha',0);
        
        % colormap(hDipPlot.hAxes,bgcolors)
        % make cortex semi-transparent
        set(hDipPlot.hCortex,'FaceAlpha',0.2);
        set(hDipPlot.hCortex,'FaceColor',[0.5 0.5 0.8]);
    elseif ishandle(HeadModelFileName)
        figure(HeadModelFileName);
        hDipPlot.hAxes = gca;
        hold on;
%         set(hDipPlot.hAxes,'nextplot','add');
    else
        figure('color','w');
        hDipPlot.hAxes = axes;
%         axis(hDipPlot.hAxes,'off');
        grid(hDipPlot.hAxes,'on');
        axis(hDipPlot.hAxes,'vis3d');
        hold on;
%         set(hDipPlot.hAxes,'nextplot','add');
    end
    
    % determine cluster colors
    centrDipColors = distinguishable_colors(M,bgcolors);
    
    for k=1:M
        % plot cluster centroids
        scatter3(hDipPlot.hAxes,S_BAR(k,1),S_BAR(k,2),S_BAR(k,3),centrDipSize,centrDipColors(k,:),'filled','sizedata',72);

        % plot cluster covariance
        mu    = squish(S_BAR(k,:));
        sigma = squish(SIGMA_S(:,:,k));
        if ~isempty(CovScaleFact)
            sigma = sigma/CovScaleFact(k);
        end
        hcov = plot_gaussian_ellipsoid(mu,sigma,1,100,hDipPlot.hAxes);
        set(hcov,'facealpha',0.3,'edgecolor','none','facecolor',centrDipColors(k,:));
    end
    
    
    if PlotAllDipoles
        switch lower(Coloring)
            case 'subjects'
                subjDipColors = distinguishable_colors(N,[bgcolors; centrDipColors]);
                subjDipColors = mat2cell(subjDipColors,ones(1,N),3);
            case 'clusters'
                for k=1:N
                    subjDipColors{k} = centrDipColors(ClusterIDs{k},:);
                end
        end
        for k=1:N
            % plot subject dipoles
            scatter3(hDipPlot.hAxes,SubjDipoles{k}(:,1),SubjDipoles{k}(:,2),SubjDipoles{k}(:,3),subjDipSize,subjDipColors{k},'filled','sizedata',72/2);
        end
    end

    if ~isempty(ClustCentrTraj)
        for k=1:M
            % plot cluster centroids for each MCMC iteration (trajectories)...
            plot3(hDipPlot.hAxes,squish(ClustCentrTraj(k,1,:)),squish(ClustCentrTraj(k,2,:)),squish(ClustCentrTraj(k,3,:)),'linewidth',0.1,'linestyle','none','marker','.','markersize',8,'color',centrDipColors(k,:));
            % plot (I) marker for initial conditions
            ctrs = squish(ClustCentrTraj(k,:,1));
            tt=text(ctrs(1),ctrs(2),ctrs(3),'I','FontSize',MarkerSize,'Color','k','FontName','Times','HorizontalAlignment','center','VerticalAlignment','middle');
            set(tt,'Parent',hDipPlot.hAxes);
            plot3(hDipPlot.hAxes,ctrs(1),ctrs(2),ctrs(3),'o','MarkerSize',MarkerSize,'LineWidth',2,'color','k');
            % plot (E) marker for mean centroid location
            ctrs = S_BAR(k,:);
            plot3(hDipPlot.hAxes,ctrs(1),ctrs(2),ctrs(3),'o','MarkerSize',MarkerSize,'LineWidth',2,'color','k');
            plot3(hDipPlot.hAxes,ctrs(1),ctrs(2),ctrs(3),'x','MarkerSize',MarkerSize,'LineWidth',2,'color','k');
        end
    end

    if ~isempty(ClustCentrTrue)
        for k=1:M
            % plot (G) marker for ground truth
            ctrs = squish(ClustCentrTrue(k,:));
            tt=text(ctrs(1),ctrs(2),ctrs(3),'G','FontSize',MarkerSize,'Color','k','FontName','Times','HorizontalAlignment','center','VerticalAlignment','middle');
            set(tt,'Parent',hDipPlot.hAxes);
            plot3(hDipPlot.hAxes,ctrs(1),ctrs(2),ctrs(3),'o','MarkerSize',MarkerSize,'LineWidth',2,'color','k');
        end
    end
    
    if ~isempty(HeadModelFileName) && ischar(HeadModelFileName)
        % cleanup
        hmObj.saveToFile(HeadModelFileName);
        delete(hmObj.surfaces);
        delete(hmObj.leadFieldFile);
    end