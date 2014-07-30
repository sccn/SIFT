classdef headModel < handle & matlab.mixin.Copyable
    properties(GetAccess=public, SetAccess=public,SetObservable)
        channelSpace = [];
        fiducials = [];
        surfaces = []; % 1) scalp, 2) skull, 3) brain (gray matter or average between gray and white matter)
        atlas
        leadFieldFile = [];
        surfNormal = [];
    end
    properties(GetAccess = public, SetAccess = public)
        label;
    end
    methods
        function obj = headModel(varargin)
            if length(varargin)==1, varargin = varargin{1};end
            
            if ~iscell(varargin) && ischar(varargin) && exist(varargin,'file')
                [obj.channelSpace,obj.label,obj.fiducials] = readMontage(varargin);
                return
            end
            
            ind = find(ismember_bc(varargin(1:2:length(varargin)-1),'channelSpace'));
            if ~isempty(ind), obj.channelSpace = varargin{ind*2};end
            
            ind = find(ismember_bc(varargin(1:2:length(varargin)-1),'surfaces'));
            if ~isempty(ind), obj.surfaces = varargin{ind*2};end
            
            ind = find(ismember_bc(varargin(1:2:length(varargin)-1),'atlas'));
            if ~isempty(ind), obj.atlas = varargin{ind*2};end
            
            ind = find(ismember_bc(varargin(1:2:length(varargin)-1),'fiducials'));
            if ~isempty(ind), obj.fiducials = varargin{ind*2};end
            
            ind = find(ismember_bc(varargin(1:2:length(varargin)-1),'leadFieldFile'));
            if ~isempty(ind), obj.leadFieldFile = varargin{ind*2};end
            
            ind = find(ismember_bc(varargin(1:2:length(varargin)-1),'label'));
            if ~isempty(ind),
                obj.label = varargin{ind*2};
            else
                N = size(obj.channelSpace,1);
                labels = num2str((1:N)');
                obj.label =num2cell(labels',[97,1])';
                for it=1:N, obj.label{it} = deblank(obj.label{it});end
            end
            
            ind = find(ismember_bc(varargin(1:2:length(varargin)-1),'surfNormal'));
            if ~isempty(ind), obj.surfNormal = varargin{ind*2};end
        end
        %%
        function h = plotHeadModel(obj,~) % do not remove the circumflex, I'm passging a second arguments when this method is called from MoBILAB's gui
            if isempty(obj.channelSpace) || isempty(obj.label) || isempty(obj.surfaces);
                error('Head model is incomplete or missing.');
            end
            h = headModelViewerHandle(obj,obj.label);
        end
        function labels = getChannelLabels(obj)
            labels = obj.label;
        end
        %%
        function h = plotMontage(obj,showNewfig)
            if isempty(obj.channelSpace) || isempty(obj.label);error('MoBILAB:noChannelSpace','Channel space is empty.');end
            if nargin < 2, showNewfig = true;end
            
            if isa(obj,'eeg')
                color = [0.93 0.96 1];
            else
                color = [0.76 0.77 1];
            end
            if showNewfig, figure('Color',color);end
            h = scatter3(obj.channelSpace(:,1),obj.channelSpace(:,2),obj.channelSpace(:,3),'filled',...
                'MarkerEdgeColor','k','MarkerFaceColor','y','parent',gca);
            
            hold on;
            N = length(obj.label);
            k = 1.1;
            for it=1:N, text('Position',k*obj.channelSpace(it,:),'String',obj.label{it});end
            mx = max(obj.channelSpace);
            k = 1.2;
            line([0 k*mx(1)],[0 0],[0 0],'LineStyle','-.','Color','b','LineWidth',2)
            line([0 0],[0 k*mx(2)],[0 0],'LineStyle','-.','Color','g','LineWidth',2)
            line([0 0],[0 0],[0 k*mx(3)],'LineStyle','-.','Color','r','LineWidth',2)
            text('Position',[k*mx(1) 0 0],'String','X','FontSize',12,'FontWeight','bold','Color','b')
            text('Position',[0 k*mx(2) 0],'String','Y','FontSize',12,'FontWeight','bold','Color','g')
            text('Position',[0 0 k*mx(3)],'String','Z','FontSize',12,'FontWeight','bold','Color','r')
            
            try %#ok
                scatter3(obj.fiducials.nasion(1),obj.fiducials.nasion(2),obj.fiducials.nasion(3),'filled','MarkerEdgeColor','k','MarkerFaceColor','K');
                text('Position',1.1*obj.fiducials.nasion,'String','Nas','FontSize',12,'FontWeight','bold','Color','k');
                
                scatter3(obj.fiducials.lpa(1),obj.fiducials.lpa(2),obj.fiducials.lpa(3),'filled','MarkerEdgeColor','k','MarkerFaceColor','K');
                text('Position',1.1*obj.fiducials.lpa,'String','LPA','FontSize',12,'FontWeight','bold','Color','k');
                
                scatter3(obj.fiducials.rpa(1),obj.fiducials.rpa(2),obj.fiducials.rpa(3),'filled','MarkerEdgeColor','k','MarkerFaceColor','K');
                text('Position',1.1*obj.fiducials.rpa,'String','RPA','FontSize',12,'FontWeight','bold','Color','k');
                
                scatter3(obj.fiducials.vertex(1),obj.fiducials.vertex(2),obj.fiducials.vertex(3),'filled','MarkerEdgeColor','k','MarkerFaceColor','K');
                text('Position',1.1*obj.fiducials.vertex,'String','Ver','FontSize',12,'FontWeight','bold','Color','k');
                
                scatter3(obj.fiducials.inion(1),obj.fiducials.inion(2),obj.fiducials.inion(3),'filled','MarkerEdgeColor','k','MarkerFaceColor','K');
                text('Position',1.1*obj.fiducials.inion,'String','Ini','FontSize',12,'FontWeight','bold','Color','k');
            end
            
            % box on;
            hold off;
            axis equal
            grid on;
        end
        %%
        function individualHeadModelFile = warpTemplate2channelSpace(obj,headModelFile,individualHeadModelFile)
            
            if nargin < 2, error('Reference head model is missing.');end
            if isempty(obj.channelSpace) || isempty(obj.label), error('Channel space or labels are missing.');end
            if ~isempty(obj.fiducials)
                individualHeadModelFile = warpTemplate2channelSpaceUsingFiducials(obj,headModelFile,individualHeadModelFile);
                return;
            end
            if ~exist(headModelFile,'file'), error('The file you''ve entered does not exist.');end
            if nargin < 3, individualHeadModelFile = tempname;end
            
            template = load(headModelFile);
            [~,loc1,loc2] = intersect_bc(obj.label,template.label);
            if isempty(loc1), error('MoBILAB:noLandmarks','No common landmarks found between the template and the individual channels.');end
            gTools = geometricTools;
            th = norminv(0.90);
            % mapping source to target spaces: S->T
            % target space: individual geometry
            T = obj.channelSpace(loc1,:);
            
            % source space: template
            S = template.channelSpace(loc2,:);;
            
            if isa(obj,'eeg'), obj.container.container.initStatusbar(1,8,'Co-registering...');end
            
            % affine co-registration
            Aff = gTools.affineMapping(S,T);
            if isa(obj,'eeg'), obj.container.container.statusbar(1);end
            
            % b-spline co-registration (only fiducial landmarks)
            options.Verbose = true;
            options.MaxRef = 2;
            surfData = template.surfData;
            Ns = length(surfData);
            for it=1:Ns
                surfData(it).vertices = gTools.applyAffineMapping(template.surfData(it).vertices,Aff);
            end
            Saff = gTools.applyAffineMapping(S,Aff);
            [Def,spacing,offset] = gTools.bSplineMapping(Saff,T,surfData(1).vertices,options);
            if isa(obj,'eeg'), obj.container.container.statusbar(2);end
            
            % b-spline co-registration (second pass)
            for it=1:Ns
                surfData(it).vertices = gTools.applyBSplineMapping(Def,spacing,offset,surfData(it).vertices);
            end
            T = obj.channelSpace;
            T(T(:,3) <= min(surfData(1).vertices(:,3)),:) = [];
            [S,d] = gTools.nearestNeighbor(surfData(1).vertices,T);
            z = zscore(d);
            S(abs(z)>th,:) = [];
            T(abs(z)>th,:) = [];
            [Def,spacing,offset] = gTools.bSplineMapping(S,T,surfData(1).vertices,options);
            if isa(obj,'eeg'), obj.container.container.statusbar(3);end
            
            % b-spline co-registration (third pass)
            for it=1:Ns
                surfData(it).vertices = gTools.applyBSplineMapping(Def,spacing,offset,surfData(it).vertices);
            end
            T = obj.channelSpace;
            T(T(:,3) <= min(surfData(1).vertices(:,3)),:) = [];
            [S,d] = gTools.nearestNeighbor(surfData(1).vertices,T);
            z = zscore(d);
            S(abs(z)>th,:) = [];
            T(abs(z)>th,:) = [];
            Tm = 0.5*(T+S);
            [Def,spacing,offset] = gTools.bSplineMapping(S,Tm,surfData(1).vertices,options);
            if isa(obj,'eeg'), obj.container.container.statusbar(4);end
            
            % apply the final transformation
            for it=1:Ns
                surfData(it).vertices = gTools.applyBSplineMapping(Def,spacing,offset,surfData(it).vertices);
            end
            
            % fixing topological defects
            if isa(obj,'eeg'), obj.container.container.statusBar.setText('Fixing topological defects...');end
            dmax = ones(Ns-1,1)*5;
            dmax(1) = 8;
            ind = fliplr(1:Ns);
            for it=1:Ns-1
                surfData(ind(it+1)).vertices = gTools.repareIntersectedSurface(surfData(ind(it)),surfData(ind(it+1)),dmax(it));
                if isa(obj,'eeg'), obj.container.container.statusbar(it+5);end
            end
            
            ind =  obj.channelSpace(:,3) > min(surfData(1).vertices(:,3));
            T = gTools.nearestNeighbor(surfData(1).vertices,obj.channelSpace);
            channelSpace = obj.channelSpace; %#ok
            channelSpace(ind,:) = T(ind,:);  %#ok
            [~,loc] = unique_bc(channelSpace,'rows');%#ok
            indInterp = setdiff_bc(1:size(obj.channelSpace,1),loc);
            if ~isempty(indInterp)
                x = setdiff_bc(channelSpace,channelSpace(indInterp,:),'rows');%#ok
                xi = gTools.nearestNeighbor(x,channelSpace(indInterp,:));%#ok
                channelSpace(indInterp,:) = 0.5*(xi + channelSpace(indInterp,:));%#ok
            end
            obj.channelSpace = channelSpace; %#ok
            
            if isfield(template,'atlas'), obj.atlas = template.atlas;end
            if exist(obj.surfaces,'file'), delete(obj.surfaces);end
            obj.surfaces = individualHeadModelFile;
            save(obj.surfaces,'surfData');
            if isa(obj,'eeg'), obj.container.container.statusbar(8);end
        end
        %%
        function individualHeadModelFile = warpTemplate2channelSpaceUsingFiducials(obj,headModelFile,individualHeadModelFile)
            %if nargin < 2, method = 'tps';end
            if nargin < 2, error('Reference head model is missing.');end
            if isempty(obj.channelSpace) || isempty(obj.label) || isempty(obj.fiducials)
                error('Channel space or fiducials are missing.');
            end
            if ~exist(headModelFile,'file'), error('The file you''ve entered does not exist.');end
            if nargin < 3, individualHeadModelFile = tempname;end
            
            template = load(headModelFile);
            gTools = geometricTools;
            th = norminv(0.90);
            % mapping source to target spaces: S->T
            % target space: individual geometry
            T = [obj.fiducials.nasion;...
                obj.fiducials.lpa;...
                obj.fiducials.rpa];
            
            % source space: template
            S = [template.fiducials.nasion;...
                template.fiducials.lpa;...
                template.fiducials.rpa;...
                template.fiducials.vertex];
            
            % estimates vertex if is missing
            if isfield(obj.fiducials,'vertex')
                if numel(obj.fiducials.vertex) == 3
                    T = [T;obj.fiducials.vertex];
                else
                    point = 0.5*(obj.fiducials.lpa + obj.fiducials.rpa);
                    point = ones(50,1)*point;
                    point(:,3) = linspace(point(3),1.5*max(obj.channelSpace(:,3)),50)';
                    [~,d] = gTools.nearestNeighbor(obj.channelSpace,point);
                    [~,loc] = min(d);
                    point = point(loc,:);
                    T = [T;point];
                end
            else
                point = 0.5*(obj.fiducials.lpa + obj.fiducials.rpa);
                point = ones(50,1)*point;
                point(:,3) = linspace(point(3),1.5*max(obj.channelSpace(:,3)),50)';
                [~,d] = gTools.nearestNeighbor(obj.channelSpace,point);
                [~,loc] = min(d);
                point = point(loc,:);
                T = [T;point];
            end
            
            if isfield(obj.fiducials,'inion')
                if numel(obj.fiducials.vertex) == 3
                    T = [T;obj.fiducials.inion];
                    S = [S;template.fiducials.inion];
                end
            end
            %             else
            %             % estimates inion if is missing
            %                 point = 0.5*(obj.fiducials.lpa + obj.fiducials.rpa);
            %                 b = regress([obj.fiducials.nasion(2); point(2)],[obj.fiducials.nasion(1) 1; point(1) 1]);
            %                 point = ones(50,1)*point;
            %                 if obj.fiducials.nasion(1)>point(1), s = -1;else s = 1;end
            %                 point(:,1) = s*linspace(point(1,1),1000,50)';
            %                 point(:,2) = b(1)*point(:,1) + b(2);
            %                 [~,d] = gTools.nearestNeighbor(obj.channelSpace,point);
            %                 [~,loc] = min(d);
            %                 point = point(loc,:);
            %                 T = [T;point];
            %             end
            if isa(obj,'eeg'), obj.container.container.initStatusbar(1,8,'Co-registering...');end
            
            % affine co-registration
            Aff = gTools.affineMapping(S,T);
            if isa(obj,'eeg'), obj.container.container.statusbar(1);end
            
            % b-spline co-registration (only fiducial landmarks)
            options.Verbose = true;
            options.MaxRef = 2;
            surfData = template.surfData;
            Ns = length(surfData);
            for it=1:Ns
                surfData(it).vertices = gTools.applyAffineMapping(template.surfData(it).vertices,Aff);
            end
            Saff = gTools.applyAffineMapping(S,Aff);
            [Def,spacing,offset] = gTools.bSplineMapping(Saff,T,surfData(1).vertices,options);
            if isa(obj,'eeg'), obj.container.container.statusbar(2);end
            
            % b-spline co-registration (second pass)
            for it=1:Ns
                surfData(it).vertices = gTools.applyBSplineMapping(Def,spacing,offset,surfData(it).vertices);
            end
            T = obj.channelSpace;
            T(T(:,3) <= min(surfData(1).vertices(:,3)),:) = [];
            [S,d] = gTools.nearestNeighbor(surfData(1).vertices,T);
            z = zscore(d);
            S(abs(z)>th,:) = [];
            T(abs(z)>th,:) = [];
            [Def,spacing,offset] = gTools.bSplineMapping(S,T,surfData(1).vertices,options);
            if isa(obj,'eeg'), obj.container.container.statusbar(3);end
            
            % b-spline co-registration (third pass)
            for it=1:Ns
                surfData(it).vertices = gTools.applyBSplineMapping(Def,spacing,offset,surfData(it).vertices);
            end
            T = obj.channelSpace;
            T(T(:,3) <= min(surfData(1).vertices(:,3)),:) = [];
            [S,d] = gTools.nearestNeighbor(surfData(1).vertices,T);
            z = zscore(d);
            S(abs(z)>th,:) = [];
            T(abs(z)>th,:) = [];
            Tm = 0.5*(T+S);
            [Def,spacing,offset] = gTools.bSplineMapping(S,Tm,surfData(1).vertices,options);
            if isa(obj,'eeg'), obj.container.container.statusbar(4);end
            
            % apply the final transformation
            for it=1:Ns
                surfData(it).vertices = gTools.applyBSplineMapping(Def,spacing,offset,surfData(it).vertices);
            end
            
            % fixing topological defects
            if isa(obj,'eeg'), obj.container.container.statusBar.setText('Fixing topological defects...');end
            dmax = ones(Ns-1,1)*5;
            dmax(1) = 8;
            ind = fliplr(1:Ns);
            for it=1:Ns-1
                surfData(ind(it+1)).vertices = gTools.repareIntersectedSurface(surfData(ind(it)),surfData(ind(it+1)),dmax(it));
                if isa(obj,'eeg'), obj.container.container.statusbar(it+5);end
            end
            
            ind =  obj.channelSpace(:,3) > min(surfData(1).vertices(:,3));
            T = gTools.nearestNeighbor(surfData(1).vertices,obj.channelSpace);
            channelSpace = obj.channelSpace; %#ok
            channelSpace(ind,:) = T(ind,:);  %#ok
            [~,loc] = unique_bc(channelSpace,'rows');%#ok
            indInterp = setdiff_bc(1:size(obj.channelSpace,1),loc);
            if ~isempty(indInterp)
                x = setdiff_bc(channelSpace,channelSpace(indInterp,:),'rows');%#ok
                xi = gTools.nearestNeighbor(x,channelSpace(indInterp,:));%#ok
                channelSpace(indInterp,:) = 0.5*(xi + channelSpace(indInterp,:));%#ok
            end
            obj.channelSpace = channelSpace; %#ok
            
            if isfield(template,'atlas'), obj.atlas = template.atlas;end
            if exist(obj.surfaces,'file'), delete(obj.surfaces);end
            obj.surfaces = individualHeadModelFile;
            save(obj.surfaces,'surfData');
            if isa(obj,'eeg'), obj.container.container.statusbar(8);end
        end
        %%
        function obj = plotOnModel(obj,J,V,figureTitle,varargin)
            if nargin < 3, error('Not enough input arguments');end
            if nargin < 3, figureTitle = '';end
            channelLabels = obj.label;
            obj = currentSourceViewer(obj,J,V,figureTitle,channelLabels,varargin{:});
        end
        
        %%
        function Aff = warpChannelSpace2Template(obj,headModelFile,individualHeadModelFile,regType)
            if nargin < 2, error('Reference head model is missing.');end
            if nargin < 3, individualHeadModelFile = ['surfaces_' num2str(round(1e5*rand)) '.mat'];end
            if nargin < 4, regType = 'bspline';end
            if isempty(obj.channelSpace) || isempty(obj.label) || isempty(obj.fiducials), error('Channel space or fiducials are missing.');end
            if ~exist(headModelFile,'file'), error('The file you''ve entered does not exist.');end
                       
            template = load(headModelFile);
            surfData = template.surfData;
            gTools = geometricTools;
            th = norminv(0.90);
            % mapping source to target spaces: S->T
            % target space: template
            T = [template.fiducials.nasion;...
                template.fiducials.lpa;...
                template.fiducials.rpa;...
                template.fiducials.vertex];
            
            % source space: individual geometry
            S = [obj.fiducials.nasion;...
                obj.fiducials.lpa;...
                obj.fiducials.rpa];
            
            % estimates vertex if is missing
            if isfield(obj.fiducials,'vertex')
                if numel(obj.fiducials.vertex) == 3
                    S = [S;obj.fiducials.vertex];
                else
                    point = 0.5*(obj.fiducials.lpa + obj.fiducials.rpa);
                    point = ones(50,1)*point;
                    point(:,3) = linspace(point(3),1.5*max(obj.channelSpace(:,3)),50)';
                    [~,d] = gTools.nearestNeighbor(obj.channelSpace,point);
                    [~,loc] = min(d);
                    point = point(loc,:);
                    S = [S;point];
                end
            else
                point = 0.5*(obj.fiducials.lpa + obj.fiducials.rpa);
                point = ones(50,1)*point;
                point(:,3) = linspace(point(3),1.5*max(obj.channelSpace(:,3)),50)';
                [~,d] = gTools.nearestNeighbor(obj.channelSpace,point);
                [~,loc] = min(d);
                point = point(loc,:);
                S = [S;point];
            end
            
            if isfield(obj.fiducials,'inion')
                if numel(obj.fiducials.vertex) == 3
                    S = [S;obj.fiducials.inion];
                    T = [T;template.fiducials.inion];
                end
            end
            if isa(obj,'eeg')
                obj.container.container.initStatusbar(1,8,'Co-registering...');
            else
                disp('Co-registering...');
            end
            
            % affine co-registration
            Aff = gTools.affineMapping(S,T);
            if isa(obj,'eeg'), obj.container.container.statusbar(1);end
            
            obj.channelSpace = gTools.applyAffineMapping(obj.channelSpace,Aff);
            obj.fiducials.lpa = gTools.applyAffineMapping(obj.fiducials.lpa,Aff);
            obj.fiducials.rpa = gTools.applyAffineMapping(obj.fiducials.rpa,Aff);
            obj.fiducials.nasion = gTools.applyAffineMapping(obj.fiducials.nasion,Aff);
            
            if ~strcmp(regType,'affine')
                % b-spline co-registration (only fiducial landmarks)
                options.Verbose = true;
                options.MaxRef = 2;
                Saff = gTools.applyAffineMapping(S,Aff);
                [Def,spacing,offset] = gTools.bSplineMapping(Saff,T,obj.channelSpace,options);
                if isa(obj,'eeg'), obj.container.container.statusbar(2);end
                
                % b-spline co-registration (second pass)
                obj.channelSpace = gTools.applyBSplineMapping(Def,spacing,offset,obj.channelSpace);
                obj.fiducials.lpa = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.lpa);
                obj.fiducials.rpa = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.rpa);
                obj.fiducials.nasion = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.nasion);
                
                T = template.surfData(1).vertices;
                S = obj.channelSpace;
                S(S(:,3) <= min(T(:,3)),:) = [];
                [S,d] = gTools.nearestNeighbor(S,T);
                z = zscore(d);
                S(abs(z)>th,:) = [];
                T(abs(z)>th,:) = [];
                [Def,spacing,offset] = gTools.bSplineMapping(S,T,obj.channelSpace,options);
                if isa(obj,'eeg'), obj.container.container.statusbar(3);end
                
                % b-spline co-registration (third pass)
                obj.channelSpace = gTools.applyBSplineMapping(Def,spacing,offset,obj.channelSpace);
                obj.fiducials.lpa = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.lpa);
                obj.fiducials.rpa = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.rpa);
                obj.fiducials.nasion = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.nasion);
                
                T = template.surfData(1).vertices;
                S = obj.channelSpace;
                S(S(:,3) <= min(T(:,3)),:) = [];
                [S,d] = gTools.nearestNeighbor(S,T);
                z = zscore(d);
                S(abs(z)>th,:) = [];
                T(abs(z)>th,:) = [];
                Tm = 0.5*(T+S);
                [Def,spacing,offset] = gTools.bSplineMapping(S,Tm,obj.channelSpace,options);
                if isa(obj,'eeg'), obj.container.container.statusbar(4);end
                
                % apply the final transformation
                obj.channelSpace = gTools.applyBSplineMapping(Def,spacing,offset,obj.channelSpace);
                obj.fiducials.lpa = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.lpa);
                obj.fiducials.rpa = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.rpa);
                obj.fiducials.nasion = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.nasion);
            end
            
            % fixing topological defects
            if isa(obj,'eeg')
                obj.container.container.statusBar.setText('Fixing topological defects...');
            else
                disp('Fixing topological defects...');
            end
            Ns = length(surfData);
            dmax = ones(Ns-1,1)*5;
            dmax(1) = 8;
            ind = fliplr(1:Ns);
            for it=1:Ns-1
                surfData(ind(it+1)).vertices = gTools.repareIntersectedSurface(surfData(ind(it)),surfData(ind(it+1)),dmax(it));
                if isa(obj,'eeg'), obj.container.container.statusbar(it+5);end
            end
            
            ind =  obj.channelSpace(:,3) > min(surfData(1).vertices(:,3));
            T = gTools.nearestNeighbor(surfData(1).vertices,obj.channelSpace);
            channelSpace = obj.channelSpace; %#ok
            channelSpace(ind,:) = T(ind,:);  %#ok
            [~,loc] = unique_bc(channelSpace,'rows');%#ok
            indInterp = setdiff_bc(1:size(obj.channelSpace,1),loc);
            if ~isempty(indInterp)
                x = setdiff_bc(channelSpace,channelSpace(indInterp,:),'rows');%#ok
                xi = gTools.nearestNeighbor(x,channelSpace(indInterp,:));%#ok
                channelSpace(indInterp,:) = 0.5*(xi + channelSpace(indInterp,:));%#ok
            end
            obj.channelSpace = channelSpace; %#ok
            
            if isfield(template,'atlas'), obj.atlas = template.atlas;end
            if exist(obj.surfaces,'file'), delete(obj.surfaces);end
            obj.surfaces = individualHeadModelFile;
            save(obj.surfaces,'surfData');
            if isa(obj,'eeg'), obj.container.container.statusbar(8);end
        end
        %%
        function computeLeadFieldBEM(obj, conductivity,normalityConstrained)
            % conductivity values taken from Valdes-Hernandez et al., 2009
            % also check Oostendrop TF, 2000; Wendel and Malmivuo, 2006
            %
            % brain and scalp: 0.33 S/m
            % skull: 0.022 S/m
            dispCommand = false;
            if nargin < 2, conductivity = [0.33 0.022 0.33];end
            if nargin < 3, normalityConstrained = true;end
            if isempty(obj.channelSpace) || isempty(obj.label) || isempty(obj.surfaces);
                error('Head model is incomplete or missing.');
            end
            if any(conductivity == -1)
                prefObj = [...
                    PropertyGridField('conductivity',[0.33 0.022 0.33],'DisplayName','Conductivity','Description',sprintf('Conductivity values are taken from Valdes-Hernandez et al., 2006, check \nalso Oostendrop TF, 2000; Wendel and Malmivuo, 2006. \nbrain and scalp: 0.33 S/m\nskull: 0.022 S/m'))...
                    PropertyGridField('normalityConstrained',true,'DisplayName','normalityConstrained','Description','If true, compute the LF matrix constrained to the normal vectors to the cortical surface resulting in a matris Nsensors X Nvertices. If false the LF matrix is Nsensors X 3*Nvertices')...
                    ];
                hFigure = figure('MenuBar','none','Name','OpenMEEG solver','NumberTitle', 'off','Toolbar', 'none','Units','pixels','Color',obj.container.container.preferences.gui.backgroundColor,...
                    'Resize','off','userData',0);
                position = get(hFigure,'position');
                set(hFigure,'position',[position(1:2) 303 250]);
                hPanel = uipanel(hFigure,'Title','','BackgroundColor','white','Units','pixels','Position',[0 55 303 175],'BorderType','none');
                g = PropertyGrid(hPanel,'Properties', prefObj,'Position', [0 0 1 1]);
                uicontrol(hFigure,'Position',[72 15 70 21],'String','Cancel','ForegroundColor',obj.container.container.preferences.gui.fontColor,...
                    'BackgroundColor',obj.container.container.preferences.gui.buttonColor,'Callback',@cancelCallback);
                uicontrol(hFigure,'Position',[164 15 70 21],'String','Ok','ForegroundColor',obj.container.container.preferences.gui.fontColor,...
                    'BackgroundColor',obj.container.container.preferences.gui.buttonColor,'Callback',@okCallback);
                uiwait(hFigure);
                if ~ishandle(hFigure), return;end
                if ~get(hFigure,'userData'), close(hFigure);return;end
                close(hFigure);
                drawnow;
                val = g.GetPropertyValues();
                conductivity = val.conductivity;
                normalityConstrained = val.normalityConstrained;
                dispCommand = true;
            end
            
            if dispCommand
                disp('Running:');
                if isa(obj,'coreStreamObject')
                    itemIndex = obj.container.findItem(obj.uuid);
                    fprintf('  mobilab.allStreams.item{%i}.computeLeadFieldBEM( [ %i %i %i ], %i );\n',itemIndex,conductivity(1),conductivity(2),conductivity(3),normalityConstrained);
                else
                    fprintf('  obj.computeLeadFieldBEM( [ %i %i %i ], %i );\n',conductivity(1),conductivity(2),conductivity(3),normalityConstrained);
                end
            end
            
            if ~exist(obj.surfaces,'file'), error('The file containing the surfaces is missing.');end
            [~,msg] = system('which om_assemble');
            existOM = ~isempty(strfind(msg,'Command not found'));
            if existOM
                try
                    mobilab = evalin('base','mobilab');
                    mobilabPath = mobilab.path;
                catch %#ok
                    mobilabPath = which('mobilabApplication');
                    if ~isempty(mobilabPath)
                        mobilabPath = fileparts(mobilabPath);
                    else
                        error('OpenMEEG is not intalled. Please download and install the sources you need from https://gforge.inria.fr/frs/?group_id=435.');
                    end
                end
                openmeegDir = [mobilabPath filesep 'dependency' filesep 'openmeeg'];
                
                %---
                % Approach taken from Brainstorm's function bst_openmeeg,
                % Francois Tadel & Alexandre Gramfort, 2011
                %---
                if ~ispc
                    if ismember_bc(computer, {'GLNX86','GLNXA64'})
                        varname = 'LD_LIBRARY_PATH';
                    else
                        varname = 'DYLD_LIBRARY_PATH';
                    end
                    libpath = getenv(varname);
                    if ~isempty(libpath), libpath = [libpath ':'];end
                    if isempty(strfind(lower(libpath),'openmeeg')), setenv(varname, [libpath openmeegDir]);end
                end
                % Set number of cores used
                try
                    numcores = feature('numcores');
                catch %#ok
                    numcores = 4;
                end
                setenv('OMP_NUM_THREADS', num2str(numcores));
                %---
            end
            
            load(obj.surfaces);
            Ns = length(surfData); %#ok
            gTools = geometricTools;
            
            rootDir = fileparts(obj.surfaces);
            if isempty(rootDir), rootDir = pwd;end
            if ispc
                binDir = fileparts(which('om_assemble.exe'));
            else
                binDir = fileparts(which('om_assemble'));
            end
            headModelGeometry = fullfile(rootDir,'head_model.geom');
            copyfile( which('head_model.geom'),headModelGeometry,'f');
            c1 = onCleanup(@()delete(headModelGeometry));
            
            headModelConductivity = fullfile(rootDir,'head_model.cond');
            fid = fopen(headModelConductivity,'w');
            fprintf(fid,'# Properties Description 1.0 (Conductivities)\n\nAir         0.0\nScalp       %.3f\nBrain       %0.3f\nSkull       %0.3f',...
                conductivity(1),conductivity(3),conductivity(2));
            fclose(fid);
            c2 = onCleanup(@()delete(headModelConductivity));
            
            dipolesFile = fullfile(rootDir,'cortex_dipoles.txt');
            normalsIn = false;
            [normals,surfData(Ns).faces] = gTools.getSurfaceNormals(surfData(Ns).vertices,surfData(Ns).faces,normalsIn);
            
            if normalityConstrained
                sourceSpace = [surfData(Ns).vertices normals];
            else
                One = ones(length(normals(:,2)),1);
                Zero = 0*One;
                sourceSpace = [surfData(Ns).vertices One Zero Zero;...
                    surfData(Ns).vertices Zero One Zero;...
                    surfData(Ns).vertices Zero Zero One];
            end
            dlmwrite(dipolesFile, sourceSpace, 'precision', 6,'delimiter',' ')
            c3 = onCleanup(@()delete(dipolesFile));
            
            electrodesFile = fullfile(rootDir,'eeg_channels_locations.txt');
            dlmwrite(electrodesFile, obj.channelSpace, 'precision', 6,'delimiter',' ')
            c4 = onCleanup(@()delete(electrodesFile));
            
            normalsIn = true;
            brain = fullfile(rootDir,'brain.tri');
            if Ns == 4
                [normals,surfData(3).faces] = gTools.getSurfaceNormals(surfData(3).vertices,surfData(3).faces,normalsIn);
                om_save_tri(brain,surfData(3).vertices,surfData(3).faces,normals)
            else
                [normals,surfData(2).faces] = gTools.getSurfaceNormals(surfData(2).vertices,surfData(2).faces,normalsIn);
                fakeSurf = surfData(2);
                fakeSurf.vertices = surfData(2).vertices + 1.05*normals;
                surfData(2).vertices = surfData(2).vertices - 1.05*normals;
                fakeSurf.vertices = gTools.repareIntersectedSurface(surfData(3),fakeSurf,3);
                [normals,fakeSurf.faces] = gTools.getSurfaceNormals(fakeSurf.vertices,fakeSurf.faces,normalsIn);
                om_save_tri(brain,fakeSurf.vertices,fakeSurf.faces,normals)
            end
            c5 = onCleanup(@()delete(brain));
            
            skull = fullfile(rootDir,'skull.tri');
            [normals,surfData(2).faces] = gTools.getSurfaceNormals(surfData(2).vertices,surfData(2).faces,normalsIn);
            om_save_tri(skull,surfData(2).vertices,surfData(2).faces,normals)
            c6 = onCleanup(@()delete(skull));
            
            head = fullfile(rootDir,'head.tri');
            [normals,surfData(1).faces] = gTools.getSurfaceNormals(surfData(1).vertices,surfData(1).faces,normalsIn);
            om_save_tri(head,surfData(1).vertices,surfData(1).faces,normals)
            c7 = onCleanup(@()delete(head));
            
            hmFile = fullfile(rootDir,'hm.bin');       c8 = onCleanup(@()delete(hmFile));
            hmInvFile = fullfile(rootDir,'hm_inv.bin');c9 = onCleanup(@()delete(hmInvFile));
            dsmFile = fullfile(rootDir,'dsm.bin');     c10 = onCleanup(@()delete(dsmFile));
            h2emFile = fullfile(rootDir,'h2em.bin');   c11 = onCleanup(@()delete(h2emFile));
            lfFile = fullfile(rootDir,'leadfield.mat');
            
            hmFile = strrep(hmFile,'\','/');
            hmInvFile = strrep(hmInvFile,'\','/');
            dsmFile = strrep(dsmFile,'\','/');
            h2emFile = strrep(h2emFile,'\','/');
            headModelGeometry = strrep(headModelGeometry,'\','/');
            
            wDir = pwd;
            cd(binDir);
            
            try
                if ispc
                    system(['om_assemble -HM "' headModelGeometry '" "' headModelConductivity '" "' hmFile '"']);
                    system(['om_minverser "' hmFile '" "' hmInvFile '"']);
                    system(['om_assemble -DSM "' headModelGeometry '" "' headModelConductivity '" "' dipolesFile '" "' dsmFile '"']);
                    system(['om_assemble -H2EM "' headModelGeometry '" "' headModelConductivity '" "' electrodesFile '" "' h2emFile '"']);
                    system(['om_gain -EEG "' hmInvFile '" "' dsmFile '" "' h2emFile '" "' lfFile '"']);
                else
                    system(['./om_assemble -HM "' headModelGeometry '" "' headModelConductivity '" "' hmFile '"']);
                    system(['./om_minverser "' hmFile '" "' hmInvFile '"']);
                    system(['./om_assemble -DSM "' headModelGeometry '" "' headModelConductivity '" "' dipolesFile '" "' dsmFile '"']);
                    system(['./om_assemble -H2EM "' headModelGeometry '" "' headModelConductivity '" "' electrodesFile '" "' h2emFile '"']);
                    system(['./om_gain -EEG "' hmInvFile '" "' dsmFile '" "' h2emFile '" "' lfFile '"']);
                end
            catch ME
                cd(wDir);
                error(ME.message);
            end
            if ~exist(lfFile,'file')
                fid = fopen('compute_leadfield.bat','w+'); c12 = onCleanup(@()delete('compute_leadfield.bat'));
                fprintf(fid,'om_assemble -HeadMat "%s" "%s" "%s"\n',headModelGeometry,headModelConductivity,hmFile);
                fprintf(fid,'om_minverser "%s" "%s"\n',hmFile,hmInvFile);
                fprintf(fid,'om_assemble -DSM "%s" "%s" "%s" "%s"\n',headModelGeometry,headModelConductivity,dipolesFile,dsmFile);
                fprintf(fid,'om_assemble -H2EM "%s" "%s" "%s" "%s"\n',headModelGeometry,headModelConductivity,electrodesFile,h2emFile);
                fprintf(fid,'om_gain -EEG "%s" "%s" "%s" "%s"\n',hmInvFile,dsmFile,h2emFile,strrep(lfFile,'\','/'));
                fclose(fid);
                
                system('compute_leadfield.bat');
            end
                
            if ~exist(lfFile,'file'), error('Cannot find OpenMEEG libraries. You can install OpenMEEG from Unix terminal with the following command: root@machine# apt-get install libopenmeeg1');end
            
            load(lfFile);
            cd(wDir);
            
            K = linop;
            clear linop;
            z = zscore(K(:));
            ind = find(z<norminv(0.01) | z>norminv(0.99));
            % removing extreme values due to numerical instability
            K(ind) = 0; %#ok
            
            if exist(lfFile,'file'), delete(lfFile);end
            if isa(obj,'coreStreamObject')
                lfFile = fullfile(obj.container.mobiDataDirectory,['lf_' obj.name '_' char(obj.uuid) '.mat']);
            end
            obj.leadFieldFile = lfFile;
            save(obj.leadFieldFile,'K');
            disp('Done.')
        end
        %%
        function [sourceSpace,K,L,rmIndices] = getSourceSpace4PEB(obj,structName,rmIndices)
            
            if isempty(obj.atlas) || isempty(obj.surfaces) || isempty(obj.leadFieldFile), error('Head model, leadfield, or atlas are missing.');end
            if nargin < 2
                structName = 'Thalamus';
                disp('Undefined structure to remove. Opening the surface by the Thalamus.')
            end
            if nargin < 3
                rmIndices = [];
            end
            load(obj.surfaces);
            sourceSpace = surfData(end); %#ok
            load(obj.leadFieldFile);
            if ~exist('L','var'),
                L = geometricTools.getSurfaceLaplacian(sourceSpace.vertices,sourceSpace.faces);
                save(obj.leadFieldFile,'K','L','-mat')
            end
            
            try
                [sourceSpace,rmIndices] = obj.removeStructureFromSourceSpace(structName,rmIndices);
            catch ME
                if strcmp(ME.identifier,'MoBILAB:noStructureMatched')
                    I = strfind(obj.atlas.label,'Thalamus');
                    ind = ~cellfun(@isempty,I);
%                     ind = false(length(obj.atlas.label),1);
%                     for it=1:length(obj.atlas.label), if ind(isempty(I{it})) = true;end
                    if any(ind)
                        indices = indices4Structure(obj,obj.atlas.label(ind));
                        mxInd = length(indices);
                        mxInd(mxInd>40) = 40;
                        rmIndices = [indices(1:mxInd) rmIndices(:)'];
                    else
                        n = size(sourceSpace.vertices,1);
                        rmIndices = [fix(n/2)-40:fix(n/2)+39 rmIndices(:)'];
                    end
                else
                  n = size(sourceSpace.vertices,1);
                  rmIndices = [fix(n/2)-40:fix(n/2)+39 rmIndices(:)'];
                end
                if ~isempty(rmIndices)
                    [nVertices,nFaces] = geometricTools.openSurface(sourceSpace.vertices,sourceSpace.faces,rmIndices);
                    sourceSpace.vertices = nVertices;
                    sourceSpace.faces = nFaces;
                end
            end
            dim = size(K); %#ok
            L(rmIndices,:) = [];
            L(:,rmIndices) = [];
            if dim(2)/3 == size(surfData(end).vertices,1) %#ok
                K = reshape(K,[dim(1) dim(2)/3 3]);
                K(:,rmIndices,:) = [];
                % K = permute(K,[1 3 2]);
                K = reshape(K,[dim(1) (dim(2)/3-length(rmIndices))*3]);
                L = kron(eye(3),L);
            else
                K(:,rmIndices) = [];
            end
        end
        %%
        function indices = indices4Structure(obj,structName)            
            ind = find(ismember_bc(obj.atlas.label,structName));
            if isempty(ind), error('MoBILAB:noStructureMatched','The structure you want to remove is not defined in this atlas.');end
            tmp = bsxfun(@eq,obj.atlas.color,ind');
            indices = find(any(tmp,2));
        end
            
        function [sourceSpace,rmIndices] = removeStructureFromSourceSpace(obj,structName,structIndices)
            if nargin<3
                structIndices = [];
            end
          
            if ~iscell(structName), structName = {structName}; end
            
            if isempty(obj.atlas) || isempty(obj.surfaces), error('Head model or atlas are missing.');end
            if nargin < 2, error('Not enough input arguments.');end
            
            load(obj.surfaces);
            sourceSpace = surfData(end);%#ok
            
            tmpIndices = indices4Structure(obj,structName);
            if isempty(tmpIndices) && isempty(structIndices)
                error('MoBILAB:noStructureMatched','The structure you want to remove is not defined in this atlas.');
            end
            
            if ~isempty(structIndices)
                % concatenate elements of structIndices into a single column vector
                structIndices = cellfun(@(x)x(:),structIndices,'UniformOutput',false)';
                structIndices = cell2mat(structIndices);
            end
            rmIndices = unique_bc([tmpIndices(:) ; structIndices]);
            
            [nVertices,nFaces] = geometricTools.openSurface(sourceSpace.vertices,sourceSpace.faces,rmIndices);
            sourceSpace.vertices = nVertices;
            sourceSpace.faces = nFaces;
        end
        %%
        % call NFT routines (doesn't work)
        function computeLeadFieldBEM_NFT(obj, conductivity, normalityConstrained)
            % conductivity values taken from Valdes-Hernandez et al., 2006
            % also check Oostendrop TF, 2000; Wendel and Malmivuo, 2006
            %
            % brain and scalp: 0.33 S/m
            % skull: 0.022 S/m
            if nargin < 2, conductivity = [0.33 0.022 0.33];end
            if isempty(obj.channelSpace) || isempty(obj.label) || isempty(obj.surfaces);
                error('Head model is incomplete or missing.');
            end
            if any(conductivity == -1)
                prompt = {'Enter conductivities (in S/m) for each layer of tissue: scalp, skull, and brain'};
                dlg_title = 'NFT BEM solver';
                num_lines = 1;
                def =  {'0.33 0.022 0.33'};
                varargin = inputdlg2(prompt,dlg_title,num_lines,def);
                if isempty(varargin), return;end
                conductivity = eval(['[' varargin{1} '];']);
                if isnan(conductivity), return;end
            end
            
            if ~exist(obj.surfaces,'file'), error('The file containing the surfaces is missing.');end
            if isempty(which('utilbem_add_mesh')), error('NFT plugin is missing, You can download it from: http://sccn.ucsd.edu/nft/');end
            
            
            session_name = 'mesh';
            if isa(obj,'eeg')
                of = obj.container.mobiDataDirectory;
                subject_name = ['BEM_' obj.name];
            else
                of = pwd;
                subject_name = 'BEM_surfaces.bec';
            end
            
            ssave.fn = '';
            ssave.eloc = 0;
            ssave.pnt = obj.channelSpace;
            ssave.ind = 1:size(obj.channelSpace,1);%#ok
            sensors_file = fullfile(of,[subject_name '_' session_name '.sensors']);
            save(sensors_file, '-STRUCT', 'ssave');
            
            
            load(obj.surfaces);
            
            nv = size(surfData(1).vertices,1);
            nf = size(surfData(1).faces,1);
            Cscalp = [(1:nv)' surfData(1).vertices];
            Escalp = [(1:nf)' surfData(1).faces];
            
            nv = size(surfData(2).vertices,1);
            nf = size(surfData(2).faces,1);
            Cskull = [(1:nv)' surfData(2).vertices];
            Eskull = [(1:nf)' surfData(2).faces];
            
            nv = size(surfData(3).vertices,1);
            nf = size(surfData(3).faces,1);
            Cbrain = [(1:nv)' surfData(3).vertices];
            Ebrain = [(1:nf)' surfData(3).faces];
            
            [Coord, Elem] = utilbem_add_mesh(Cscalp,Escalp,Cskull,Eskull);
            [Coord, Elem] = utilbem_add_mesh(Coord,Elem,Cbrain,Ebrain);%#ok
            
            bec_file = fullfile(of,[subject_name '.bec']);
            save(bec_file, 'Coord', '-ascii');
            
            nl = length(surfData);
            Info(1,1) = nl;
            Info(1,2) = size(Escalp,1)+size(Eskull,1)+size(Ebrain,1);
            Info(1,3) = size(Cscalp,1)+size(Cskull,1)+size(Cbrain,1);
            
            Info(1,4) = size(Escalp,2)-1; % number of nodes per element
            Info(2,1) = 1;
            Info(2,2) = size(Escalp,1);
            Info(2,3:4) = [1 0];
            Info(3,1) = 2;
            Info(3,2) = size(Eskull,1);
            Info(3,3:4) = [2 1];
            Info(4,1) = 3;
            Info(4,2) = size(Ebrain,1);
            Info(4,3:4) = [3 2];
            
            bei_file = fullfile(of,[subject_name '.bei']);
            fid = fopen(bei_file,'w');
            fprintf(fid, '%d %d %d %d\r\n', Info');
            fclose(fid);
            
            bee_file = fullfile(of,[subject_name '.bee']);
            fid = fopen(bee_file,'w');
            if size(Elem,2) == 4
                fprintf(fid, '%d %d %d %d\r\n', Elem');
            elseif size(Elem,2) ==7
                fprintf(fid, '%d %d %d %d %d %d %d\r\n', Elem');
            end
            fclose(fid);
            
            Ns = size(surfData(3).vertices,1);
%             so = zeros(Ns, 6);
%             so(:, 1:3) = surfData(3).vertices;
%             so(:, 4:6) = geometricTools.getSurfaceNormals(surfData(3).vertices,surfData(3).faces,false); %#ok
            
            so = zeros(3*Ns, 6);
            so(1:Ns, 1:3) = surfData(3).vertices;
            so(1:Ns, 4) = 1;
            so(1+Ns:2*Ns, 1:3) = surfData(3).vertices;
            so(1+Ns:2*Ns, 5) = 1;
            so(1+Ns*2:3*Ns, 1:3) = surfData(3).vertices;
            so(1+Ns*2:3*Ns, 6) = 1; %#ok
            %
            % save source space
            dip_file = fullfile(of,[subject_name '_sourcespace.dip']);
            save(dip_file, 'so', '-ascii');
            
            try
                nft_wrapper(subject_name, session_name, of, 'cond', conductivity)
                % nft_forward_problem_solution(subject_name, session_name, of, 'cond', conductivity);
            catch ME
                % cleaning up
                delete(sensors_file);
                delete(bec_file);
                delete(bei_file);
                delete(bee_file);
                delete(dip_file);
                ME.rethrow;
            end
            
            % cleaning up
            delete(sensors_file);
            delete(bec_file);
            delete(bei_file);
            delete(bee_file);
            delete(dip_file);
            lfFile = fullfile(of,[session_name '_LFM.mat']);
            
            load(lfFile);
            cd(wDir);
            
            K = LFM;
            clear LFM;
            z = zscore(K(:));
            ind = find(z<norminv(0.01) | z>norminv(0.99));
            % removing extreme values due to numerical instability
            K(ind) = 0; %#ok
            
            if exist(lfFile,'file'), delete(lfFile);end
            if isa(obj,'coreStreamObject')
                lfFile = fullfile(obj.container.mobiDataDirectory,['lf_' obj.name '_' char(obj.uuid) '.mat']);
            end
            obj.leadFieldFile = lfFile;
            save(obj.leadFieldFile,'K');
        end
        %         %%
        %         function computeLeadFieldFEM(obj)
        %
        %         end
        function saveToFile(obj,filename)
            if exist(obj.surfaces,'file')
                load(obj.surfaces,'-mat');
            else
                surfData = [];
            end
            metadata.surfaces = surfData;
            if ~isempty(obj.surfaces)
                [~,metadata.surfacesFilename,ext] = fileparts(obj.surfaces);
                metadata.surfacesFilename = [metadata.surfacesFilename ext];
            else
                metadata.surfacesFilename = '';
            end
            metadata.atlas = obj.atlas;
            metadata.fiducials = obj.fiducials;
            metadata.channelSpace = obj.channelSpace;
            metadata.label = obj.label;
            
            if exist(obj.leadFieldFile,'file')
                load(obj.leadFieldFile,'-mat');
            else
                K = [];
            end
            metadata.leadField = K;
            if exist('L','var'), metadata.L = L;end
            if ~isempty(obj.leadFieldFile);
                [~,metadata.leadFieldFile,ext] = fileparts(obj.leadFieldFile);
                metadata.leadFieldFile = [metadata.leadFieldFile ext];
            else
                metadata.leadFieldFile = '';
            end
            if isprop(obj,'surfNormal')
                metadata.surfNormal = obj.surfNormal;
            end
            save(filename,'metadata','-mat');
        end
        
        function metadata = saveobj(obj)
            error('Please use saveToFile() to save this object');
        end
        
    end
    methods(Static)
        function obj = loadFromFile(filename)
            
            if isempty(dir(filename));
                error('File %s does not exist ',filename);
            end
                
            filePath = fileparts(filename);
            
            % if filePath is the name of a file in current directory...
            if isempty(filePath)
                % ... use curdir as filePath
                filePath = pwd;
            end
            
            load(filename,'-mat');
            surfData = metadata.surfaces; %#ok
            if ~isempty(metadata.leadFieldFile)
                % update folder pointers for data structures 
                % leadfield
                [~, leadFieldFile, ext] = fileparts(metadata.leadFieldFile);
                leadFieldFile = [leadFieldFile ext];
                leadFieldFile = [filePath filesep leadFieldFile];
                metadata.leadFieldFile = leadFieldFile;
            else
                metadata.leadFieldFile = '';
            end
            % surfaces 
            [~, surfacesFilename, ext] = fileparts(metadata.surfacesFilename);
            surfacesFilename = [surfacesFilename ext];
            surfacesFilename = [filePath filesep surfacesFilename];
            metadata.surfacesFilename = surfacesFilename;
            
            save(metadata.surfacesFilename,'-mat','surfData');
            K = metadata.leadField; %#ok
            if ~isempty(metadata.leadFieldFile)
                if isfield(metadata,'L')
                    L = metadata.L;
                    save(metadata.leadFieldFile,'-mat','K','L');
                else
                    save(metadata.leadFieldFile,'-mat','K');
                end
            end
            if isfield(metadata,'surfNormal')
                surfNormal = metadata.surfNormal;
            else
                surfNormal = [];
            end
            obj = headModel('surfaces',metadata.surfacesFilename,'atlas',metadata.atlas,'fiducials',metadata.fiducials,...
                'channelSpace',metadata.channelSpace,'label',metadata.label,'leadFieldFile',[metadata.leadFieldFile],'surfNormal',surfNormal);
        end
        
        function obj = loadobj(metadata)
            error('Please use loadFromFile() to load this object');
        end
    end
end

%--
function [elec,labels,fiducials] = readMontage(file)
[eloc, labels] = readlocs(file);
elec = [cell2mat({eloc.X}'), cell2mat({eloc.Y}'), cell2mat({eloc.Z}')];
Nl = length(labels);
count = 1;
for it=1:Nl
    if ~isempty(strfind(labels{it},'fidnz')) || ~isempty(strfind(labels{it},'nasion'))
        fiducials.nasion = elec(it,:);
        count = count+1;
    elseif ~isempty(strfind(labels{it},'fidt9')) || ~isempty(strfind(labels{it},'lpa'))
        fiducials.lpa = elec(it,:);  
        count = count+1;
    elseif ~isempty(strfind(labels{it},'fidt10')) || ~isempty(strfind(labels{it},'rpa'))
        fiducials.rpa = elec(it,:);
        count = count+1;
    elseif ~isempty(strfind(labels{it},'fidt10')) || ~isempty(strfind(labels{it},'vertex'))
        fiducials.vertex = elec(it,:);
        count = count+1;
    end
    if count > 4, break;end
end
end