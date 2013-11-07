classdef currentSourceViewer < handle
    properties
        hFigure
        hAxes
        streamObj
        hLabels
        hSensors
        hScalp
        hCortex
        hVector
        dcmHandle
        scalpInterpMat
        hTimeLabel
        normals
    end
    methods
        function obj = currentSourceViewer(streamObj,J,V,figureTitle,channelLabels,objIn,surfData,colorLimits,currentTime)
            
            persistent oldColorLimits
            
            c = onCleanup(@() clear('currentSourceViewer'));
            
            if exist('objIn','var') && isprop(objIn,'hFigure')
                obj = objIn;
            else
                objIn = [];
                obj.streamObj = streamObj;
            end
            
            if nargin < 3, V = [];end
            if nargin < 4, figureTitle = '';end
            if nargin < 5
                channelLabels = cell(length(V),1);
                for it=1:length(channelLabels), channelLabels{it} = num2str(it);end
            end
            if ~exist('colorLimits','var')
                colorLimits = [];
            end
            if ~exist('surfData','var')
                load(obj.streamObj.surfaces);
            end
            if ~exist('currentTime','var')
                currentTime = '';
            end
            
            if ~isempty(objIn);
                % update the figure only
                
                % update normal vectors
                if ~exist('surfData','var')
                    load(obj.streamObj.surfaces);
                end
                [~, Jm] = drawNormalVectors(J,surfData(3),obj.hVector,obj.normals);
                
                % update cortex colors
                set(obj.hCortex,'FaceVertexCData',Jm);
            
                % update scalp colors
                if ~isempty(V)
                    %Vi = geometricTools.splineInterpolator(obj.streamObj.channelSpace,V,surfData(1).vertices);
                    if isempty(obj.scalpInterpMat)
                        obj.scalpInterpMat = geometricTools.localGaussianInterpolator(obj.streamObj.channelSpace,surfData(1).vertices,6);
                        obj.scalpInterpMat = obj.scalpInterpMat/max(obj.scalpInterpMat(:));
                    end
                    set(obj.hScalp,'FaceVertexCData',obj.scalpInterpMat*double(V));
                end
                
                if ~isempty(currentTime)
                    % update time label
                    set(obj.hTimeLabel,'string',sprintf('%5.5g sec\n',currentTime));
                end
                
                if isempty(colorLimits)
                    mx = max(Jm(:));
                    set(obj.hAxes,'Clim',[-mx mx]);
                elseif ~isequal(colorLimits,oldColorLimits)
                    set(obj.hAxes,'Clim',colorLimits);
                    oldColorLimits = colorLimits;
                end
            
                return;
            end
                
            % create a new figure
            color = [0.93 0.96 1];

            mPath = which('mobilabApplication');
            path = fullfile(fileparts(mPath),'skin');
            labelsOn  = imread([path filesep 'labelsOn.png']);
            labelsOff = imread([path filesep 'labelsOff.png']);
            sensorsOn = imread([path filesep 'sensorsOn.png']);
            sensorsOff = imread([path filesep 'sensorsOff.png']);
            scalpOn = imread([path filesep 'scalpOn.png']);
            scalpOff = imread([path filesep 'scalpOff.png']);
            vectorOn = imread([path filesep 'vectorOn.png']);
            vectorOff = imread([path filesep 'vectorOff.png']);

            if isa(streamObj,'struct'), visible = 'off';else visible = 'on';end
            obj.hFigure = figure('Menubar','figure','ToolBar','figure','renderer','opengl','Visible',visible,'Color',color);
            obj.hAxes = axes('Parent',obj.hFigure);

            toolbarHandle = findall(obj.hFigure,'Type','uitoolbar');

            hcb(1) = uitoggletool(toolbarHandle,'CData',labelsOff,'Separator','on','HandleVisibility','off','TooltipString','Labels On/Off','userData',{labelsOn,labelsOff},'State','off');
            set(hcb(1),'OnCallback',@(src,event)rePaint(obj,hcb(1),'labelsOn'),'OffCallback',@(src, event)rePaint(obj,hcb(1),'labelsOff'));

            hcb(2) = uitoggletool(toolbarHandle,'CData',sensorsOff,'Separator','off','HandleVisibility','off','TooltipString','Sensors On/Off','userData',{sensorsOn,sensorsOff},'State','off');
            set(hcb(2),'OnCallback',@(src,event)rePaint(obj,hcb(2),'sensorsOn'),'OffCallback',@(src, event)rePaint(obj,hcb(2),'sensorsOff'));

            hcb(3) = uitoggletool(toolbarHandle,'CData',scalpOff,'Separator','off','HandleVisibility','off','TooltipString','Scalp On/Off','userData',{scalpOn,scalpOff},'State','off');
            set(hcb(3),'OnCallback',@(src,event)rePaint(obj,hcb(3),'scalpOn'),'OffCallback',@(src, event)rePaint(obj,hcb(3),'scalpOff'));

            hcb(4) = uitoggletool(toolbarHandle,'CData',vectorOff,'Separator','off','HandleVisibility','off','TooltipString','Scalp On/Off','userData',{vectorOn,vectorOff},'State','off');
            set(hcb(4),'OnCallback',@(src,event)rePaint(obj,hcb(4),'vectorOn'),'OffCallback',@(src, event)rePaint(obj,hcb(4),'vectorOff'));

            obj.dcmHandle = datacursormode(obj.hFigure);
            obj.dcmHandle.SnapToDataVertex = 'off';
            set(obj.dcmHandle,'UpdateFcn',@(src,event)showLabel(obj,event));
            obj.dcmHandle.Enable = 'off';

            hold(obj.hAxes,'on');

            obj.hSensors = scatter3(obj.hAxes,obj.streamObj.channelSpace(:,1),obj.streamObj.channelSpace(:,2),...
                obj.streamObj.channelSpace(:,3),'filled','MarkerEdgeColor','k','MarkerFaceColor','y');
            set(obj.hSensors,'Visible','off');

            N = length(channelLabels);
            k = 1.1;
            obj.hLabels = zeros(N,1);
            for it=1:N, obj.hLabels(it) = text('Position',k*obj.streamObj.channelSpace(it,:),'String',channelLabels{it},'Parent',obj.hAxes);end
            set(obj.hLabels,'Visible','off');
                        
            % vectors
            [obj.hVector Jm obj.normals] = drawNormalVectors(J,surfData(3),obj.hVector,obj.normals);
            set(obj.hVector,'Color','k','Visible','off');
            
            % cortex
            obj.hCortex = patch('vertices',surfData(3).vertices,'faces',surfData(3).faces,'FaceVertexCData',Jm,...
                'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
                'SpecularExponent',50,'SpecularStrength',0.5,'Parent',obj.hAxes);
            h=camlight(0,180);
            set(h,'Parent',obj.hAxes);
            h=camlight(0,0);
            set(h,'Parent',obj.hAxes);
            
            % scalp
            if isempty(V)
                skinColor = [1,.75,.65];
                obj.hScalp = patch('vertices',surfData(1).vertices,'faces',surfData(1).faces,'facecolor',skinColor,...
                    'facelighting','phong','LineStyle','none','FaceAlpha',.85,'Parent',obj.hAxes,'Visible','off');
            else
                %Vi = geometricTools.splineInterpolator(obj.streamObj.channelSpace,V,surfData(1).vertices);
                if isempty(obj.scalpInterpMat)
                    obj.scalpInterpMat = geometricTools.localGaussianInterpolator(obj.streamObj.channelSpace,surfData(1).vertices,6);
                    obj.scalpInterpMat = obj.scalpInterpMat/max(obj.scalpInterpMat(:));
                end
                Vi = obj.scalpInterpMat*double(V);
                obj.hScalp = patch('vertices',surfData(1).vertices,'faces',surfData(1).faces,'FaceVertexCData',Vi,...
                    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',0.85,'SpecularColorReflectance',0,...
                    'SpecularExponent',50,'SpecularStrength',0.5,'Parent',obj.hAxes,'Visible','off');
            end
            
            if ~isempty(currentTime)
                % create time lael
                obj.hTimeLabel = annotation('textbox',[0.0189, 0, 0.1, 0.1],'string',sprintf('%5.5g sec\n',currentTime),'edgecolor','none');
            end
                
            view(obj.hAxes,[90 0]);
            colorbar('peer',obj.hAxes);
            % box on;
            title(texlabel(figureTitle,'literal'));
            hold(obj.hAxes,'off');
            axis(obj.hAxes,'equal');
            axis(obj.hAxes,'off')
            if isempty(colorLimits)
                mx = max(Jm(:));
                set(obj.hAxes,'Clim',[-mx mx]);
            else
                set(obj.hAxes,'Clim',colorLimits);
            end
            set(obj.hFigure,'Visible',visible,'userData',obj);
            rotate3d(obj.hFigure);
            axis(obj.hAxes,'vis3d');
            
%             function hContour = drawRoiContour(obj,
                
            function [hVector Jm normals] = drawNormalVectors(J,surfData,hVector,normals)
                % calculate and plot (or update) normal vectors
                % if hVector is supplied and is a handle to a vector field 
                % then the existing vector field is updated.                
                if ~exist('hVector','var')
                    hVector = [];
                end
                if ~exist('normals','var') || isempty(normals)
                    normalsIn = false;
                    normals = geometricTools.getSurfaceNormals(surfData.vertices,surfData.faces,normalsIn);
                end
%                 normals = hlp_microcache('getnormals',@geometricTools.getSurfaceNormals,surfData.vertices,surfData.faces,normalsIn);
                if size(J,1) == 3*size(surfData.vertices,1)  
                    J = reshape(J,[size(J,1)/3 3]);
                    Jm = sqrt(sum(J.^2,2));
                    s = sign(dot(normals,J,2));
                    Jm = s.*Jm;
                    hVector = drawQuiverPlot(surfData.vertices,J,hVector);
                else
                    Jm = J;
                    hVector = drawQuiverPlot(surfData.vertices,normals,hVector);
                end
            end
            
            function hVector = drawQuiverPlot(startpoints,endpoints,hVector)
                % if hVector is supplied and is a handle to a vector field 
                % then the existing vector field is updated.
                if nargin < 3 || isempty(hVector) || ~ishandle(hVector)
                    % create new vector field
                    hVector = quiver3(startpoints(:,1),startpoints(:,2),startpoints(:,3), ...
                                      endpoints(:,1)  ,endpoints(:,2)  ,endpoints(:,3) ,2);
                else
                    % update current vector field
                    set(hVector,'Xdata',startpoints(:,1),'YData',startpoints(:,2),'ZData',startpoints(:,3), ...
                                'UData',endpoints(:,1),'VData',endpoints(:,2),'WData',endpoints(:,3));
                end
            end
                
            
        end
        %%
        function rePaint(obj,hObject,opt)
            CData = get(hObject,'userData');
            if isempty(strfind(opt,'Off'))
                set(hObject,'CData',CData{2});
            else
                set(hObject,'CData',CData{1});
            end
            switch opt
                case 'labelsOn'
                    set(obj.hLabels,'Visible','on');
                case 'labelsOff'
                    set(obj.hLabels,'Visible','off');
                case 'sensorsOn'
                    set(obj.hSensors,'Visible','on');
                case 'sensorsOff'
                    set(obj.hSensors,'Visible','off');
                case 'scalpOn'
                    val = get(obj.hScalp,'FaceVertexCData');
                    mx = max(abs([min(val) max(val)]));
                    %set(obj.hCortex,'Visible','off');
                    set(obj.hScalp,'Visible','on');
                    set(get(obj.hScalp,'Parent'),'Clim',[-mx mx]);
                case 'scalpOff'
                    val = get(obj.hCortex,'FaceVertexCData');
                    mx = max(val(:));
                    set(obj.hScalp,'Visible','off');
                    set(obj.hCortex,'Visible','on');
                    set(get(obj.hCortex,'Parent'),'Clim',[-mx mx]);
                case 'vectorOn'
                    set(obj.hVector,'Visible','on');
                    set(obj.hCortex,'FaceAlpha',0.75);
                case 'vectorOff'
                    set(obj.hVector,'Visible','off');
                    set(obj.hCortex,'FaceAlpha',1);
            end
        end
        %%
        function output_txt = showLabel(obj,event_obj)
            persistent DT
            if strcmp(obj.dcmHandle.Enable,'off'),return;end
            if isempty(DT)
                load(obj.streamObj.surfaces);
                vertices = surfData(3).vertices;
                DT = DelaunayTri(vertices(:,1),vertices(:,2),vertices(:,3));
            end
            pos = get(event_obj,'Position');
            loc = nearestNeighbor(DT, pos);
            output_txt = obj.streamObj.atlas.label{obj.streamObj.atlas.color(loc)};
            %updateCursor(obj.dcmHandle,pos);
        end
    end
end