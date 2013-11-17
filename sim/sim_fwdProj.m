function [scalpData srcData LFM chaninds centroids_LFM_idx roiVertices sourceCoords] = sim_fwdProj(varargin)
% Project source amplitudes through a forward model to generate channel data.
%
% Author: Tim Mullen, SCCN/INC/UCSD, 2013

headmodel_default = 'resources:/headmodels/standard-Colin27-385ch.mat';
if ~onl_isonline
    hmObj = arg_extract(varargin,{'hmObj','HeadModelObject'},[],headmodel_default);
    hmObj = hlp_validateHeadModelObject(hmObj);
    if isempty(hmObj)
        ROINames = {};
        defChlbl = {};
    else
        % get the unique ROI names
        [tmp, idx] = unique_bc(hmObj.atlas.label);
        ROINames = hmObj.atlas.label(sort(idx))';
        defChlbl = hmObj.getChannelLabels();
    end
else
    ROINames = {};
    defChlbl = {};
end

g=arg_define([0 Inf],varargin,...
    arg_norep({'sourceAmps','SourceAmplitudes'},mandatory,[],'Current density. This is a [num_sources x num_points] matrix of current density for each source'), ...
    arg({'hmObj','HeadModelObject'},headmodel_default,[],'Head model object generated by MOBILAB. See MOBILAB''s headModel class.'), ...
    arg({'sourceAtlasLabels','SourceAtlasLabels','SourceAtlasROI'},ROINames,ROINames,'Source regions of interest (atlas labels). This is a cell array of strings corresponding to a subset of the labels stored in hmObj.atlas.label. Current source density will be estimated only for these ROIs.','type','logical'), ...
    arg_nogui({'LFM','LeadFieldMatrix'},[],[],'Lead field matrix. Optional. Overrides LFM from HeadModelObject. Otherwise we use hmObj.leadFieldMatrix'), ...
    arg({'channels','Channels'},defChlbl,[],'Channels to keep','type','cellstr','shape','row'), ...
    arg_subswitch({'sourceShape','SourceShape'},{'gausspatch'}, ...
    { ...
    'gausspatch', ...
        { ...
        arg({'sourceCoords','SourceCoordinates'},[],[],'Source coordinates. This is a [num sources x 3] matrix of coordinates of sources. This will overwrite roiAtlasLabels','shape','matrix'), ...
        arg({'roiAtlasLabels','RoiAtlasLabels'},{},ROINames,'ROI Labels','type','logical'), ...
        arg({'roiOrdered','RoiOrdered'},[],[],'ROI labels in specific order. Overrides roiAtlasLabels','type','cellstr','shape','row'), ...
        arg({'nearestNeighbor','NearestNeighbor'},true,[],'Project centroid onto ROI surface'), ...
        arg({'sigma','Sigma'},10,[eps Inf],'Patch standard deviation') ...
        }, ...
    'dipole', ...
        { ...
        arg({'sourceCoords','SourceCoordinates'},[],[],'Source coordinates. This is a [num sources x 3] matrix of coordinates of sources. This will overwrite roiAtlasLabels','shape','matrix'), ...
        arg({'roiAtlasLabels','RoiAtlasLabels'},{},ROINames,'ROI Labels','type','logical'), ...
        arg({'roiOrdered','RoiOrdered'},[],[],'ROI labels in specific order. Overrides roiAtlasLabels','type','cellstr','shape','row'), ...
        arg({'nearestNeighbor','NearestNeighbor'},true,[],'Project centroid onto ROI surface'), ...
        }, ...
    'atlasroi', ...
        { ...
        arg({'roiAtlasLabels','RoiAtlasLabels'},{},ROINames,'ROI Labels','type','logical') ...
        }}, 'Source shape'), ...
    arg_subtoggle({'addNoise','AddNoise'},{},...
        {...
        arg({'snr','SignalToNoise','SNR'},5,[0 Inf],'Signal to noise ratio'), ...
        },'Add noise to scalp data'), ...
    arg({'verb','VerbosityLevel'},int32(2),{int32(0) int32(1) int32(2)},'Verbose','type','int32') ...
    );
%     arg({'cortexSurface','CortexSurface'},[],[],'Cortical surface mesh'), ...
%     arg_nogui({'elocs','ElectrodeLocations'},[],[],'Electrode locations. Optional. [num_channels x 3] (X,Y,Z). Overrides locations from HeadModelObject.'), ...

if islogical(g.sourceAtlasLabels)
    g.sourceAtlasLabels = ROINames(g.sourceAtlasLabels);
end
if isempty(g.sourceShape.roiAtlasLabels)
    g.sourceShape.roiAtlasLabels = ROINames;
end
if isempty(g.sourceAtlasLabels)
    g.sourceAtlasLabels = ROINames;
end

if ~isempty(g.sourceShape.roiOrdered) && ~isempty(g.sourceShape.roiOrdered{1})
    if ~all(ismember_bc(g.sourceShape.roiOrdered,ROINames))
        error('The following ROIs from roiOrdered do not appear in the Atlas: %s',hlp_tostring(g.sourceShape.roiOrdered(~ismember_bc(g.sourceShape.roiOrdered,ROINames))));
    end
    g.sourceShape.roiAtlasLabels = g.sourceShape.roiOrdered;
end

% validate the head model
g.hmObj = hlp_validateHeadModelObject(g.hmObj);

% initialize outputs
centroids_LFM_idx = nan(1,length(g.sourceShape.roiAtlasLabels));
roiVertices       = {};
chaninds          = 1:length(g.hmObj.getChannelLabels());

% Initiate progress bar
waitbarTitle = 'Generating Forward Projection...';
if g.verb
    fprintf(waitbarTitle); end
if g.verb==2
    multiWaitbar(waitbarTitle, ...
                 'Color', [1.0 0.4 0.0], ...
                 'CanCancel','off');
end

% set up model
% -------------------------------------------------------------------------

% load the surfaces
load(g.hmObj.surfaces)
cortexSurface = surfData(3);
nvert = size(cortexSurface.vertices,1);
vertinds = 1:nvert;

% load the lead field matrix
if isempty(g.LFM)
    tmp=load(g.hmObj.leadFieldFile);
    g.LFM=tmp.K; clear tmp;
    
    if ~isempty(g.channels)
        % prune the lead field matrix to contain only desired channels
        hmChanlabels  = lower(g.hmObj.getChannelLabels());
        % use only selected channels that are in the head model
        chaninds = ismember_bc(hmChanlabels,lower(g.channels));
        if nnz(chaninds)~=length(g.channels)
            error('Some channels could not be matched to the headmodel');
        end
        g.LFM = g.LFM(chaninds,:);
    end
end

if g.verb
    fprintf('Generating source coordinates...\n'); end
if g.verb==2
    multiWaitbar(waitbarTitle,1/3); end

% construct the source coordinates
% -------------------------------------------------------------------------
if strcmpi(g.sourceShape.arg_selection,'dipole')
    % handle the dipole case
    % spatial distribution is a gaussian patch with ~0 variance
    % (point-process)
    g.sourceShape.sigma = eps;
    g.sourceShape.arg_selection = 'gausspatch';
end

switch lower(g.sourceShape.arg_selection)
    case 'gausspatch'
        if isempty(g.sourceShape.sourceCoords) && ~isempty(g.sourceShape.roiAtlasLabels)
            numSources = length(g.sourceShape.roiAtlasLabels);
            % postion sources in centroids of each ROI
            for k=1:numSources
                % get the vertices for this ROI
                roiVertices{k} = indices4Structure(hmObj,g.sourceShape.roiAtlasLabels{k})';
                [v,f]          = geometricTools.getSurfaceROI(cortexSurface.vertices,...
                                                           cortexSurface.faces,   ...
                                                           roiVertices{k});
                posxyz = mean(meshcentroid(v,f));
                if size(v,1)==0
                    fprintf('ROI ''%s'' contains no vertices. Skipping...\n', ...
                            g.sourceShape.roiAtlasLabels{k});
                elseif ~isfinite(posxyz)
                    fprintf('Could not compute the centroid of ROI ''%s''. Selecting median vertex instead.\n', ...
                            g.sourceShape.roiAtlasLabels{k});
                    posxyz = v(fix(length(v)/2),:);
                    [~,idx] = ismember_bc(posxyz,cortexSurface.vertices,'rows');
                    centroids_LFM_idx(k) = idx(1);
                elseif g.sourceShape.nearestNeighbor
                    % map the centroid of the ROI onto the ROI surface
                    dt     = DelaunayTri(v(:,1),v(:,2),v(:,3));
                    loc    = nearestNeighbor(dt, posxyz);
                    posxyz = v(loc,:);
                    [~,idx] = ismember_bc(posxyz,cortexSurface.vertices,'rows');
                    centroids_LFM_idx(k) = idx(1);
                end
                g.sourceShape.sourceCoords(k,:) = posxyz;
            end
        end
    case 'atlasroi'
        error('atlasroi option not yet implemented');
end

if g.verb
    fprintf('Generating current density...\n'); end
if g.verb==2
    multiWaitbar(waitbarTitle,2/3); end

% contruct the source patches and generate current density
% -------------------------------------------------------------------------
if size(g.sourceShape.sourceCoords,2)~=3
    error('The number of columns of sourceCoords must be 3 (X, Y, Z)');
end

numSources = size(g.sourceShape.sourceCoords,1);

if isempty(g.sourceAmps)
    g.sourceAmps = ones(numSources,1);
elseif isscalar(g.sourceAmps)
    g.sourceAmps = repmat(g.sourceAmps,[numSources,1]);
elseif isvector(g.sourceAmps)
    if size(g.sourceAmps,2)==numSources && size(g.sourceAmps,1)~=numSources
        % make column vector
        g.sourceAmps = g.sourceAmps(:); 
    end
end
    
if size(g.sourceAmps,1)~=numSources
    error('number of rows of sourceAmps must equal the number of sources');
end

if isscalar(g.sourceShape.sigma)
    g.sourceShape.sigma = repmat(g.sourceShape.sigma,[1,numSources]);
end
if length(g.sourceShape.sigma)~=numSources
    error('number of sigma entries must equal number of sources');
end

nchs  = size(g.LFM,1);
npnts = size(g.sourceAmps,2);
ntr   = size(g.sourceAmps,3);
scalpData = zeros(nchs,npnts*ntr);
srcData   = zeros(nvert,npnts*ntr);

if ntr > 1
    % reshape data to 2D (if 3D)
    g.sourceAmps = g.sourceAmps(:,:);
end

for k=1:numSources
    if any(isnan(g.sourceShape.sourceCoords(k,:)))
        % skip this source
        continue;
    end
    sigma = g.sourceShape.sigma(k);
    P       = zeros(nvert,1);           % source spatial distribution
    J       = zeros(nvert,npnts*ntr);       % current density estimate
    x0 = g.sourceShape.sourceCoords(k,:);
    d = sqrt(sum((cortexSurface.vertices(vertinds,:) - ones(length(vertinds),1)*x0).^2,2));
    P(vertinds) = normpdf(d,0,sigma)*sigma*sqrt(2*pi);  % unit amplitude
    if ~isempty(g.sourceAmps)
        % distribute sourceAmps time series according to source distribution
        J(vertinds,:) = P(vertinds)*g.sourceAmps(k,:);
    else
        J = P;
    end
    % sum projected current over sources
    scalpData = scalpData + g.LFM*J(vertinds,:);
    srcData   = srcData   + J;
end

% adding noise
if g.addNoise.arg_selection
    scalpData = scalpData + (std(scalpData(:))./g.addNoise.snr).*randn(nchs,npnts*ntr);
end

% reshape data back to 3D
if ntr>1
    scalpData = reshape(scalpData,[nchs,npnts,ntr]);
    srcData   = reshape(srcData,[nvert,npnts,ntr]);
    J         = reshape(J,[nvert,npnts,ntr]);
end

% prepare outputs
if nargout > 2
    LFM = g.LFM;
end
if nargout > 6
    sourceCoords = g.sourceShape.sourceCoords;
end

if g.verb==2
    multiWaitbar(waitbarTitle,'Close'); 
end
