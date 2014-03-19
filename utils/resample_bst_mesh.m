function SurfaceDS = resample_bst_mesh(CortexSurfacePath,NewNumVertices,OutputPath,PromptAtlas,RepairMesh)
% Resample a BrainStorm3 cortical surface object using iso2mesh
% This preserves (resamples) the atlases stored in the surface object
%
% OutputFile will be of the form [OutputPath/tess_cortex_xxV.mat] where xx
% is the number of new vertices in the mesh
%
% Author: Tim Mullen, SCCN/INC/UCSD, 10-6-2013

if ~exist('RepairMesh','var')
    RepairMesh = true;
end

CortexSurface = load(CortexSurfacePath);

if exist('PromptAtlas','var') && PromptAtlas
    % get the desired atlas
    atlasNames = {CortexSurface.Atlas.Name};
    qst = '';
    for k=1:length(atlasNames)
        qst = [qst sprintf('[%d] %s \n',k,atlasNames{k})]; end
    qst = [qst 'Enter Atlas Number: '];
    res = input(sprintf('%s',qst));
    % atlasname = atlasNames{res};
    atlasList = res;
else
    atlasList = 1:length(CortexSurface.Atlas);
end

%% resample the mesh and atlas
dsfactor = NewNumVertices/size(CortexSurface.Vertices,1);

SurfaceDS = [];
Scouts    = [];

% resample the surface mesh (requires iso2mesh)
fprintf('Resampling surface mesh (ratio=%0.2f)...',dsfactor);
[SurfaceDS.Vertices,SurfaceDS.Faces] = hlp_microcache('resample_bst_mesh',@meshresample,CortexSurface.Vertices,CortexSurface.Faces,dsfactor);
fprintf('done\n');

if RepairMesh
    fprintf('Checking and repairing surface mesh...');
    [SurfaceDS.Vertices,SurfaceDS.Faces]=meshcheckrepair(SurfaceDS.Vertices,SurfaceDS.Faces);
    fprintf('done\n');
end

% resample each atlas
for ai = 1:length(atlasList)
    bstAtlas = CortexSurface.Atlas(atlasList(ai));
    fprintf('Resampling atlas %s...',bstAtlas.Name);
    if isempty(bstAtlas.Scouts)
        fprintf('empty atlas\n');
        continue;
    end
    % assign each vertex to its unique atlas index (roi)
    atlasidx = zeros(size(CortexSurface.Vertices,1),1);
    for k=1:length(bstAtlas.Scouts)
        atlasidx(bstAtlas.Scouts(k).Vertices) = k; 
    end
    % find a nearest-neighbor mapping from the atlas labels for the
    % original mesh vertices to the atlas labels for the resampled mesh
    % vertices
    F = TriScatteredInterp(CortexSurface.Vertices(:,1), ...
                           CortexSurface.Vertices(:,2), ...
                           CortexSurface.Vertices(:,3), ...
                           atlasidx,'nearest');
    atlasidx = F(SurfaceDS.Vertices(:,1),   ...
                 SurfaceDS.Vertices(:,2),   ...
                 SurfaceDS.Vertices(:,3));
    % convert atlas back to brainstorm format
    for k=1:length(bstAtlas.Scouts)
        bstAtlas.Scouts(k).Vertices = hlp_vec(find(atlasidx==k));
        if isempty(bstAtlas.Scouts(k).Vertices)
            bstAtlas.Scouts(k).Seed     = [];
            continue;
        end
%       % create new seeds (set them to the centroid of the mesh)
        [v,f] = geometricTools.getSurfaceROI(...
                    SurfaceDS.Vertices,...
                    SurfaceDS.Faces,   ...
                    bstAtlas.Scouts(k).Vertices);
        if isempty(f)
            % use first vertex in scout as Seed
            bstAtlas.Scouts(k).Seed = bstAtlas.Scouts(k).Vertices(1);
        else
            % use scout centroid (projected onto cortical surface) as Seed
            posxyz = mean(meshcentroid(v,f));
            dt     = DelaunayTri(v(:,1),v(:,2),v(:,3));
            if isempty(dt.Triangulation)
                % use first vertex in scout as Seed
                bstAtlas.Scouts(k).Seed = bstAtlas.Scouts(k).Vertices(1);
            else
                loc    = nearestNeighbor(dt, posxyz);
                posxyz = v(loc,:);
                [~,idx] = ismember_bc(posxyz,SurfaceDS.Vertices,'rows');
                bstAtlas.Scouts(k).Seed = idx(1);
            end
        end
    end
    SurfaceDS.Atlas(ai) = bstAtlas;
    fprintf('done\n');
end

SurfaceDS.iAtlas  = length(atlasList);
SurfaceDS.Comment = sprintf('cortex_%dV',size(SurfaceDS.Vertices,1));
%SurfaceDS.Comment = regexprep(CortexSurface.Comment,'_[0-9]+V',sprintf('_%dV',size(SurfaceDS.Vertices,1)));

if exist('OutputPath','var') && ~isempty(OutputPath)
    outfile = fullfile(OutputPath,['tess_' SurfaceDS.Comment]);
    fprintf('Writing surface file to: %s\n',outfile);
    save(outfile,'-struct','SurfaceDS');
end
