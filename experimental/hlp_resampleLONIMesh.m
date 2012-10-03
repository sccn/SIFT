function SurfaceDS = hlp_resampleLONIMesh(Surface,decimationRate)
% Author: Alejandro Ojeda and Tim Mullen, SCCN/ICN, UCSD 2012
% requires iso2mesh Matlab toolbox
% Surface is a LONI probablistic atlas mesh containing fields
%   .faces
%   .vertices
%   .labels
%   .colortable

if nargin<2
    decimationRate = 0.01; % reducing the mesh by 1%
end

% downsample the surface mesh
[SurfaceDS.vertices,SurfaceDS.faces]=meshresample(Surface.vertices,Surface.faces,decimationRate);

% update the color labeling
F = TriScatteredInterp(Surface.vertices(:,1),Surface.vertices(:,2),Surface.vertices(:,3),Surface.label,'nearest');
SurfaceDS.label = F(SurfaceDS.vertices(:,1),SurfaceDS.vertices(:,2),SurfaceDS.vertices(:,3));
SurfaceDS.colortable = Surface.colortable;

