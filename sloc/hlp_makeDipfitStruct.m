function dipfit = hlp_makeDipfitStruct(sourceSpace,roiVertices,CSD)
% create dipfit structure containing the locations (.posxyz) and moments 
% (.momxyz) of a set of dipoles located at each vertex of the source space.  
% Alternately, one can provide a cell array of integers (roiVertices) where
% each cell contains the indices of vertices within a given region of 
% interest. If this is provided, then dipfit.model(k).posxyz contains the
% centroid location of the kth ROI and momxyz contains its mean moment.
    
% Inputs:
%   sourceSpace:    Surface structure containing fields 
%                     .vertices [num_vertices x 3]
%                     .faces [num_faces x 3]
% Optional:
%   roiVertices:    A cell array of integers (roiVertices) where each cell 
%                   contains the indices of vertices within a given ROI
%   CSD:            [num_vertices x num_samples] current source density
%
% Outputs:
%   dipfit:         dipfit model (struct array) containing fields
%                   .model.posxyz     dipole location
%                   .model.momxyz     dipole moment
%
% Author: Tim Mullen, 2013, SCCN/INC/UCSD

if nargin<1
    error('You must provide a sourceSpace');
end
if nargin<2
    roiVertices = {};
end
if nargin<3
    CSD = [];
end

% determine number of dipoles
if isempty(roiVertices)
    numDipoles = size(sourceSpace.vertices,1);
else
    numDipoles = length(roiVertices);
end


% compute normal vectors
if ~isempty(CSD)
    normals = geometricTools.getSurfaceNormals(sourceSpace.vertices,sourceSpace.faces,false);
    if size(CSD,1) == 3*size(sourceSpace.vertices,1)  
        CSD = reshape(CSD,[size(CSD,1)/3 3]);
        momxyz  = sqrt(sum(CSD.^2,2));
        s   = sign(dot(normals,CSD,2));
        momxyz  = s.*momxyz; % momxyz are the moment vectors (scaled by power)
%         obj.hVector = quiver3(sourceSpace.vertices(:,1),sourceSpace.vertices(:,2),sourceSpace.vertices(:,3),CSD(:,1),CSD(:,2),CSD(:,3),2);
    else
        momxyz = CSD;
%         obj.hVector = quiver3(sourceSpace.vertices(:,1),sourceSpace.vertices(:,2),sourceSpace.vertices(:,3),normals(:,1),normals(:,2),normals(:,3),2);
    end
else
    momxyz = zeros(numDipoles,3);
end

% compute ROI averages
if ~isempty(roiVertices)
    
    % average vertex locations (spatial centroid)
    % obtain surfaces for each ROI
    for k=1:length(roiVertices)
        [nVertices,nFaces] = geometricTools.getSurfaceROI(sourceSpace.vertices,sourceSpace.faces,roiVertices{k});
        
%         posxyz = cellfun(@(x) mean(sourceSpace.vertices(x,:),1)', ...
%                          roiVertices, 'UniformOutput',false);
        posxyz = cell2mat(posxyz)';
    
    % average vertex moments
    if ~isempty(CSD)
        momxyz = cellfun(@(x) mean(momxyz(x,:),1)', ...
                     roiVertices, 'UniformOutput',false);
        momxyz = cell2mat(momxyz)';
    end
end
       
% construct the dipole model
dipfit.hdmfile  = '';
dipfit.mrifile  = '';
dipfit.chanfile = '';
dipfit.chansel  = [];
dipfit.coordformat = 'MNI';
dipfit.coord_transform = [0 0 0 0 0 0 1 1 1];

for i=1:numDipoles
    dipfit.model(i).posxyz      = posxyz(i,:);
    dipfit.model(i).momxyz      = momxyz(i,:);
    dipfit.model(i).rv          = 0;
    dipfit.model(i).select      = 1;
    dipfit.model(i).diffmap     = [];
    dipfit.model(i).sourcepot   = [];
    dipfit.model(i).datapot     = [];
end


function centroid=meshcentroid(v,f)
%
% centroid=meshcentroid(v,f)
% 
% compute the centroids of a mesh defined by nodes and elements
% (surface or tetrahedra) in R^n space
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% input:
%      v: surface node list, dimension (nn,3)
%      f: surface face element list, dimension (be,3)
%
% output:
%      centroid: centroid positions, one row for each element
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

ec=reshape(v(f(:,1:size(f,2))',:)', [size(v,2) size(f,2) size(f,1)]);
centroid=squeeze(mean(ec,2))';
