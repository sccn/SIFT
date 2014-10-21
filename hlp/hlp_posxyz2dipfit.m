function dipfit = hlp_posxyz2dipfit(posxyz,coordformat)
% store [x y z] source locations in a "minimal" EEGLAB dipfit structure
% 
% Inputs: 
%   posxyz:      [N x 3] matrix of [X Y Z] coordinates for N sources
%   coordformat: (optional) coordinate format string {def: MNI}
% Output:
%   dipfit:      EEGLAB dipfit structure
%
% Tim Mullen, SCCN/INC 2014

if nargin<2
    coordformat = 'MNI';
end

% construct the dipfit model
dipfit.hdmfile  = '';
dipfit.mrifile  = '';
dipfit.chanfile = '';
dipfit.chansel  = [];
dipfit.coordformat = coordformat;
dipfit.coord_transform = [0 0 0 0 0 0 1 1 1];

numDipoles = size(posxyz,1);
for k=1:numDipoles
    dipfit.model(k).posxyz      = posxyz(k,:);
    dipfit.model(k).momxyz      = eps;
    dipfit.model(k).rv          = 0;
    dipfit.model(k).select      = 1;
    dipfit.model(k).diffmap     = [];
    dipfit.model(k).sourcepot   = [];
    dipfit.model(k).datapot     = [];
end