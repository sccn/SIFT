function ok = StartSIFT(interactive,nosplash)
% This function initializes the Source Information Flow Toolbox (SIFT)
% Author: Tim Mullen, 2011, SCCN/INC/UCSD

ok = false;

if nargin<1
    interactive = false;
end
if nargin<2
    nosplash = false;
end

fprintf('Initializing SIFT...\n');
    
% temporary hack to remove measure projection toolbox from path 
% (MPT interferes with SIFT)
% if exist('measure_projection')
%     hlp_removeMeasureProjection;
% end

if ispc, SLASH = '\'; else SLASH = '/'; end
siftroot    = fileparts(which('StartSIFT.m'));
ArfitURL    = 'http://www.gps.caltech.edu/~tapio/arfit/arfit.zip';
ARfitTargetPath  = [siftroot SLASH 'external' SLASH 'arfit'];
        
% add subfolders to path (if not already present)
if ~exist('vis_TimeFreqGrid.m','file')
%     p = which('StartSIFT.m');
%     p = p(1:strfind(p,'StartSIFT.m')-1);
    addpath(genpath(siftroot));
end
    
% call up the splash screen
if interactive && ~nosplash
    gui_splashscreen;
end

% optionally download arfit (if not already present)
if ~exist(ARfitTargetPath,'file') && interactive
    res = input('SIFT: Would you like to download and install the ARFIT toolbox as a SIFT plugin (recommended)? ''y''/''n'': ','s');
    if strcmpi(res,'y')
       
        fprintf('SIFT: Downloading and installing ARFIT from %s ...\n',ArfitURL);
        
        try 
            outdir = unzip(ArfitURL,ARfitTargetPath);
        catch e
            switch e.identifier
                case 'MATLAB:unzip:urlwriteError'
                    fprintf('%s Aborting download.\n',e.message);
                case 'MATLAB:unzip:invalidZipFile'
                    fprintf('Unable to unpack file %s Aborting installation.\n', e.message);
                    res=rmdir(ARfitTargetPath,'s'); % clean up
            end
            
            return;
        end
        
        addpath(genpath(ARfitTargetPath));
        
        fprintf('SIFT: ARFIT installed to %s%c\n',fileparts(outdir{1}),SLASH);
    end
end

ok = true;

fprintf('Start SIFTing!\n');

