% This scriptlet processes a collection of SIFT BrainMovie frames 
% (saved in .fig format). 
% The brain, network, and labels are extracted from the figure, background
% and text color are converted, some optional camlighting modified and 
% figure is saved using export_fig().
%
% Note: this requires the export_fig() package from the Matlab FEX:
% http://www.mathworks.com/matlabcentral/fileexchange/23629
%
% Author: Tim Mullen, 2012, SCCN/INC/UCSD.

% User Options:
dataExportFormat    = 'png';    % any image formats supported by export_fig(). Can also be 'fig' for Matlab figure export
bgColor             = 'w';      % background figure color
textColor           = 'k';      % text color
exportFigArgs       = {};       % additional args for export_fig(). See 'doc export_fig'
suffix = ''; 'image';

% get the path to a collection of brainmovie .fig files
% each frame in the set should be labeled in the format
% '<some-prefix>image<frame-number>.fig'
% e.g. 'MyBrainMovie_image001.fig'
% This is default export from SIFT's vis_causalBrainMovie3D
importFigPath = uigetdir(pwd,'Select Import Folder');
% get path to export folder (where images will be saved)
exportFigPath = uigetdir(pwd,'Select Export Folder');

% get all filenames
dd=dir(fullfile(importFigPath, sprintf('*%s*.fig',suffix)));
fnames = {dd.name};

if isempty(fnames)
    fnames = {''};
end

% load each figure and convert background to white and text to black (.fig)
for k=1:length(fnames)
    
    if importFigPath
        % load figure
        uiopen([importFigPath filesep fnames{k}],1);
    end
    oriFig = gcf;
    
    
    % extract the brain
    newFig=isolate_axes(findall(oriFig,'tag','brain1'),true);
    
    % remove title and make bg white and text black, bold
    ax = findall(newFig,'tag','brain1');
    h  = get(ax,'Title');
    set(h,'String','');
    set(newFig,'color',bgColor);
    set(findall(newFig,'type','text'),'color',textColor);
    set(findall(newFig,'type','text'),'fontweight','bold');
    
    % (optionally) you can replace the lighting with a headlight 
    % (has a nice effect for white background)
    delete(findall(newFig,'type','light'));
    axes(ax);
    camlight headlight;
    
    drawnow;
        
    % save figure
    if exportFigPath
        switch lower(dataExportFormat)
            case 'fig'
                saveas(newFig,fullfile(exportFigPath,strtok(fnames{k},'.')));
            otherwise
                export_fig(fullfile(exportFigPath,strtok(fnames{k},'.')), ...
                    ['-' dataExportFormat], ...
                    newFig,exportFigArgs{:});
        end
        close([newFig oriFig]);
    end
end

