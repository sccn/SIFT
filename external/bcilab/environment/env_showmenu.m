function env_showmenu(varargin)
% Links the BCILAB menu into another menu, or creates a new root menu if necessary.
% env_showmenu(Options...)
%
% In:
%   Options... : optional name-value pairs; names are:
%                 'parent': parent menu to link into
%
%                 'shortcuts': whether to enable keyboard shortcuts
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-10-29

if ~exist('varargin','var')
    varargin = {}; end

% parse options...
hlp_varargin2struct(varargin,'parent',[], 'shortcuts',true);

if ~isempty(findobj('Tag','bcilab_menu'))
    return; end

if isempty(parent)
    % create new root menu
    mp = get(0,'MonitorPositions');
    parent = figure('DockControls','off','NumberTitle','off','Name','BCILAB 0.9','Resize','off','MenuBar','none','Position',[100 mp(1,4)-100 400 1]);
end

% set up the menu

source = uimenu(parent, 'Label','Data Source','Tag','bcilab_menu');
uimenu(source,'Label','Load recording(s)...','Accelerator',char(shortcuts*'l'),'Callback','gui_loadset');
uimenu(source,'Label','Load study...','Callback','gui_loadstudy','Enable','off');
uimenu(source,'Label','Define marker transform...','Separator','on','Accelerator',char(shortcuts*'d'),'Callback','gui_eventinsertion');

offline = uimenu(parent, 'Label','Offline Analysis');
uimenu(offline,'Label','New approach...','Accelerator',char(shortcuts*'n'),'Callback','gui_newapproach');
uimenu(offline,'Label','Modify approach...','Accelerator',char(shortcuts*'m'),'Callback','gui_configapproach');
uimenu(offline,'Label','Review/edit approach...','Accelerator',char(shortcuts*'r'),'Callback','gui_reviewapproach');
uimenu(offline,'Label','Save approach...','Accelerator',char(shortcuts*'s'),'Callback','gui_saveapproach');
uimenu(offline,'Label','Train new model...','Accelerator',char(shortcuts*'t'),'Callback','gui_calibratemodel','Separator','on');


online = uimenu(parent,'Label','Online Analysis');
process = uimenu(online,'Label','Process data within...');
uimenu(process,'Label','DataRiver...','Callback','gui_rundatariver');
uimenu(process,'Label','BCI2000...','Callback','gui_runbci2000','Enable','off');
uimenu(process,'Label','OpenVibe...','Callback','gui_runopenvibe','Enable','off');
receive = uimenu(online,'Label','Receive input from...');
uimenu(receive,'Label','DataRiver...','Callback','gui_receivedatariver','Enable','off');
uimenu(receive,'Label','BrainVision Recorder...','Callback','gui_receivebva','Enable','off');
uimenu(receive,'Label','TCP...','Callback','gui_receivetcp','Enable','off');
produce = uimenu(online,'Label','Provide output to...');
uimenu(produce,'Label','DataRiver...','Callback','gui_producedatariver','Enable','off');
uimenu(produce,'Label','TCP...','Callback','gui_producetcp','Enable','off');
uimenu(produce,'Label','EPrime...','Callback','gui_produceeprime','Enable','off');
uimenu(produce,'Label','Paradigm...','Callback','gui_produceparadigm','Enable','off');
uimenu(produce,'Label','FLXlab...','Callback','gui_produceflx','Enable','off');