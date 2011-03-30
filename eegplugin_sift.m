% eegplugin_sift() - EEGLAB plugin for Source Information Flow Toolbox
%
% Usage:
%   >> eegplugin_sift(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Notes:
%   This plugins consist of the following Matlab files:
%
% Create a plugin:
%   For more information on how to create an EEGLAB plugin see the
%   help message of eegplugin_besa() or visit http://www.sccn.ucsd.edu/eeglab/contrib.html
%
% Author: Tim Mullen, SCCN, INC, UCSD

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2010 Tim Mullen
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function vers = eegplugin_sift(fig, trystrs, catchstrs)

    vers = 'sift';
    if nargin < 3
        error('eegplugin_sift requires 3 arguments');
    end;
    
    % add folder to path
    % ------------------
    if ~exist('vis_TimeFreqGrid')
        p = which('eegplugin_sift.m');
        p = p(1:findstr(p,'eegplugin_sift.m')-1);
        addpath(genpath(p));
    end;
    
    % find import data menu
    % ---------------------
    highlevelmenu = findobj(fig, 'tag', 'tools');
    
    % menu callbacks
    % --------------
    cmd = 'EEG = pop_pre_prepData(EEG);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    PreProc_callback        = [finalcmd catchstrs.new_and_hist];
    
    cmd = 'EEG = pop_est_mvarConnectivity(EEG);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    Connectivity_callback   = [finalcmd catchstrs.store_and_hist];
    
    
%     cmd = 'pop_vis_TimeFreqGrid(EEG,cfg);';
    cmd = ['try, cfg = EEG(1).CAT.configs.TimeFreqGrid; catch, cfg = []; end;' ...
               'EEG(1).CAT.tmpcfg = pop_vis_TimeFreqGrid(EEG,cfg);' ...
               'if ~isempty(EEG(1).CAT.tmpcfg), ' ...
               ' for cnd = 1:length(EEG), EEG(cnd).CAT.configs.TimeFreqGrid = EEG(1).CAT.tmpcfg; end; end; ' ...
               ' EEG(1).CAT = rmfield(EEG(1).CAT,''tmpcfg''); '];
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' strrep(cmd,'''','''''') ''';' ];
    TFGrid_callback         = [finalcmd catchstrs.store_and_hist];
    
    
    CausalProjection_callback = 'warndlg2(''Coming Soon!'')';
    
%     cmd = 'pop_vis_causalBrainMovie3D(EEG,cfg);';
    cmd = ['try, cfg = EEG(1).CAT.configs.BrainMovie3D; catch, cfg = []; end;' ...
               'EEG(1).CAT.tmpcfg = pop_vis_causalBrainMovie3D(EEG,cfg);' ...
               'if ~isempty(EEG(1).CAT.tmpcfg), ' ...
               ' for cnd = 1:length(EEG), EEG(cnd).CAT.configs.BrainMovie3D = EEG(1).CAT.tmpcfg; end; end; ' ...
               ' EEG(1).CAT = rmfield(EEG(1).CAT,''tmpcfg''); '];
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' strrep(cmd,'''','''''') ''';' ];
    BranMovie_callback      = [finalcmd catchstrs.store_and_hist];
    
    BootStat_callback       = 'warndlg2(''Coming Soon!'')';
    AnalyticStat_callback   = 'warndlg2(''Coming Soon!'')';
    SimpleStat_callback     = 'warndlg2(''Coming Soon!'')';
    
    cmd = 'EEG = pop_est_fitMVAR(EEG,0);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    FitModel_callback       = [finalcmd catchstrs.store_and_hist];
    
    ValidateModel_callback  = 'pop_est_validateMVAR(EEG,0);';
    
    
    % create menus
    % ------------
    menu = uimenu( highlevelmenu, 'label', 'SIFT', 'separator', 'on' );
    uimenu( menu, 'label', 'Pre-processing', 'callback', PreProc_callback);
    modelmenu   = uimenu( menu, 'label', 'Model fitting and validation');
    connectmenu = uimenu( menu, 'label', 'Connectivity'  ,'callback',Connectivity_callback);
    statmenu    = uimenu( menu, 'label', 'Statistics');
    vismenu     = uimenu( menu, 'label', 'Visualization' );
    uimenu( vismenu , 'label', 'Time-Frequency Grid', 'callback', TFGrid_callback );
    uimenu( vismenu , 'label', 'BrainMovie3D', 'callback', BranMovie_callback );
    uimenu( vismenu , 'label', 'Causal Projection', 'callback', CausalProjection_callback, 'enable','off' );
    
    uimenu( statmenu, 'label', 'Bootstrap', 'callback', BootStat_callback ,'enable','off');
    uimenu( statmenu, 'label', 'Analytic', 'callback', AnalyticStat_callback, 'enable','off' );
    uimenu( statmenu, 'label', 'Simple statistics', 'callback', SimpleStat_callback, 'separator', 'on' ,'enable','off');
    
    uimenu( modelmenu, 'label', 'Fit AMVAR Model', 'callback',FitModel_callback);
    uimenu( modelmenu, 'label', 'Validate model', 'callback', ValidateModel_callback );

 

 
    