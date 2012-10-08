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
% Author: Tim Mullen and Arnaud Delorme, SCCN, INC, UCSD

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

    fid = fopen([fileparts(which('StartSIFT.m')) filesep 'resources' filesep 'version.txt']);
    version = fscanf(fid,'%s');
    fclose(fid);

    vers = num2str(version);
    
    if nargin < 3
        error('eegplugin_sift requires 3 arguments');
    end;
    
    % run SIFT startup routines
    ok=StartSIFT(false);
    
    if ~ok
        fprintf('SIFT initialization failed!\n');
        pause(2);
        return;
    end
        
    % find import data menu
    % ---------------------
    highlevelmenu = findobj(fig, 'tag', 'tools');
    
    % menu callbacks
    % --------------
    cmd = 'EEG = pop_pre_prepData(EEG);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    PreProc_callback        = [finalcmd catchstrs.new_and_hist];
    
    cmd = 'EEG = pop_est_fitMVAR(EEG,0);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    FitModel_callback       = [finalcmd catchstrs.store_and_hist];
    
    cmd  = 'EEG = pop_est_selModelOrder(EEG,0);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    SelectModelOrder_callback   = [finalcmd catchstrs.store_and_hist];
    
    cmd  = 'pop_est_validateMVAR(EEG,0);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    ValidateModel_callback   = [finalcmd catchstrs.store_and_hist];
    
    cmd = 'EEG = pop_est_mvarConnectivity(EEG);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    Connectivity_callback   = [finalcmd catchstrs.store_and_hist];
    
    cmd = 'EEG = pop_vis_TimeFreqGrid(EEG)';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    TFGrid_callback         = [finalcmd catchstrs.store_and_hist];
    
    
    CausalProjection_callback = 'warndlg2(''Coming Soon!'')';
    
    cmd = 'pop_vis_causalBrainMovie3D(EEG);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    BranMovie_callback      = [finalcmd catchstrs.store_and_hist];
    
    cmd = 'EEG = pop_stat_surrogate(EEG,0);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    SurrogateDistrib_callback = [finalcmd catchstrs.store_and_hist];
    
    cmd = 'EEG = pop_stat_surrogateStats(EEG,0);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    SurrogateStats_callback = [finalcmd catchstrs.store_and_hist];
    
    cmd = 'EEG = pop_stat_analyticStats(EEG,0);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    AnalyticStat_callback   = [finalcmd catchstrs.store_and_hist];
    
    SimpleStat_callback     = 'warndlg2(''Coming Soon!'')';
    
    AboutSIFT_callback = 'gui_splashscreen;';

    % create menus
    % ------------
    menu = uimenu( highlevelmenu, 'label', 'SIFT', 'separator', 'on' ,'userdata', 'startup:off;study:on');
    uimenu( menu, 'label', 'Pre-processing', 'callback', PreProc_callback,'userdata', 'startup:off;study:on');
    modelmenu   = uimenu( menu, 'label', 'Model fitting and validation','userdata', 'startup:off;study:on');
    connectmenu = uimenu( menu, 'label', 'Connectivity'  ,'callback',Connectivity_callback,'userdata', 'startup:off;study:on');
    statmenu    = uimenu( menu, 'label', 'Statistics','userdata', 'startup:off;study:on');
    vismenu     = uimenu( menu, 'label', 'Visualization' ,'userdata', 'startup:off;study:on');
    aboutmenu   = uimenu( menu, 'label', 'About SIFT' ,'callback',AboutSIFT_callback,'userdata', 'startup:off;study:on');
    
    uimenu( vismenu , 'label', 'Time-Frequency Grid', 'callback', TFGrid_callback ,'userdata', 'startup:off;study:on');
    uimenu( vismenu , 'label', 'BrainMovie3D', 'callback', BranMovie_callback ,'userdata', 'startup:off;study:on');
    uimenu( vismenu , 'label', 'Causal Projection', 'callback', CausalProjection_callback, 'enable','off' );
    
    uimenu( statmenu, 'label', 'Surrogate Distributions', 'callback', SurrogateDistrib_callback ,'enable','on','userdata', 'startup:off;study:on');
    uimenu( statmenu, 'label', 'Surrogate Statistics', 'callback', SurrogateStats_callback ,'enable','on','userdata', 'startup:off;study:on');
    uimenu( statmenu, 'label', 'Analytic Statistics', 'callback', AnalyticStat_callback,'separator', 'on', 'enable','on','userdata', 'startup:off;study:on' );
    uimenu( statmenu, 'label', 'Simple Statistics', 'callback', SimpleStat_callback, 'separator', 'on' ,'enable','off');
    
    uimenu( modelmenu, 'label', 'Model Order Selection', 'callback', SelectModelOrder_callback ,'userdata', 'startup:off;study:on');
    uimenu( modelmenu, 'label', 'Fit AMVAR Model', 'callback',FitModel_callback,'userdata', 'startup:off;study:on');
    uimenu( modelmenu, 'label', 'Validate model', 'callback', ValidateModel_callback ,'userdata', 'startup:off;study:on');

 

 
    