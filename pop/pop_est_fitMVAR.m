function varargout = pop_est_fitMVAR(ALLEEG,typeproc,varargin)
%
% Fit an (adaptive) multivariate autoregressive model (MVAR) to the data.
% If only two inputs are provided, an input GUI will popup.
%
% Input:
%
%   EEG                Preprocessed EEG structure.
%   typeproc           Reserved for future use. Use 0
%
% optional:
%
%   ('name', value) pairs to pass to est_fitMVAR. See help est_fitMVAR. If
%   any name,value pair is provided, GUI will be supressed and est_fitMVAR
%   directly called.
%
% Output:
%
%   ALLEEG:         EEG structure with .CAT.MODEL substructure containing
%                   fitted VAR[p] model
%   cfg:            arguments structure containing parameters used to fit
%                   model
%
% See Also: est_fitMVAR()
%
% References: 
% 
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
% 
% Author: Tim Mullen, 2010, SCCN/INC, UCSD. 
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
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


varargout{1} = ALLEEG;

if nargin<2
    typeproc = 0;
end

popup = nargin<3;

if nargin<3
    % use default values
    g = hlp_getDefaultParams('est');
else
    % user must input ALL values
    var = hlp_mergeVarargin(varargin{:});
    g = finputcheck(var, hlp_getDefaultArglist('est'), 'pop_est_fitMVAR','ignore','quiet');
    if ischar(g), error(g); end
    if isempty(g.epochTimeLims), g.epochTimeLims = [0 ALLEEG(1).pnts/ALLEEG(1).srate]; end
    if isempty(g.morder) || length(g.morder)>2, error('invalid entry for field ''morder'''); end
end

% if a new method is added, add here and in hlp_getDefaultArglist.m
mvarAlgorithms = {'vieira-morf'};
% names of algorithms in exact order they appear in mvarAlgorithms
mvarAlgsFullNames = {'Vieira-Morf'};

if exist('arfit')==2
    mvarAlgorithms    = [mvarAlgorithms 'arfit'];
    mvarAlgsFullNames = [mvarAlgsFullNames 'ARFIT'];
end

% check if Jacket is installed
has_jacket = 1;
try
    ginfo;
catch
    has_jacket =0;
end

if ischar(ALLEEG)
    command = ALLEEG;
    fig     = typeproc;
    ALLEEG  = get(fig, 'userdata');
    
    if strcmpi(command, 'redraw')
        %         hld_cmdmodelorder = findobj(fig, 'tag', 'cmdModelOrder');
        %         hld_winlen = findobj(fig, 'tag', 'edtWindowLength');
        %         hld_winstep = findobj(fig, 'tag', 'edtStepSize');
        %
        %         g = get(hld_cmdmodelorder,'userdata');
        %         g.winlen = str2num(get(hld_winlen,'string'));
        %         g.winstep = str2num(get(hld_winstep,'string'));
        %         set(hld_cmdmodelorder, 'userdata',g);
        
    elseif strcmpi(command, 'modorder')
        %         hld_cmdmodelorder = findobj(fig, 'tag', 'cmdModelOrder');
        hld_winlen = findobj(fig, 'tag', 'edtWindowLength');
        hld_winstep = findobj(fig, 'tag', 'edtStepSize');
        hld_algorithm = findobj(fig,'tag','popMVARalgs');
        hld_morder = findobj(fig,'tag','edtModelOrder');
        
        g.winlen = str2num(get(hld_winlen,'string'));
        g.winstep = str2num(get(hld_winstep,'string'));
        g.algorithm = mvarAlgorithms{get(hld_algorithm,'value')};
        g.morder = str2num(get(hld_morder,'string'));
        
        %         if strcmpi(g.algorithm,'10') && ALLEEG(1).trials>1
        %             errordlg2('ARFIT cannot currently be used with multiple-trial data. Please select another algorithm.');
        %             return;
        %         end
        
        
        if strcmpi(g.algorithm,'vieira-morf-pll') && ~has_jacket
            errordlg2('Parallel version only compatible with Jacket installation. Please select another algorithm.');
            return;
        end
        
        
        if isempty(g.winlen) || isempty(g.winstep)
            errordlg2('Window length and step size must be specified');
            return;
        end
        
        % check parameters
        for cond=1:length(ALLEEG)
            fprintf('===================================================\n');
            fprintf('MVAR PARAMETER SUMMARY FOR CONDITION: %s\n',ALLEEG(cond).condition);
            fprintf('===================================================\n');
            checkcode = checkMVARParams(ALLEEG(cond),g);
            fprintf('\n\n')
            
            if isobject(checkcode)
                errordlg2(checkcode.message,'Error in MVAR configuration!');
                return;
            end
        end

        switch checkcode
            case 'error'
                % generate error
                errordlg2('One or more parameters are invalid (see command-window for details)','Checking MVAR parameters...');
                return;
                % go back to main input GUI
            case 'warning'
                % if OK is pressed continue onward, otherwise, go back to main input GUI
                res=questdlg2('Some warnings were generated (see command-window for details), Continue?','Checking MVAR parameters', 'Cancel', 'OK', 'OK');
                if strcmpi(res,'cancel')
                    return
                end
        end
        
        % all-clear, continue to model order selection
        pop_est_selModelOrder(ALLEEG,0,g);
    end;
    
    return;
end


if popup
    
    cb_winlen = 'pop_est_fitMVAR(''redraw'', gcbf);';
    cb_winstep = 'pop_est_fitMVAR(''redraw'',gcbf);';
    
    
    geomhoriz = {1 1 [3 2] 1 [3 2] [3 2] 1};
    uilist = { ...
        { 'Style', 'text',       'string', '1. Select MVAR algorithm'}...
        { 'Style', 'popup',      'string', mvarAlgsFullNames, 'tag', 'popMVARalgs','Value',1} ...
        { 'Style', 'text',       'string', '2. Window length (sec) '}...
        { 'Style', 'edit',       'string', num2str(ALLEEG(1).xmax-ALLEEG(1).xmin) ,'tag', 'edtWindowLength','callback',cb_winlen}...
        { 'Style', 'pushbutton', 'string', 'Start Window Length Assistant...' ,'tag', 'cmdWindowLength','callback','warndlg2(''Coming soon!'')'}...
        { 'Style', 'text',       'string', '3. Step size (sec) '}...
        { 'Style', 'edit',       'string', num2str(0.06*(ALLEEG(1).xmax-ALLEEG(1).xmin)) ,'tag', 'edtStepSize','callback',cb_winstep}...
        { 'Style', 'text',       'string', '4. Model order '}...
        { 'Style', 'edit',       'string', '5' ,'tag', 'edtModelOrder'}...
        { 'Style', 'pushbutton', 'string', 'Start Model Order Assistant...' ,'tag', 'cmdModelOrder','callback','pop_est_fitMVAR(''modorder'', gcbf);'} ...
        };
    
    
    options = { 'geometry', geomhoriz, 'geomvert',[1 1 1 1 1 1 1], 'uilist',uilist, 'helpcom','pophelp(''pop_est_fitMVAR'');', ...
        'title','Fit AMVAR Model', 'mode', 'plot', 'userdata', ALLEEG };
    inputgui( options{:} );
    fig = gcf;
    
    cont = true;
    fig = gcf;
    while cont
        waitfor( findobj('parent', fig, 'tag', 'ok'), 'userdata');
        try findobj(fig); % figure still exist ?
        catch, return; end;
        
        
        [tmp1 tmp2 strhalt result] = inputgui('getresult', fig, options{:} );
        g.winlen = str2double(result.edtWindowLength);
        g.winstep = str2double(result.edtStepSize);
        g.algorithm = mvarAlgorithms{result.popMVARalgs};
        g.morder = str2double(result.edtModelOrder);
        
        set(findobj('parent', fig, 'tag', 'ok'), 'userdata', '');
        
        if strcmpi(g.algorithm,'vieira-morf-pll') && ~has_jacket
            errordlg2('Parallel version only compatible with Jacket installation. Please select another algorithm.');
            return;
        end
        
        if any(isnan(g.morder)) || length(g.morder)>=2 || rem(g.morder,1)
            errordlg2('Please specify a single, positive, integer model order');
            continue;
        end
        
        % check parameters
        for cond=1:length(ALLEEG)
            fprintf('===================================================\n');
            fprintf('MVAR PARAMETER SUMMARY FOR CONDITION: %s\n',ALLEEG(cond).condition);
            fprintf('===================================================\n');
            checkcode = checkMVARParams(ALLEEG(cond),g);
            fprintf('\n\n')
            
            if isobject(checkcode)
                errordlg2(checkcode.message,'Error in MVAR configuration!');
                cont = true;
                continue;
            end
        end

        switch checkcode
            case 'error'
                % generate error
                errordlg2('One or more parameters are invalid (see command-window for details)','Checking MVAR parameters...');
                % go back to main input GUI
            case 'warning'
                % if OK is pressed continue onward, otherwise, go back to main input GUI
                res=questdlg2('Some warnings were generated (see command-window for details), Continue?','Checking MVAR parameters', 'Cancel', 'OK', 'OK');
                if strcmpi(res,'OK')
                    cont=false; % exit loop
                end
            case 'ok'
                % no warnings, exit loop;
                cont=false;
        end

    end;
    
    close(fig);
    
    if isempty( tmp1 ), return; end;
end


% fit the MVAR model
for cond=1:length(ALLEEG)
    fprintf('analyzing condition %s...\n',ALLEEG(cond).condition);
    [ALLEEG(cond).CAT.MODEL] = est_fitMVAR(ALLEEG(cond),typeproc,g);
end


varargout{1} = ALLEEG;
if nargout<2
    for cond=1:length(ALLEEG)
        ALLEEG(cond).CAT.params = g;
    end
else
    varargout{2} = g;
end




% perform sanity checks on MVAR parameters
function checkcode = checkMVARParams(ALLEEG,g)

    checkcode = 'ok';
    
    % check that parameters are OK
    try
        [infostring warnstring errstring] = est_checkMVARParams(ALLEEG,g);
    catch err
%         errordlg2(err.message,'Error in MVAR configuration!');
        checkcode = err;
        return;
    end


    if g.verb>0
        % display summary on command line
        for str=1:length(infostring)
            if ~isempty(errstring{str})
                fprintf('ERROR\t');
            elseif ~isempty(warnstring{str})
                fprintf('WARNING\t');
            else
                fprintf('OK\t');
            end
            fprintf(infostring{str});
            fprintf(warnstring{str});
            fprintf(errstring{str});
            fprintf('\n');
        end
    end
    
    if ~all(ismember(errstring,''))
        checkcode = 'error';
    elseif ~all(ismember(warnstring,''))
        checkcode = 'warning';
    end
    
    
    