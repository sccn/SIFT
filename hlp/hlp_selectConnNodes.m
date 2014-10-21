function CAT = hlp_selectConnNodes(CAT,nodeIndsToKeep,connmethods,fields,keepNodeOrder)
% remove some nodes from N-D tensors stored in a SIFT structure
% can also be used to permute the order of nodes in a datastructure
%
% Inputs:
%
%     CAT:            SIFT data object or EEG data structure with 'CAT' field
%     nodesIdsToKeep: (optional) indices of nodes to keep {default: keep all}
%     connmethods:    (optional) a list of measures to prune {default: prune all}
%     fields:         (optional) a specific field of CAT to prune (Conn, PConn, Pnull, Stats, ...) {default: prune all of the above)
%     keepNodeOrder:  (optional) preserve the node ordering. If false, order will be returned in order of nodesIndsToKeep {default: true}
%
% Outputs:
% 
%     CAT:            Pruned dataset
%
% Author: Tim Mullen, SCCN/INC/UCSD 2014

if nargin<2 || isempty(nodeIndsToKeep)
    return;  end
if nargin<3
    connmethods = []; end
if nargin<4 || isempty(fields)
    fields = {}; end
if nargin<5
    keepNodeOrder = true;
end

if isfield(CAT,'CAT') && isfield(CAT,'data')
    % this is an EEGLAB datastructure
    CAT.CAT = hlp_selectConnNodes(CAT.CAT,nodeIndsToKeep,connmethods,fields);
    return;
end

if isempty(fields)
    fields = {};
    % ensure we have a Conn struct
    res = hlp_checkeegset(struct('CAT',CAT),'conn');
    if isempty(res), fields = [fields 'Conn']; end

    % check if we have a PConn struct
    res = hlp_checkeegset(struct('CAT',CAT),'pconn');
    if isempty(res), fields = [fields 'PConn']; end

    % check if we have a Pnull struct
    res = hlp_checkeegset(struct('CAT',CAT),'pnull');
    if isempty(res), fields = [fields 'Pnull']; end

    % check if we have a Stats struct
    res = hlp_checkeegset(struct('CAT',CAT),'stats');
    if isempty(res), fields = [fields 'Stats']; end

    if isempty(fields)
        error('Please ensure you have Conn, PConn, Pnull, or Stats structures'); 
    end
else
    % sanity check
    for k=1:length(fields)
        if ~isfield(CAT,fields{k})
            error('SIFT:hlp_selectConnNodes','%s is not a valid field of CAT object',fields{k}); 
        end
    end
end

if keepNodeOrder
    % ensure we preserve the original ordering
    nodeIndsToKeep = sort(nodeIndsToKeep,'ascend');
end

for fn = fields
    fnk = fn{1};
    % determine measures to prune
    if ~isempty(connmethods)
        if ~iscell(connmethods)
            connmethods = {connmethods};
        end
    else
        connmethods = hlp_getConnMethodNames(CAT.(fnk));
    end

    try
        for m=1:length(connmethods)
            cm = connmethods{m};
            % prune this conn method
            if strcmpi(fnk,'stats')
                % we have to process stats a bit differently
                % since it has it's own sub-fields
                subfn = fieldnames(CAT.(fnk).(cm)); % get sub-fields
                for smi = 1:length(subfn)
                    sm = subfn{smi};
                    if strcmpi(sm,'ci')
                        % one caveat is the Stats.(method).ci array which stores [upper lower]
                        % bounds in first dimension
                        CAT.(fnk).(cm).(sm) = CAT.(fnk).(cm).(sm)(:,nodeIndsToKeep,nodeIndsToKeep,:,:,:,:,:);
                    else
                        CAT.(fnk).(cm).(sm) = CAT.(fnk).(cm).(sm)(nodeIndsToKeep,nodeIndsToKeep,:,:,:,:,:);
                    end
                end
            else
                CAT.(fnk).(cm) = CAT.(fnk).(cm)(nodeIndsToKeep,nodeIndsToKeep,:,:,:,:,:);
            end
        end
    catch err
        if strcmp(err.identifier,'MATLAB:badsubscript')
            error('SIFT:hlp_selectConnNodes', ...
                  'The list of nodes you provided does not agree with the number of nodes available for CAT.%s.%s',fnk,cm);
        else
            error('SIFT:hlp_selectConnNodes',err.message);
        end
    end
end

% clean up the rest of the CAT structure
CAT.curComponentNames = CAT.curComponentNames(nodeIndsToKeep);
CAT.curComps = CAT.curComps(nodeIndsToKeep);
CAT.nbchan = length(nodeIndsToKeep);
try CAT.srcdata = CAT.srcdata(nodeIndsToKeep,:,:,:); catch; end


