function cancel = hlp_confirmWaitbarCancel(waitbarTitle)
% Confirm cancellation of a waitbar with specified title
% If cancellation is confirmed, waitbar will be closed. otherwise, cancel
% status is reset
% Author: Tim Mullen, 2012, SCCN/INC, UCSD

if strcmpi('yes',questdlg2( ...
        'Are you sure you want to cancel?', ...
        'Confirmation','Yes','No','No'));
    multiWaitbar(waitbarTitle,'Close');
    cancel = true;
else
    multiWaitbar(waitbarTitle,'ResetCancel',true);
    cancel = false;
end

