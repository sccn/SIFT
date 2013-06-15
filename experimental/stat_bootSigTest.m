function Stats = stat_bootSigTest(PConn,multipleComparisonCorrection)

% compute the difference between two connectivity distributions. Stats are
% performed across the 5th dimension (boostrap samples or subjects)


connmethods = setdiff(fieldnames(PConn),{'winCenterTimes','erWinCenterTimes','freqs'});

for m=1:length(connmethods)
    
    if length(PConn)==1
        % univariate test
        [statval, df, Stats.(connmethods{m}).pval] = statcond( PConn.(connmethods{m}) , 'mode','param' );
    else
        % difference test
        [statval, df, Stats.(connmethods{m}).pval] = statcond( {PConn.(connmethods{m})} , 'mode','param');
    end
    
    
    switch multipleComparisonCorrection
        case 'fdr'
            Stats.(connmethods{m}).pval = fdr(Stats.(connmethods{m}).pval);
        case 'bonferonni'
            Stats.(connmethods{m}).pval = Stats.(connmethods{m}).pval./numel(Stats.(connmethods{m}).pval);
    end
end

