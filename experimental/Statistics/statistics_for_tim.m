% parameters to the pop_statparams and std_stat functions

% The format of each line is { key GUI_txt } default [range] help_message

{ 'condstats'  'Indep. var. 1 statistics' } 'off'    { 'on' 'off' } 'Compute (or not) statistics for the second independent variable'
{ 'groupstats' 'Indep. var. 2 statistics' } 'off'    { 'on' 'off' } 'Compute (or not) statistics for the first independent variable'
{ 'mode'       'Toolbox for computing statistics' } 'eeglab' { 'eeglab' 'fieldtrip' } 'Statistical framework to use. ''eeglab'' uses EEGLAB statistical functions and ''fieldtrip'' uses Fieldtrip statistical funcitons'
% eeglab options
	{ 'method'   'Type of statistics'       } 'parametric' { 'parametric' 'permutation' 'bootstrap' } 'Type of statistics to use. ''parametric'' uses parametric statistics. ''permutation'' uses permutation statistics and ''bootstrap'' uses bootstrap statistics'
        % only if 'permutation' or 'bootstrap'
        { 'naccu'    'Number of accumulation'   } []       [10 Inf]       'Number of surrogate averages (naccu) fo accumulate when computing permutation-based statistics. For example, to test p<0.01 use naccu>=200; for p<0.001, use naccu>=2000. If a non-NaN ''alpha'' is set (see below) and ''naccu'' is too low, it will be automatically increased. An empty entry means that the number of permutation will be computed automatically.'
    % common to 'parametric', 'permutation', 'bootstrap' (I know you told me to duplicate them but I thought it would be simplier if I was writing it in this way) 
    { 'mcorrect' 'Correction for mult. comparisons' } 'none' { 'none' 'fdr' } 'When the ''fdr'' option is selected, we will apply the False Discovery Rate method for correcting for multiple comparisons'
    { 'alpha'    'Alpha'                    } NaN      []             'Optional statistical threshold for computing p-value'
% fieldtrip options
   { 'fieldtripmethod' 'Type of statistics'       } 'analytic' { 'analytic' 'montecarlo' } 'Type of statistics to use. ''analytic'' uses parametric statistics. ''montecarlo'' uses permutation statistics.'
       % if analitc is selected
       { 'fieldtripalpha'  'Alpha'                    } NaN 'alpha' 'Optional statistical threshold for computing p-value'ture for cluster correction method, see function std_prepare_neighbors for more information.'
       { 'fieldtripmcorrect' 'Correction for mult. comparisons' } 'no' { 'no' 'bonferoni' 'holms' 'fdr' } 'Fieldtrip correction for mutliple comparison method. See help ft_statistics_analytics for more information. Default is ''no''.'
       % if montecarlo is selected
       { 'fieldtripalpha'  'Alpha'                    } NaN 'alpha' 'Optional statistical threshold for computing p-value'ture for cluster correction method, see function std_prepare_neighbors for more information.'
       { 'fieldtripnaccu'  'Number of accumulation'   } [] { [10 Inf] 'all' } 'Number of surrogate averages (naccu) fo accumulate when computing montecarlo-based statistics. This correspond to the Fieldtrip parameter ''numrandomization''. An empty entry means that the number of permutation will be computed automatically. Using ''all'' allows to perform all possible permutations.'
       { 'fieldtripmcorrect' 'Correction for mult. comparisons' } 'no' { 'no' 'bonferoni' 'holms' 'fdr' 'max' 'cluster' } 'Fieldtrip correction for mutliple comparison method. See help ft_statistics_montecarlo for more information. Default is ''no''.'
       { 'fieldtripclusterparam' 'Cluster ''key'', val parameters' } '''clusterstatistic','maxsum''' {} 'String or cell array for optional parameters for cluster correction method, see function ft_statistics_montecarlo for more information.'
       { 'fieldtripchannelneighbor' 'Channel neighbors ''key'', val parameters' } '''method','triangulation''' {} 'Fieldtrip channel neighbour structure. See help ft_statistics_montecarlo for more information.'

% execute pop_statparams to modify STUDY structure
% STUDY = pop_statparams(STUDY, ''key'', val ...)