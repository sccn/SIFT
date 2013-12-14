%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% information outflow %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmpRPDC = tmpEEG.CAT.Conn.RPDC;
    % process inflow
    for m = 1:size(tmpRPDC,1)
        tmpIcInflow = squeeze(mean(tmpRPDC(m,setdiff(1:size(tmpRPDC,1),m),:,:),2));
        tmpIcInflow = squeeze(mean(tmpIcInflow,2));
        tmpEEG.inflow(:,m) = tmpIcInflow(1:195,:);
    end
EEG.inflow = tmpEEG.inflow;


% make sure that you have loaded STUDY in which all datasets have SIFT
% measures
dipoleAndMeasure = pr.dipoleAndMeasureOfStudyCustom(STUDY, ALLEEG);

% build 'linearizedMeasure' under dipoleAndMeasure 
for i=1:length(dipoleAndMeasure.datasetIdAllConditions)   
    outflow = ALLEEG(dipoleAndMeasure.datasetIdAllConditions(i)).outflow(:,dipoleAndMeasure.numberInDataset(i));
    % outflow = outflow(goodFreqIdx);
    % outflow = log(outflow);
    % outflow = outflow - groupMeanOutflow'; % subtract mean
    outflow = bsxfun(@minus, outflow, mean(outflow)); % subtract mean
    dipoleAndMeasure.linearizedMeasure(:,i) = outflow;
end;

% compute optimal Gaussian width
[optimalGaussianWidth meanPredictionSImilarity gaussianStdValues ] = pr.find_optimal_gaussian_width(dipoleAndMeasure);

% compute group-difference
outflowProjection = pr.meanProjection(dipoleAndMeasure, dipoleAndMeasure.getPairwiseFishersZSimilarity, pr.headGrid, 'stdOfDipoleGaussian', optimalGaussianWidth);

% project measures
outflowProjection = outflowProjection.createDomain(dipoleAndMeasure, 0.9, 0.05);
save /data/cta/Anthony/Rest64_Biosemi/mpOutflowLogMeansub.mat dipoleAndMeasure outflowProjection optimalGaussianWidth
