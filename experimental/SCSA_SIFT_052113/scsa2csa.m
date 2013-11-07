function W = csa2scsa(B, H)
%SCSA2CSA converts the demixing matrix B and MVAR coefficients H of SCSA
% into convolutive parameters W of CSA 
%
%   input  : B - M x M demixing matrix
%            H - M x M x K MVAR coefficient matrices
%   outputs: W - M x M x K+1 deconvolution matrix 

W(:, :, 1) = B;
for k = 1:size(H, 3)
  W(:, :, k+1) = -H(:, :, k)*B;
end

end

