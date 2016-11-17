function [B H] = csa2scsa(W)
%CSA2SCSA converts convolutive parameters W of CSA into demixing matrix B
%and MVAR coefficients H of SCSA
%
%   input  : W - M x M x K+1 deconvolution matrix 
%   outputs: B - M x M demixing matrix
%            H - M x M x K MVAR coefficient matrices

B = W(:, :, 1);
for k = 1:size(W, 3)-1
  H(:, :, k) = -W(:, :, k+1)/B;
end

end

