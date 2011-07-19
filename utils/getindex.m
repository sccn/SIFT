% function [idx value] = getindex(A,els)
% return the indices of values in vector A that are nearest to those in vector els
% outputs are the length of els
%
% Tim Mullen 2010, SCCN/INC, UCSD

function [idx value] = getindex(A,els)

if isempty(els)
    idx = [1 length(A)];
    return;
end

% if any(els<min(A)) || any(els>max(A))
%     fprintf('warning! some elements of els are not within range of A -- selecting nearest indices\n');
% end

L = length(els);
[value idx] = deal(zeros(1,L));
for i=1:L
    [value(i) idx(i)] = min(abs(A-els(i)));
end
