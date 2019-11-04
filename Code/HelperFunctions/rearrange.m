function O = rearrange(ints,cols)
O = nan(1,cols);
if isempty(ints)
    return
elseif length(ints)>cols
    c = cols;
else
    c = length(ints);
end
O(1,1:c) = ints(1:c);