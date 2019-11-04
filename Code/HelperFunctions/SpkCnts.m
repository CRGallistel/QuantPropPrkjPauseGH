function O = SpkCnts(T)
O = [];
if isempty (T) || T(end)-T(1)<2 % only continue if spike train at least 2s long
    return
end
LV1 = T>T(1) & T<=T(1)+1;
LV2 = T>T(end)-1;
O = [sum(LV1)+1 sum(LV2)];