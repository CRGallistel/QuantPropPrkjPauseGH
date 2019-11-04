function [O1,O2] = SpkTmsIspI(spktms)
O1 = spktms(1); % latency of 1st spike following US offset
O2 = [(1:size(spktms,1)-1)' spktms(1:end-1) diff(spktms)];