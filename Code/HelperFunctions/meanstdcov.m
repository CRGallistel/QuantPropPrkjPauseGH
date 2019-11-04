function O = meanstdcov(D)
O = [nanmean(D) nanstd(D) nanstd(D)/nanmean(D)];