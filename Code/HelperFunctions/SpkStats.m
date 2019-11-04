function O = SpkStats(D)
ispi = diff(D);
O = [mean(ispi) std(ispi) std(ispi)/mean(ispi)];