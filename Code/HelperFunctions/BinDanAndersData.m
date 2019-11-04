function BV = BinDanAndersData(tsd)
rTon = find(tsd(:,2)==80); % row #s for trial on
rToff = find(tsd(:,2)==90); % row #s for trial off
rcs = find(tsd(:,2)==50); % row #s for CS on
BV = tsd;
for r = 1:length(rToff)
    BV(rTon(r):rToff(r),1) = BV(rTon(r):rToff(r),1)-BV(rcs(r),1);
end
    