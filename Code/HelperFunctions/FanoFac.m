function O = FanoFac(Cnts)
if isempty(Cnts) || size(Cnts,1)<10
    O=[];
else
    O = [var(Cnts(:,1))/mean(Cnts(:,1)) var(Cnts(:,2))/mean(Cnts(:,2)) size(Cnts,1)];
end % 2 Fano Factors and the n (sample size)