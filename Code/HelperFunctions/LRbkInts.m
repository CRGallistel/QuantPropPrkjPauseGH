function O = LRbkInts(D)
O=double.empty(0,3); %initializing output
if length(D)<5 || min(D(:,1))>-3
    return
else
    [FO G] = fit(D(:,1),D(:,2),'poly1');
    CI = confint(FO);
    O = [FO.p1 CI(1,1) G.rsquare]; % slope of regression, lower limit on
    % slope confidence interval; amt variance accounted for
end