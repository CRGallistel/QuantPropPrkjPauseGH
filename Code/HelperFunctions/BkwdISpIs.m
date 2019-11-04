function O=BkwdISpIs(bsv,PsOns)
O = double.empty(0,2); % initializing output array
rws = find(bsv(:,1)>-.0005 & bsv(:,1)<.0005); % CS onset rows
rnv = (1:length(bsv))';
LVs = bsv(:,2)>0; % flags spike rows
i=1;
for br = rws' % stepping through the trials
    if PsOns(i)>0
        er = br+round(PsOns(i)/.001); % row # at pause onset
        LVt = rnv>br&rnv<er; % flags rows btw CS onset and pause onset
        spktms = bsv(LVt&LVs,1);
        n = length(spktms);
        if n>1 % more than 1 spike btw CS onset & pause onset
            O = [O;[(-1:-1:-n+1)' flipud(diff(spktms))]];
        else
            i = i+1;
            continue
        end 
    else
        i=i+1;
        continue
    end
    i=i+1;
end
    