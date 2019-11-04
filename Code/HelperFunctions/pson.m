function OT = pson(tsd,pdf)
bins = -.3:.001:pdf(end-1,1); % bins for binary vector
LVs = tsd(:,2)==40; % flags spikes
ron = find(tsd(:,2)==80); % rows where trials start
roff = find(tsd(:,2)==90); % rows where trials end
OT = nan(length(ron),3);
pdf(end,:)=[]; % deleting last entry, which is always 0
for r = 1:length(roff)
    %%
    Dd=tsd(ron(r):roff(r),1); % times btw trial start & trial end
    LVd=LVs(ron(r):roff(r)); % flags spike times in this same stretch
    Ntst=histc(Dd(LVd),bins); % counts spike times into successive
    %% 1ms wide bins btw -.3 and end of CS
    bv = [bins' Ntst]; % binarized spike vector
    bv(bv(:,2)>1,2)=1; % eliminating the (rare) integers>1
    LVcs = pdf(:,1)>0; % flags portions of pdf within CS
    Mn = min(pdf(LVcs,2)); % smallest p_s within CS
    alphaB = 1; % setting hyperparameters
    betaB = ceil((1-Mn)/Mn); % setting hyperparameters
    alphaA = 1; % setting hyperparameters
    AvCSp = mean(pdf(~LVcs,2)); % setting hyperparameters
    betaA = ceil((1-AvCSp)/AvCSp);% setting hyperparameters
    [CP,Odds] = BernCP(flipud(bv(:,2)),alphaB,betaB,alphaA,betaA,1/length(bv));
    W = log10(Odds);
    cp = bv(length(bins)-CP,1);
    LVps = bv(:,1)>cp;
    lamPs = sum(bv(LVps,2))/(pdf(end,1)-cp);
    lamPre = sum(bv(~LVps,2))/(cp-pdf(1,1));
    OT(r,:) = [cp lamPre-lamPs W];    
end    
