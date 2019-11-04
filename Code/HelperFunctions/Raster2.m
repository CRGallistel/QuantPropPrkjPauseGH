function Raster2(tsd,pre,post,ax)
LVcs = tsd(:,2)==50; % flags CS onsets
LVus = tsd(:,2)==70; % flags US onsets
Db = [tsd(LVcs,1)-pre 45*ones(size(tsd(LVcs,1)))]; % to be inserted onset data
De = [tsd(LVus,1)+post 85*ones(size(tsd(LVus,1)))];% to be inserted offset data
tsd = sortrows([tsd;Db;De]);
[~,b] = TSmatch(tsd,{[80 90]});
for tt = 1:length(b)
    Dt = tsd(b{tt}(1):b{tt}(2),:); % data for this trial
    on = Dt(find(Dt(:,2)==50,1),1);
    Dt(:,1) = Dt(:,1)-on; % refered to a 0 at CS onset
    LVs = Dt(:,2)==40; % flags spikes
    plot(ax,Dt(LVs,1),tt*ones(sum(LVs),1),'k.')
    hold on
end        
off = Dt(find(Dt(:,2)==70,1),1);
xlim([-.3 .6])
plot([0 0],ylim,'g--',[off off],ylim,'r--')
title(['Cell' num2str(evalin('caller','sub'))])
