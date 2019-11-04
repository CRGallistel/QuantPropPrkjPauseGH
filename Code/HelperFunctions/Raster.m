function Raster(tsd,pre,post)
LVcs = tsd(:,2)==50; % flags CS onsets
LVus = tsd(:,2)==70; % flags US onsets
Db = [tsd(LVcs,1)-pre 45*ones(size(tsd(LVcs,1)))]; % to be inserted onset data
De = [tsd(LVus,1)+post 85*ones(size(tsd(LVus,1)))];% to be inserted offset data
tsd = sortrows([tsd;Db;De]);
TSraster(tsd,[45 85],[50 0;70 0;40 0],['b+';'b+';'k.'])
% 50 0 = CSon, marked by green pluses; 70 0 = USon, marked by blue pluses;
% 60 0 = CSoff, marked by red inverted triangles; 40 0 = spikes, marked by
% black dots
xlim([0 .1])
plot([on on],ylim,'g--',[off off],ylim,'r--'))
title(['Cell' num2str(evalin('caller','sub'))])
