function [spikes,Mn,Mx] = Bindata2(tsd,BW,spike,CSon)
% A Benoulli vector digitization of the spike train with bins of  width BW
% Syntax     spikes = Bindata2(tsd,BW,spike,CSon)
LV = tsd(:,2)==spike; % flags spikes
LVCSon = tsd(:,2)==CSon; % flags CS on
spks = tsd(LV,1)-tsd(LVCSon,1); % spike times referenced to 0 at CS on
Mn = min(spks);
Mx = max(spks);
bins = Mn:BW:Mx; % creating the bin edges
BV = histc(spks,bins);
t = bins'+BW/2;
spikes = [t BV];