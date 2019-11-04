function hst = PeriStimHist(D,bw,tr)
% Computes a peristimulus spike histogram from 2-col binary vector (col 1 =
% times, col 2 flags spikes) for bin width bw over time range tr
% (e.g., [-.2 1.1]) and generates the histogram
LVs = D(:,2)>0; % flags rows with spikes
LV = D(:,1)>tr(1) & D(:,1)<tr(2) & LVs; % flags rows w spikes within the
% the specified time range
edges = tr(1):bw:tr(2);
SpkCnts= histc(D(LV,1),edges);
figure
bar(edges,SpkCnts,'histc')
hst = [edges' SpkCnts];