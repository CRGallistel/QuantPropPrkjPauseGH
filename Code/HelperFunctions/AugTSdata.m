function tsdA = AugTSdata(tsd,CSdur)
CSdur = CSdur/1000; % CS duration passed in in milliseconds, but time
% stamps are in seconds
Mtms=tsd(tsd(:,2)==20,1); % session times for the "mark" pseudo-event
CSonTms = Mtms+.2; % the mark precedes CS onset by .2s
CSoffTms = CSonTms+CSdur;
ITIs = CSonTms(2:end)-CSoffTms(1:end-1); % the intertrial intervals
MidITItms = [tsd(1,1)-.001;CSoffTms(1:end-1)+round(1000*ITIs/2)/1000;tsd(end,1)+.001];
% The points midway between the CSoffs and the CSons, that is, in the
% middle of the ITIs--with first time 1 ms before the first spike and a
% final time 1 ms after the last spike
MidITIs = [MidITItms 60*ones(length(MidITItms),1)]; % adding event codes 
CSons = [CSonTms 30*ones(length(CSonTms),1)]; 
CSoffs=[CSoffTms 50*ones(length(CSoffTms),1)];
tsdA = sortrows([tsd;CSons;CSoffs;MidITIs]);