function tsdA = AddEvents(tsd,ISI)
% tsd is the original data; ISI is the CS-US interval from the Phase field
LVto = diff(tsd(:,1))>5; % Logical vector flagging 1st spike in each trial,
% but displaced one event too early because of differencing
LVtoff = [LVto;true]; % logical vector flagging last spike in
% each trial
LVto = [true;LVto]; % adding onset of 1st trial to logical vector and
% displacing flags one event forward to where they belong
TrlOnTms = tsd(LVto,1)-.001; % trial onsets set to 1 ms before 1st spike time
TrlOnD = [TrlOnTms 80*ones(size(TrlOnTms))]; % trial onset data; 80 is the
% event code for trial on
TrlOffTms = tsd(LVtoff,1)+.001; % trial onsets set to 1 ms after last spike
TrlOffD = [TrlOffTms 90*ones(size(TrlOffTms))]; % trial offset data; 90 is
% the event code for trial off

CSonTms = nan(size(TrlOnTms));
USonTms = nan(size(TrlOnTms));
CSoffTms = nan(size(TrlOnTms));
i = 1;
for t = TrlOnTms'
    CSonTms(i) = tsd(find(tsd(:,1)>t & tsd(:,2)==30,1),1)-.001; % subtract
    % 1 ms from 1st pulse time to insure onset precedes 1st pulse
    LV = tsd(:,1)>CSonTms(i) & tsd(:,1)<CSonTms(i)+1; % compassing current CS
    d = tsd(LV,:); 
    r = find(d(:,2)==30,1,'last'); % last pulse of current CS  
    CSoffTms(i) = d(r,1)+.001; % add 1 ms to insure that this event occurs
    % after the final pulse
    USonTms(i) = CSonTms(i)+ISI/1000;
    i=i+1;
end
    
CSonD = [CSonTms 50*ones(size(CSonTms))]; % CS on data; 50 is the event
% code for CS on
CSoffD = [CSoffTms 60*ones(size(CSoffTms))]; % CS off data; 60 is the event
% code for CS off
USonD = [USonTms 70*ones(size(USonTms))]; % US on data; 70 is the event code

tsdA = sortrows([tsd;TrlOnD;TrlOffD;CSonD;CSoffD;USonD]);