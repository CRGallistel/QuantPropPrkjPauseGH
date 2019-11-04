function tsd=AddUSoff(tsd,USdurs,ONcode,OFFcode)
S = evalin('caller','Experiment.Subject(sub).SubId'); % current subject
LVon = tsd(:,2)==ONcode; % flags USons
dur = USdurs(USdurs(:,1)==S,2)/1000; % duration of US for current subject, in ms
tsd = [tsd;[tsd(LVon,1)+dur OFFcode*ones(sum(LVon),1)]]; % tsd(LVon,1)+dur is a
    % the vector of session times for the US offs; OFFcode*ones(sum(LVon),1) is
    % a vector of the same length whose elements are the event code for USoff
tsd = sortrows(tsd);