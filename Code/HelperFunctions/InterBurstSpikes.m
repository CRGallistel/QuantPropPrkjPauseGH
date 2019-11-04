function O = InterBurstSpikes(tsd,spikecode)
starttime = tsd(1,1);
LVspk = tsd(:,2)==spikecode; % flags spikes
if sum(LVspk)<1 % no spikes
    O=[];
else
    O = tsd(LVspk,1) - starttime - .008;
    % last pulse in 1st burst of stimulation occurs 8 ms after US onset
    O = O(O>0); % deleting negative spike latencies, which come from spikes
        % recorded during stimulation burst
end