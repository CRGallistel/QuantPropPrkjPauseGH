function O = IspkInts(tsd)
if isempty(tsd)
    O=[];
else
    O = [tsd(1,1);diff(tsd(:,1))];
    % latency to 1st spike measured from last pulse in 1st US burst and any
    % subsequent inter-spike intervals
end