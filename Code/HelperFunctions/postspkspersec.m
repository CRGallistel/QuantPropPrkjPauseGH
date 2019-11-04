function r=postspkspersec(spktms)
if isempty(spktms)
    r = 0;
else
    r=numel(spktms)/spktms(end); % spktms(end) gives duration of trial up
    % to point where recording stopped
end