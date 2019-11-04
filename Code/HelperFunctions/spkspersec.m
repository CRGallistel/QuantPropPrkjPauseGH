function r=spkspersec(spktms,td)
if isempty(spktms)
    r=0;
    S=evalin('caller','sub');s=evalin('caller','ses');t=evalin('caller','tri');
    fprintf('\nNo spikes in S%d,s%d,t%d\n',S,s,t)
else
    N = numel(spktms);
    D=td-spktms(1); % because recorder not turned on for some seconds
    r=N/D;
end