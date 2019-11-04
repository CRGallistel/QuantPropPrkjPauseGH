function Hst = PeriCShist(D,CPs,bw)
%%
S = evalin('caller','sub');
%
switch S
    case 6
        nb = 5; % number of steps back in the CPs
    case 7
        nb = 3;
    case {(1) (3) (4) (8) (9) (10)}
        nb = 1;
    otherwise
        nb = 2;
end
CPs(1,1)=1;
ron = find(D(:,2)==80); % row #s at trial starts
 % row # where clear pause has emerged
D(1:ron(CPs(end-nb,1)),:)=[]; % deleting data prior to clear pause trials
nt = CPs(end,1) - CPs(end-nb,1); % number of trials from which counts computed
LVs = D(:,2)==40; % flags spikes
EndCSUS = D(find(D(:,2)== 70,1),1)-.005; % trial time when CS terminates
LVdh = D(:,1)>-.3 & D(:,1)<EndCSUS; % flags all the stretches between -0.3s
% and the end of the CS-US interval
edges = -.3:bw:EndCSUS;
N = histc(D(LVdh&LVs,1),edges);
Hst = [edges' N/(nt*bw/.001)]; % p_s in a 1ms bin is the count divided by
% the product of the number of trials over which histogram was computed
% (nt) and the # of 1 ms bins in a histogram bin (bw/.001) 
subplot(5,2,S)
% bar(edges,N,'histc')
bar(edges,N/(nt*bw/.001),'histc')
xlim([-.3 .3])
if S>8
    xlabel('Trial Time (s)')
end
if mod(S,2)>0
    ylabel('p_s')
end