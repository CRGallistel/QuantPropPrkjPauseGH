function [Mx,Min] = MxMin(tsd)
ispkis = diff(tsd(:,1));
Mx = max(ispkis);
Min = min(ispkis);