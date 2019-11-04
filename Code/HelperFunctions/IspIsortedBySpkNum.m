function O = IspIsortedBySpkNum(A,HSN)
% O is an array with HSN columns and as many rows as there were trials
% Each column, c, gives the interspike interval commenced by spike#c in
% Trial r
%
% A is a 4-col array: 1st col is spike # within a trial; 2nd col is spike 
% time; 3rd col is interspike interval; 4th col is Trial number

O = nan(A(end),HSN); % Initializing output array; there are as many rows as
% trials and HSN columns, where HSN is the highest spike count to be
% considered
for sp=1:HSN % stepping through successive spikes (columns of O)
    for t = 1:A(end) % stepping through the trials (rows of O)
        LVt = A(:,4)==t; % flags rows in A belonging to Trial t
        nsp = sum(LVt); % # of spikes in Trial t
        if nsp < sp % spike count in Trial t < nsp
            continue % go on to next trial
        else % at least nsp spikes in Trial t
            LVsp = A(:,1)==sp; % flags rows with spike count sp
            try
                O(t,sp) = A(LVt&LVsp,3); % the interspike interval commenced by
                % spike # sp in Trial t
            catch ME
                keyboard
            end
        end
    end
end