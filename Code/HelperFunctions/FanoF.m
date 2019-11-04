function O = FanoF(T)
%%
edges = T(1):.1:T(end);
N = histc(T,edges);
lam = numel(T)/(T(end)-T(1)); % estimated firing rate
n = numel(edges); % number of samples
FF = var(N)/mean(N); % estimated Fano Factor
if FF<1
    p = gamcdf(FF,(n-1)/2,2/(n-1));
else
    p = 1-gamcdf(FF,(n-1)/2,2/(n-1));
end    
O = [FF lam n p]; % [Fano Factor lambda n (p if 1)]