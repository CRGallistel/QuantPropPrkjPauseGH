function [On,Won,Off,Woff,lb,ld,la,MxIpI,LatMx] = BinPsOnOff2(D,dur)
% Computes the onsets and offsets of the pauses
D(D(:,2)>1,2)=1; % This corrects a single bin with a count of 2 in the
% BinSpkVecS field for S10 in the parallel fiber data. That triggers a
% crash of this function when BernCP is called
r = find(D(:,1)>-.0005 & D(:,1)<.0005); % vector of row numbers for CS onset
NmTrls = length(r);
On = nan(NmTrls,1); % initializing On vector
Off = nan(NmTrls,1); % initializing Off vector
Won = zeros(NmTrls,1); % initializing vector for weight of the evidence for
% pause onset
Woff = zeros(NmTrls,1); % inti
lb = nan(NmTrls,1);
ld = nan(NmTrls,1);
la = nan(NmTrls,1);
MxIpI = nan(NmTrls,1);
LatMx = nan(NmTrls,1);
LVs = D(:,2)>0; % flags rows with spikes
LV = D(:,1)>-dur & D(:,1)<3*dur & LVs; % flags rows w spikes btw -dur and
% 3dur
edges = -dur:.2*dur:3*dur;
SpkCnts= histc(D(LV,1),edges);
[Nsdur,loc] = min(SpkCnts(1:end-1)); % lowest spike count within an interval
% equal to 1/5th of CS duration and the bin that contains it 
% 
LVpre = D(:,1)<0 & D(:,1)>-2*dur; % flags pre CS rows over an immediately
% preceding interval equal in duration to the CS
LVpst = D(:,1)>dur & D(:,1)<2*dur; % flags post CS rows over an immediately
% following interval equal in duration to the CS
LVb = D(:,1)>edges(loc); % flags rows greater than lower edge of lowest-count 
%  bin in the histogram
LVe = D(:,1)<edges(loc)+.2*dur; % flags rows less than the upper edge of 
% the lowest-count bin in the histogram
LVmid = LVb&LVe;
Npre = sum(LVpre); % number of pretrial bins 
Nspre = sum(LVpre&LVs); % number of pretrial spikes
if Nspre<1
    Nspre = 1;
end
HypPre = [1 round((Npre-Nspre)/Nspre)]; % hyperparameters of beta prior on
    % pre-CS probability of a spike in a 1 ms bin
if isnan(HypPre(2)) || ~beta(HypPre(1),HypPre(2))>0 
    keyboard
end
Ndur = sum(LVmid); % number of binarize bins in lowest bin of histogram 
if Nsdur<1
    Nsdur = 1;
end
HypDur = [1 round((Ndur-Nsdur)/Nsdur)]; % hyperparameters on the probability
% of a spike in any one ms in the middle 40% of the CS bins
if isnan(HypDur(2)) || ~beta(HypDur(1),HypDur(2))>0
    keyboard
end
Npst = sum(LVpst); % # of post-CS bins
Nspst = sum(LVpst&LVs); % # of spikes in post-CS bins
if Nspst<1
    Nspst = 1;
end
HypPst = [1 round((Npst-Nspst)/Nspst)];
if isnan(HypPst(2)) || ~beta(HypPst(1),HypPst(2))>0
    keyboard
end
nbpers = 1000; % number of bins per second
%
NbinsBack = floor(dur*1000); % number of bins prior to CS onset for pause On
NbinsFwd = ceil(edges(loc+1)*1000); % positions start of backward search for CS
% onset at forward edge of minimum bin in peristimulus histogram
NbinsToStrt = floor(edges(loc)*1000); % positions start of forward search for pause
% offset at backward edge of minimum bin in peristimulus histogram
NbinsPst = ceil(3000*dur); % number of bins to end of data for pause Off
%%
for tt = 1:NmTrls
    %%
    if r(tt)-NbinsBack < 1 % happens in some data sets
        NbinsBack=r(tt)-1;
    end
    Don = flipud(D(r(tt)-NbinsBack:r(tt)+NbinsFwd,2)); % the binary spike
    % vector running backwards in time from 70% of the way through the CS
    % (whence the 700 times a time equal to CS onset minus CS duration in
    % seconds)
    try
    [CPon,OddsD] = BernCP_local(Don,HypDur(1),HypDur(2),HypPre(1),HypPre(2),1/length(Don));
    % the onset cp in # of bins counting backwards
    catch ME
        keyboard
    end
    %
    if r(tt)+NbinsPst<length(D) 
        Doff = D(r(tt)+NbinsToStrt:r(tt)+NbinsPst,2); % data for finding 
        % pause offset run from backward edge of minimum bin to three
        % times the duration of the CS
    else % don't go past end of data
        Doff = D(r(tt)+NbinsToStrt:end,2);
    end
    
    [CPoff,OddsU] = BernCP_local(Doff,HypDur(1),HypDur(2),HypPst(1),HypPst(2),1/length(Doff));
%     if tt == 12
%         keyboard
%     end
    if sum(Don(1:CPon))>sum(Don(CPon:end))...
            ||sum(Doff(1:CPoff))>sum(Doff(CPoff:end))...
            ||sum(Don)+sum(Doff)<10
        % more spikes after CPon than before (remember: it's backwards here)
        % or during CP than after or too few spikes to make valid decision
        continue % don't record pause onsets & offsets
    end
%     if tt == 15;keyboard;end
    On(tt) = edges(loc+1)-CPon/1000; % the estimated pause onset time in
    % seconds relative to the CS onset time (remember CP is backwards from
    % forward edge of minimum bin
    Won(tt) = log10(OddsD); % the weight of the evidence for the onset
    Off(tt) = edges(loc)+CPoff/1000; % the estimated pause offset time in
    % seconds relative to a 0 at CS onset
    Woff(tt) = log10(OddsU); % weight of the evidence for the pause offset
    lb(tt) = sum(D(r(tt)-NbinsBack:r(tt),2))/dur;
    ld(tt) = sum(D(r(tt)+NbinsToStrt:r(tt)+NbinsFwd,2))/dur;
    la(tt) = sum(D(r(tt)+dur*nbpers:r(tt)+2*dur*nbpers,2))/dur;
    % Following code added 04/23/18
    LVtt = D(:,3)==tt; % flags data from current trial
    LVps = D(:,1)>On(tt) & D(:,1)<Off(tt); % flags times where trial time
    % falls within the PsOn and PsOff times for this trial (but across all
    % trials)
    SpkTms = D(LVtt&LVps&LVs,1); % spike times within the pause boundaries on
    % THIS trial
    if isempty(SpkTms)
        MxIpI(tt) = Off(tt)-On(tt);
        LatMx(tt) = 0;
    else
        Tms = [On(tt);SpkTms;Off(tt)];
        IspIs = diff(Tms); % intervals
        [MxIpI(tt),indx] = max(IspIs);
        LatMx(tt) = Tms(indx)-On(tt);
        if LatMx(tt)<0;LatMx(tt)=0;end
    end        
end

function [CP,Odds,LkFun,CL,CPmean] = BernCP_local(Data,AlphaB,BetaB,AlphaA,BetaA,pc)
% Computes maximally likely change point (CP) and odds in favor of a change
% given a binary data string, Data, the parameters of the beta prior
% on the value of p before the change (AlphaB & BetaB), the parameters
% of the beta prior on the p after the change (AlphaA,BetaA), and the
% prior probability of a change after any given trial (pc). This is called
% from Line 515 with AlphaB & BetaB evolving and from Line 597 with them
% set to their initial agnostic values (.5 or 1)
%
%Syntax [CP,Odds,LkFun,CL] = BernCP(Data,AlphaB,BetaB,AlphaA,BetaA,pc)
%
% CL is a cell array giving lower (L) and upper (U) relative likelihood
% limits as follows:
% CL = {[.01L] [.05L] [.1L] [.2L] [.5L] [.5U] [.2U] [.1U] [.05U] [.01U]}
% CPmean is the mean of the marginal likelihood distribution on CP. The CP
% output is this value rounded to the nearest integer

Data = reshape(Data,length(Data),1); % making sure Data is a column vector

Data = Data+0; % converts logic vectors to 0's & 1's

if ~isempty(setdiff(Data,[0 1]))
    disp('Error: Input is not a binary vector')
    return
end

CP = nan; % intializing

Odds = nan;% intializing

LkFun = zeros(length(Data),1);% intializing

CL(1,1:10) = {[]};% intializing

CPmean = nan;

if nargin<6
    display({'Not enough input arguments;';'6 required: ';...
        'Vector of binary data';...
        'Hyperparameter alpha for prior on p before change';...
        'Hyperparameter beta for prior on p before change';...
        'Hyperparameter alpha for prior on p after change';...
        'Hyperparameter beta for prior on p after change';...
        'Probability of a change after any one trial'})
    return
end
    
if ~(beta(AlphaB,BetaB)>0)||~(beta(AlphaA,BetaA)>0)
    
    display('Prior too strong; makes beta(a,b)=0')
    
    return
    
end
%%

p = .01:.01:.99; % row vector of p values

dp = p(2)-p(1); % delta p

L = (1:length(Data))'; % cumulative number of observations (col vector)

Sf = cumsum(Data); % cumulative successes up to and including the nth datum
% (col vector w length = length(Data))

Ff = L - Sf; % cumulative failures up to and including the nth datum
% (col vector w length = the length of the data)

Sr = sum(Data) - Sf; % succeses after nth datum (col vector w length =
% the length of the data))

Lr = ((L(end)-1):-1:0)'; % number of observations after the nth datum
% (length = the length of the data)

Fr = Lr - Sr; % failures after the nth datum (col vector, w length =
% the length of the data)

%% The following relies on the analytic result that
% given a beta prior with parameters, alpha & beta, and observed
% numbers of successes and failures, nS & nF, the posterior marginal
% likelihood of a change hypotheses, given the data and the assumption of
% a beta prior with parameters alpha & beta) is beta(alpha+nS,beta+nF).

LgPstLkf = log(beta(AlphaB + Sf, BetaB + Ff)); % col vector
% giving the log likelihood for the data up to and including datum(n) under
% the prior distribution on p specified by hyperparameters AlphaB and BetaB

LgPstLkr = log(beta(AlphaA + Sr, BetaA + Fr)); % col vector
% giving the log likelihood for the data after datum(n) under the prior
% distribution specified by hyperparameters AlphaA & BetaA


%%
LgLkFun = LgPstLkf + LgPstLkr; % the log likelihood function giving, as a
% function of cp, the log likelihood of a model that assumes a change in p 
% after observation cp.

LgLkFunN = LgLkFun + abs(max(LgLkFun)); % setting peak value of log
% likelihood function to 0, so that when it is transformed into a
% likelihood function prior to integration its values are reasonable. This
% is equivalent to rescaling the likelihood function, which has no effect
% on the relative likelihood or the expectation, which are the two
% quantities of interest

LkFun = exp(LgLkFunN); % the likelihood function
% The last value in this vector is the likelihood of a model
% that assumes a change after the end of the data, in other words, no
% change within the span of the observations. Thus, assuming a uniform
% prior on cp (change equally likely after any trial), the marginal
% likelihood of a change is the mean value of this likelihood function

CPmean = sum((1:length(Data)).*LkFun')/sum(LkFun); % expectation of the
% likelihood function (NB, not its mean value)

RelLkhd = (sum(LkFun(1:end))/(length(Data)))/LkFun(end); % the Bayes
% Factor, the ratio of the posterior likelihoods of the two models (change
% vs no-change). (sum(LkFun(1:end))/(length(Data)) is the mean likelihood,
% which, assuming a uniform prior, is the the marginal likelihood of a
% change; LkFun(end) is the likelihood of no change.
%%

pNC = 1-pc; % Probabilty that p does not change after any given trial

PriorOdds = L(end)*pc/pNC; % prior odds that there is a change
% in a sequence as long as Data = Np_c/(1-p_c)

Odds = PriorOdds*RelLkhd; % Posterior odds equal relative likelihood of the
% two models (change: no-change) times the prior odds in favor of one or
% the other

BF = RelLkhd;
%% Computing CP & relative likelihood limits on its location

CP = round(CPmean);

Mx = max(LkFun(1:end-1)); % maximum of the likelihood function for
% change points before the last datum and location of this maximum
% (the estimate of CP)

CL{1} = find(LkFun>.01*Mx,1)-1;

CL{2} = find(LkFun>.025*Mx,1)-1;

CL{3} = find(LkFun>.1*Mx,1)-1;

CL{4} = find(LkFun>.2*Mx,1)-1;

CL{5} = find(LkFun>.5*Mx,1)-1;

CL{6} = find(LkFun>.5*Mx,1,'last')+1;

CL{7} = find(LkFun>.2*Mx,1,'last')+1;

CL{8} = find(LkFun>.1*Mx,1,'last')+1;

CL{9} = find(LkFun>.025*Mx,1,'last')+1;

CL{10} = find(LkFun>.01*Mx,1,'last')+1;

