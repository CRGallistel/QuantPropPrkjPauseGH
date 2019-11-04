function [Pre,CS,Post,E]=PrbDst(tsd,Edges,midITI,CSon,CSoff)
% computes three empirical discrete probability distribution functions, one 
% from the preCS data, one from the CS data, and one from the postCS data, 
% using the vector Edges to define the intervals
%%
mid = find(tsd(:,2)==midITI); % rows where trials notionally begin & end
CSb=find(tsd(:,2)==CSon); % CSon rows
CSe=find(tsd(:,2)==CSoff); % CSoff rows
PreIspI = [];
CSispi = [];
PostIspI = [];
%
for t=1:length(CSb) % stepping through the 20 trials
    ispi1 = diff(tsd(mid(t)+1:CSb(t)-1,1)); % pre ispki's
    ispi2 = diff(tsd(CSe(t):mid(t+1)-1,1)); % post ispki's
    ispi3 = diff(tsd(CSb(t):CSe(t),1));
    PreIspI(end+1:end+length(ispi1),1)=ispi1;
    CSispi(end+1:end+length(ispi3),1) = ispi3;
    PostIspI(end+1:end+length(ispi2),1)=ispi2;
end
Pre = histc(PreIspI,Edges);
CS = histc(CSispi,Edges);
Post = histc(PostIspI,Edges);
Pre = Pre/sum(Pre); % normalized
CS = CS/sum(CS); % normalized
Post = Post/sum(Post); % normalized
E=Edges;

