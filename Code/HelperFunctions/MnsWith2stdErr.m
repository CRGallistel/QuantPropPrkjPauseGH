function MnsWith2stdErr(A,xlm)
mu = mean(A); % bin means
stder = std(A)/sqrt(size(A,1)); % standard errors
b = 1:size(A,2);
figure
errorbar(b,mu,2*stder)
xlabel('Successive 25 ms Bins')
ylabel('Mean # Spikes +/-2 std er')
xlim(xlm)