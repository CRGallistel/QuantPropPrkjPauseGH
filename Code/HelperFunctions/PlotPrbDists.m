function PlotPrbDists(preD,csD,postD,Hv)
axes(Hv(1))
h=bar(preD,1); % bar plot of prob dist in upper right panel
ylim([0 .7])
set(h,'FaceColor',[1 1 1],'LineWidth',2)
ylabel('Probability','FontSize',12)
set(gca,'XTick',.5:1:8.5,'XTickLabel',{'0' '5' '10' '15' '20' '25' '30' '40' '>40'})
text(.2,.6,'PreCS')

axes(Hv(2))
h=bar(csD,1);
ylim([0 .7])
set(h,'FaceColor',[1 1 1],'LineWidth',2)
ylabel('Probability','FontSize',12)
set(gca,'XTick',.5:1:8.5,'XTickLabel',{'0' '5' '10' '15' '20' '25' '30' '40' '>40'})
text(.2,.6,'CS')

axes(Hv(3))
h=bar(postD,1);
ylim([0 .7])
set(h,'FaceColor',[1 1 1],'LineWidth',2)
ylabel('Probability','FontSize',12)
set(gca,'XTick',.5:1:8.5,'XTickLabel',{'0' '5' '10' '15' '20' '25' '30' '40' '>40'})
xlabel('InterspikeInterval (ms)','FontSize',12)
text(.2,.6,'PostCS')