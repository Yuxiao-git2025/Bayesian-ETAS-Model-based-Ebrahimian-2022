fore=poissrnd(12.22,1e4,1);
obs=16;
[d1,d2]=N_test(fore,obs);
fprintf('d1= %.4f d2= %.4f Same= %d \n',d1,d2,length(find(fore==obs)));

histogram(fore,'LineWidth',1,'FaceAlpha',0.6,'FaceColor','r','Normalization', ...
    'pdf');
% SmpKDE_Plot(fore);
hold on;
xline(obs,'LineWidth',2,'Color','k','LineStyle','--');
Fun_defaultAxes;
xlabel('Nums'); ylabel('PDF');
%%
N_Poitest(12.22,16);
