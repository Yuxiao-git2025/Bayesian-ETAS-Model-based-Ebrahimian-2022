function Bayesian_CorrMapping(samples)
[m,n]=size(samples);
if m==8 && n>=1
    fprintf('# 8 parameter have been loaded >> \n');
    samples=samples'; % Transform col
    samples=samples(:,[1,3:end]);
else
    error('Please check dimension \n');
end
ncolumn = size(samples, 2);
figure('Color','w');
BigAx = axes('Position', [0.1, 0.1, 0.8, 0.8], 'Visible', 'off','Color','w');
hold(BigAx, 'on');

AX = gobjects(ncolumn);
h = gobjects(ncolumn, 1);
s = gobjects(ncolumn);

gap = 0.05;
width = (1 - gap*(ncolumn+1)) / ncolumn;

for i = 1:ncolumn
    pos = [gap*i + width*(i-1), gap*i + width*(i-1), width, width];
    AX(i,i) = axes('Position', pos);

%     h(i) = histogram(AX(i,i), X(:,i),'Normalization','pdf');
%     [cp1, cp2] = calculate_CI(X(:,i), 2, 98);
%     hold(AX(i,i),"on");
%     xline(cp1,'Color',[1 0 0 ],'LineWidth',1.8,'LineStyle','--');
%     xline(cp2,'Color',[1 0 0],'LineWidth',1.8,'LineStyle','--');
%     h(i).FaceColor = [1 0 0];
%     h(i).FaceAlpha = 0.5;
%     h(i).EdgeColor = 'k';
%     h(i).EdgeAlpha = 0.8;
%     h(i).LineWidth = 1.2;
% 除了直方图，还有核密度图（KDE）以及等值线图可以使用
    SmpKDE_Plot(samples(:,i));
    hold(AX(i,i),"on");
    [cp1, cp2] = calculate_CI1(samples(:,i), 2, 98);
%     [cp3, cp4] = calculate_CI2(X(:,i), 16, 84);
    xline(cp1,'Color','k','LineWidth',1.8,'LineStyle','--');
    xline(cp2,'Color','k','LineWidth',1.8,'LineStyle','--');
end

% 创建三角散点图
for i = 1:ncolumn
    for j = (i+1):ncolumn
        pos = [gap*j + width*(j-1), gap*i + width*(i-1), width, width];
        AX(i,j) = axes('Position', pos);
        s(i,j)=scatter(AX(i,j), samples(:,j), samples(:,i),6);
        s(i,j).Marker = "o";
        %             s(i,j).MarkerSize=6;
        s(i,j).MarkerFaceColor = [0.5020 0.5020 0.5020];
        s(i,j).MarkerEdgeColor = [0.6510 0.6510 0.6510];
        % 计算相关系数
        r = corr(samples(:,j), samples(:,i));
        text(AX(i,j), 0.05, 0.90, sprintf('r=%.4f', r), ...
            'Units', 'normalized', ...
            'FontSize', 15, ...
            'BackgroundColor', 'none','Color','k');
%         set(AX(i,j), 'XTick', [], 'YTick', []);

    end
end

for i = 1:ncolumn
    for j = i:ncolumn
        AX(i,j).LineWidth = 1.2;
        AX(i,j).FontWeight = "bold";
        AX(i,j).MinorGridLineStyle = "none";
        AX(i,j).GridLineStyle = "none";
        AX(i,j).Box = 'on';
    end
end

end