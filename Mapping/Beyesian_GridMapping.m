function Beyesian_GridMapping(Gridfore,colmap)
global lonMin lonMax latMin latMax Xcgrid Ycgrid Ggrid
if isempty(colmap)
    colmap='jet';
end
row=length(Ycgrid);
col=length(Xcgrid);
fprintf('Grind lines: %d X %d \n',row,col);
fprintf('longitude-Bound: [%.4f\t %.4f\t] \n',lonMin,lonMax);
fprintf('latitude -Bound: [%.4f\t %.4f\t] \n',latMin,latMax);
% observe data
[~,magin,lonin,latin]=generate_Observe();
Choose=2;
switch Choose
    case 1
%%                            Drawing without map
% figure;
% Fun_defaultAxes;
% hold on;
% if length(linesx)~=size(Gridfore,2) || length(linesy)~=size(Gridfore,1)
%     fprintf('Notice the dimension of Grid is error \n');
% else
%     % forecast data
%     ima=imagesc(linesx,linesy,(Gridfore),'HandleVisibility','off');
%     colormap(slanCM(colmap));   colorbar;
% 
%     markersize=2*4.^(magin-min(magin))+45;
%     scatter(lonin,latin,markersize,'Marker','square','MarkerFaceColor',[1 0.4118 0.1608], ...
%         'MarkerEdgeColor','k','LineWidth',0.8);
%     axis xy;
%     axis('tight');
% end

% figure;
% Fun_defaultAxes;
% [lonGrid, latGrid]=meshgrid(linesx, linesy);
% geoshow(latGrid, lonGrid, Gridfore, 'DisplayType', 'mesh');
% colormap(slanCM(colmap));
% c = colorbar();
% c.Label.String = 'Density';
% c.Label.FontWeight = 'bold';
    case 2
%%                              Drawing with map
figure('Color','w','Position',[200, 200, 800, 500]);
ax=geoaxes;
hold(ax,"on");
% You can try other method to get grids mapping data
% Here making matirx become an array
mapw=[];
mapx=Bayesian_ExpandArray(Xcgrid,row);
mapy=Bayesian_RepeatArray(Ycgrid,col);
markersize=2*4.^(magin-min(magin))+45;

for i=1:col
    mapw_copy=Gridfore(:,i);
    mapw_copy=mapw_copy(:)';
    mapw=[mapw, mapw_copy];
end
geodensityplot(ax,mapy,mapx,(mapw),"FaceAlpha","interp","FaceColor","interp");
colormap(ax,slanCM(colmap));
cb=colorbar(ax,'eastoutside','Box','on');
geobasemap(ax,'none');

geoplot(ax,[latMin,latMin,latMax,latMax,latMin],[lonMin,lonMax,lonMax,lonMin,lonMin], ...
    'LineWidth', 2,'Color','r','LineStyle','--');
cb.FontSize=15;
cb.FontWeight='bold';
cb.Label.String='Density';
cb.LineWidth=2;
% ax2=geoaxes;
geoscatter(ax,latin,lonin,'Marker','square','MarkerFaceColor',[1 0.4118 0.1608], ...
        'MarkerEdgeColor','k','LineWidth',0.8);
% Here can not change the range cause it is density
% ticks=linspace(min(mapw(:)), max(mapw(:)), 5);
% tickLabels=arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
% cb.Ticks=(ticks);
% cb.TickLabels=tickLabels;
ax.FontWeight='bold';
ax.LineWidth=2;
ax.FontSize=15;
ax.LatitudeLabel.String='latitude';
ax.LongitudeLabel.String='longitude';
ax.LatitudeLabel.FontSize=15;
ax.LongitudeLabel.FontSize=15;
ax.GridLineStyle="none";
end
end