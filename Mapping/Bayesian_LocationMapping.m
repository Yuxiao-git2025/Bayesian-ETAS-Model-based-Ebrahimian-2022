%     'satellite';     % 卫星影像（彩色）-1
%     'colorterrain';  % 彩色地形图-2
%     'bluegreen';     % 蓝绿色海洋和绿色陆地-3
%     'landcover';     % 土地覆盖图（彩色）-4
%     'streets';       % 街道地图（彩色）-5
%     'streets-light'; % 浅色街道地图（default）-6
%     'streets-dark';  % 深色街道地图-7
%     'topographic';   % 地形图（彩色）-8
%     'grayterrain';   % 灰色地形图-9
function ax=Bayesian_LocationMapping(latbound,lonbound,colmap,idxmap)
global latitude longitude M
if isempty(latbound) && isempty(lonbound)
    latbound=[min(latitude) max(latitude)];
    lonbound=[min(longitude) max(longitude)];
end
if isempty(colmap)
    colmap='jet';
end
if isempty(idxmap)
    idxmap=2 ;
end
if idxmap>9
    error('Please try smaller idx');
end
f=figure("Color",'w');
ax=geoaxes;
geolimits(ax,latbound,lonbound);
basemaps = {
    'satellite';     % 卫星影像（彩色）-1
    'colorterrain';  % 彩色地形图-2
    'bluegreen';     % 蓝绿色海洋和绿色陆地-3
    'landcover';     % 土地覆盖图（彩色）-4
    'streets';       % 街道地图（彩色）-5
    'streets-light'; % 浅色街道地图-6
    'streets-dark';  % 深色街道地图-7
    'topographic';    % 地形图（彩色）-8
    'grayterrain';   %灰色地形图-9
};

geobasemap(ax,basemaps{idxmap}); 
hold(ax,'on');

ax.FontWeight='bold';
ax.LineWidth=2;
ax.FontSize=15;
ax.LatitudeLabel.String='latitude';
ax.LongitudeLabel.String='longitude';
ax.LatitudeLabel.FontSize=15;
ax.LongitudeLabel.FontSize=15;
ax.GridLineStyle="none";

% Main plot
markerSize=min(1*4.^(M-min(M))+10,60);
% Colors=Fun_Value2Color(M,length(M),colmap,1);
geoscatter(ax, latitude, longitude, markerSize,M,'o','filled', ...
    'MarkerEdgeColor',[0 0 0],'LineWidth',0.6);
% geobubble(ax,latitude,longitude,markerSize);
c =colorbar(ax);
c.LineWidth=1;
c.Label.FontWeight='bold';
c.Label.String=' ';
colormap(ax,slanCM(colmap));
id=find(M==max(M));
geoscatter(ax,latitude(id),longitude(id),350,'Marker','pentagram','MarkerFaceColor','r', ...
    'MarkerEdgeColor','k','LineWidth',0.6)
% geodensityplot(ax,latitude,longitude,M,'FaceColor',[ 1.0000    0.0745    0.6510]);


end