function x=topgeo(c1,c2,lar,lor,opz)
% Transform [lat,lon] Coordinates to Map Coordinate (Cartesian coordinate system)
Rt=6371;
if nargin < 5, opz=0; end
if opz == 0
   x=[Rt*tan(c2-lor).*cos(c1)./cos(c1-lar),Rt*tan(c1-lar)];
else
   x1=lar+atan(c2/Rt);
   y1=lor+atan(c1.*cos(x1-lar)./cos(x1)/Rt);
   x=[x1,y1];
end

return