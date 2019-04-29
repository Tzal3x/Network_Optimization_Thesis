% akala

r = 10;
thetainit = 45;
phiinit = 30;

t = 0:1:360;
R1 = [cosd(thetainit) -sind(thetainit) 0;...
      sind(thetainit) cosd(thetainit)  0;...
      0               0                1];
  
R2 = [1 0 0;...
      0 cosd(phiinit) -sind(phiinit);...
      0 sind(phiinit) cosd(phiinit)];
%R1 kai R2 einai sxedon idia metatopismena 
R = R1 * R2;

x2 = r * cosd(t);
y2 = r * sind(t);
z2 = zeros(1,length(t));
xyz = R * [x2; y2; z2];
scatter3(xyz(1,:), xyz(2,:), xyz(3,:))
