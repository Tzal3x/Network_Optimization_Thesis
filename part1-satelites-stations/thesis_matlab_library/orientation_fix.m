% akala

r = 10;
thetainit = 45;
phiinit = 30;
% 
R1 = [cosd(thetainit) -sind(thetainit) 0;...
      sind(thetainit) cosd(thetainit)  0;...
      0               0                1];
R2 = [1 0 0;...
      0 cosd(phiinit) -sind(phiinit);...
      0 sind(phiinit) cosd(phiinit)];
R = R1 * R2;
t = 0:1:360; % lifetime coordinates
% Metatropi se sferikes sydetagmenes
x2 = r * cosd(t); 
y2 = r * sind(t);
z2 = zeros(1,length(t));
% Convert to degrees to rads
x2 = deg2rad(x2);
y2 = deg2rad(y2);
z2 = deg2rad(z2);
% Apply transform
xyz = R * [x2; y2; z2];
% Plot 3D
scatter3(xyz(1,:), xyz(2,:), xyz(3,:))
