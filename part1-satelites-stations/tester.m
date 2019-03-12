% --------- HELPING SCRIPT FOR EXPERIMENTATION AND DEBUGING. IGNORE IT -----------------------
%% test satelite
axis([-earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200]); %setting up figure size
earth_radius = 300;
step_bound = 10000;
sat1 = satelite(30,earth_radius+60,30,10);
sat2 = satelite(60,earth_radius+100,70,10);

hold on; % keep plotting on the existing figure
title("Satelites orbiting the earth and some stations placed on the planet 's surface.");

for i=1:step_bound
    earth = viscircles([0 0],earth_radius,'Color','b'); %earth
    
    vis_sat1 = viscircles(sat1.lifetime_coordinates(:,i)',8,'Color','g');
    vis_sat2 = viscircles(sat2.lifetime_coordinates(:,i)',8,'Color','g');

    
    pause(0.01) %WARNING: all 'delete' (for graphic objects) functions must be after the 'pause' function 
    delete(vis_sat1) 
    delete(vis_sat2) 
    delete(earth)
end
hold off;



%% Space set up:
% Constants
step_bound = 1000;
earth_radius = 300;

sat1_angles = linspace(30*8, 30*8 + 2*360, step_bound); %(arxikh gwnia, arxikh gwnia+-2kykloi, bhma)
sat2_angles = linspace(-26*19, -26*19 - 2*360, step_bound);
stat1_angles = linspace(2*10, 2*10 + 2*360, step_bound);
stat2_angles = linspace(2*81, 2*81 + 2*360, step_bound);

scoordinates1 = [(earth_radius+60)*cosd(sat1_angles); (earth_radius+60)*sind(sat1_angles)]; %2x1000 pinakes
scoordinates2 = [(earth_radius+30)*cosd(sat2_angles); (earth_radius+30)*sind(sat2_angles)];

tcoordinates1 = [earth_radius*cosd(stat1_angles); earth_radius*sind(stat1_angles)];
tcoordinates2 = [earth_radius*cosd(stat2_angles); earth_radius*sind(stat2_angles)];

%% For Plotting
hold on; % keep plotting on the existing figure
title("Satelites orbiting the earth and some stations placed on the planet 's surface.");
axis([-earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200]); %setting up figure size
for i=1:step_bound
    earth = viscircles([0 0],earth_radius,'Color','b'); %earth
   
    vis_sat1 = viscircles(scoordinates1(:,i)',8,'Color','g'); %8elw iosth sthlh, kai ton ' gia anastrofo
    vis_sat2 = viscircles(scoordinates2(:,i)',8,'Color','g');
    vis_stat1 = viscircles(tcoordinates1(:,i)',12,'Color','r');
    vis_stat2 = viscircles(tcoordinates2(:,i)',12,'Color','r');
       
    pause(0.02) %WARNING: all 'delete' (for graphic objects) functions must be after the 'pause' function 
    delete(vis_sat1) 
    delete(vis_sat2) 
    delete(vis_stat1) 
    delete(vis_stat2)
    delete(earth)
end
hold off;





%% Earth 3D
[x,y,z] = sphere();
r = 408;
surf( r*x, r*y, r*z ) 