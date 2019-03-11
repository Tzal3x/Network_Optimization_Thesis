% --------- HELPING SCRIPT FOR EXPERIMENTATION AND DEBUGING. IGNORE IT -----------------------
step_bound = 1000;
hold on % keep plotting on the existing figure
title('Satelites orbiting the earth and some stations fixed in place')
rad = 408;

axis([-rad-100 rad+100 -rad-100 rad+100]) %setting up figure size
% viscircles([0 0],408-88) %earth 2D

%Earth 3D
[x,y,z] = sphere();
r = 408;
surf( r*x, r*y, r*z ) 

% creating satelites:
sat1 = satelite(110,100,50);
sat2 = satelite(550,200,-60);
sat3 = satelite(150,300,-70);
sat4 = satelite(740,500,80);
sat5 = satelite(110,400,60);

sats = [sat1,sat2,sat3,sat4,sat5]; %list containing every satelite

for step=1:0.1:step_bound
    %Move every satelite:
    for sati = sats
        sati.next_pos(step);
%         plot(sati.get_lat, sati.get_long,'r*'); %2D
    plot3(sati.get_lat, sati.get_long,408,'r*');
    end
    drawnow
    pause(0.3) % 1 frame per second for animation purposes
end
hold off

%% test stations
stat1 = station(2,earth_radius,10)
stat2 = station(2,earth_radius,16)
