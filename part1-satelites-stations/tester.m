%% Space 3D:--------------------------------------------------------
cd C:\Users\User\Documents\GitHub\Network_Optimization_Thesis\part1-satelites-stations
figure('Name','3D Simulation of satelite orbits');
title("Satelites orbiting the earth and some stations placed on the planet 's surface. -3D ");
hold on; % keep plotting on the existing figure
earth_radius = 200;
axis([-earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200]); %setting up figure size
axis equal

NUMBER_OF_SATELITES = 12;
NUMBER_OF_STATIONS = 2;
DIFFERENT_VELOCITIES = true; % boolean

%%%% ---------[s s   s s s s s s st st] Iverse velocities
nodes = create_nodes(NUMBER_OF_SATELITES,NUMBER_OF_STATIONS,3,80,DIFFERENT_VELOCITIES);
%%%% WARNING! 't' must be always < min(velocities)(lesser from linspace) of all satelites and stations (avoiding index out of bounds error)
rounds = 10;
stop = 30;
times = 10000;
[x,y,z] = sphere();%get coordinates in order to create 3D spheres. A sphere can be the earth, a satelite or a station.
for i = 1:times
    
   %Coordinates ---------------------------
   coords = [];
   for j = 1:length(nodes)
       coords = [coords; nodes(j).lifetime_coordinates(:,i)' ];
   end
   
   %Visualize ---------------------------
   earth = surf( earth_radius*x, earth_radius*y, earth_radius*z );
   displays = [];
   for j = 1:length(nodes)
      if j <= NUMBER_OF_SATELITES
        displays = [displays, surf(coords(j,1)+10*x,coords(j,2)+10*y,coords(j,3)+10*z)]; 
      else
        displays = [displays, surf(coords(j,1)+20*x,coords(j,2)+20*y,coords(j,3)+20*z)]; 
      end
   end
   
   if i == stop
       break
   end
   
   %Delete ---------------------------
   pause(0.01) %WARNING: all 'delete' (for graphic objects) functions must be after the 'pause' function 
   
   for j = 1:length(displays)
      delete(displays(j))
   end
   delete(earth)
   disp(i)
end
% hold off %not necessary(?)


function out = create_nodes(num_satelites, num_stations, sat_inverse_vel, stat_inverse_vel, random_factor)
% Creates a list of nodes in the form of [satelite_objects .... station_objects]
%
% - num_satelites: number of satelites (set default as 9)
% - num_stations: number of stations (set default as 2)
% - sat_inverse_vel: inverse velocity of satelites (smaller values imply high
% velocity) use values around 3
% - stat_inverse_vel: inverse velocity of satelites (every station has the
% same value) (smaller values imply high velocity) use values around 80
% - random_factor: boolean, if true it adds some random noise from uniform distribution in [-10,10] on the satelite velocities

    earth_radius = 200;
    rounds = 10;
    matr = [ ];
    inverse_velocities = [ones(1,num_satelites)*sat_inverse_vel, ones(1,num_stations)*stat_inverse_vel]*1000; % the bigger the value the slower the velocity
    
    if random_factor
        rng(10) %setting seed
        rands = [(rand(1,num_satelites)*40)-20, zeros(1,num_stations)];
        %(rand(1,num_satelites)*20)-10 -> rand() returns values in (0,1) so
        %I do this in order to turn it into (-10,10)
        inverse_velocities = inverse_velocities + rands;
    end

    % Remember: satelite3D(arg_theta,arg_phi, arg_alt, arg_init_pos, arg_periods, arg_vel, arg_name)
    % Constructing satelites & stations:
    initializer_difs = (0:(num_satelites+num_stations-1))*40; %used to initialiaze satelites at different positions in order to not overlap each other
    for ii = 1:(num_satelites + num_stations)
        if ii <= num_satelites
            matr = [matr, satelite3D(20, 20, earth_radius+50, 20 + initializer_difs(ii), rounds, inverse_velocities(ii), 'satelite')];
        else
            matr = [matr, satelite3D(20, 20, earth_radius, 200 + initializer_difs(ii)*4, rounds, (num_satelites+num_stations)*1000, 'station')];
        end
    end
    
    
    disp('[~Report:] Satelites and stations successfully created')
    out = matr;
end
