%% Space 3D:--------------------------------------------------------
%{
 Note that this is not a very robust implementation of the creation of the
 network.
- - - - - - - - - - - - - - - - -
 My Debug tools: 
1) disp('[~DEBUG:]'+string(object_of_interest))% debug
2) pause(1000)% debug 
%}

cd C:\Users\User\Documents\GitHub\Network_Optimization_Thesis\part1-satelites-stations

% Main parameters: - - - - -- - - - -- - - - - - - - - - - 
% I should consider adding THETA, PHI and ALTITUDE vectors on the
% "create_nodes" function.
NUMBER_OF_SATELITES = 9; % integer, default 9
NUMBER_OF_STATIONS = 2; % integer, default 2
RANDOM_VELOCITIES = false; % boolean, default false
INVERSE_VELOCITIES_SATEL = ones(1,NUMBER_OF_SATELITES) * 3; % smaller value -> faster
INVERSE_VELOCITIES_STATIONS = ones(1,NUMBER_OF_STATIONS) * 80; % larger value -> slower
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

nodes = create_nodes(NUMBER_OF_SATELITES, NUMBER_OF_STATIONS, INVERSE_VELOCITIES_SATEL, ... 
                     INVERSE_VELOCITIES_STATIONS, RANDOM_VELOCITIES); 

% Setting up figure display options: - - - - - - - - - - - - - - - - -
figure('Name','3D Simulation of satelite orbits');
title("satelites:"+string(NUMBER_OF_SATELITES)+", stations:"+string(NUMBER_OF_STATIONS)+", random velocities:"+string(RANDOM_VELOCITIES));
hold on; % keep plotting on the existing figure
earth_radius = 200;
axis([-earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200]); %setting up figure size
axis equal

%%%% WARNING! 't' must be always < min(velocities)(lesser from linspace) of all satelites and stations (avoiding index out of bounds error)
rounds = 10;
stop = 10; % 30
times = 10000;
[x,y,z] = sphere();%get coordinates in order to create 3D spheres. A sphere can be the earth, a satelite or a station.
for i = 1:times

   %Coordinates ---------------------------
   coords = [];
   for j = 1:length(nodes)
       disp('start----'+string(j))
       disp(nodes(j).lifetime_coordinates(:,i)')
       disp('middle----'+string(j))
       coords = [coords; nodes(j).lifetime_coordinates(:,i)' ];
       disp('end----'+string(j))
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
   
   % Placing here the 'stop' break for the figure to stay active with all
   % objects:
   if i == stop
       break
   end
   
   %Delete (to create the animation illusion)
   pause(0.01) %WARNING: all 'delete' (for graphic objects) functions must be after the 'pause' function 
   for j = 1:length(displays)
      delete(displays(j))
   end
   delete(earth)
   disp(i) % debug
end
% hold off %not necessary(?)


function out = create_nodes(num_satelites, num_stations, sat_inverse_vel, stat_inverse_vel, random_factor)
% Creates a list of nodes in the form of [satelite_objects .... station_objects]
%
% - num_satelites: number of satelites (set default as 9)
% - num_stations: number of stations (set default as 2)
% - sat_inverse_vel: vector, inverse velocity of satelites (smaller values imply high
% velocity) use values around [3 3 3 ... 3]
% - stat_inverse_vel: vector, inverse velocity of satelites (every station has the
% same value) (smaller values imply high velocity) use values around [80 80 80 ... 80]
% - random_factor: boolean, if true it adds some random noise from uniform distribution in [-10,10] on the satelite velocities

    earth_radius = 200;
    rounds = 10;
    matr = [ ];
    inverse_velocities = [sat_inverse_vel, stat_inverse_vel]; % the bigger the value the slower the velocity
    
    if random_factor
        rng(79) %setting seed
        random_numbers_sats = (rand(1,num_satelites)*8+1)-(4-1);%(rand(1,num_satelites)*number*2)-number -> rand() returns values in (0,1) so
                                                                % I do this in order to turn it into (-number,number)
        rands = [random_numbers_sats, zeros(1,num_stations)];
        
        inverse_velocities = inverse_velocities + rands;
    end
    inverse_velocities = inverse_velocities * 1000;
    
    % Remember: satelite3D(arg_theta,arg_phi, arg_alt, arg_init_pos, arg_periods, arg_vel, arg_name)
    % Constructing satelites & stations:
    initializer_difs = (0:(num_satelites+num_stations-1))*27; %used to initialiaze satelites at different positions in order to not overlap each other
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
