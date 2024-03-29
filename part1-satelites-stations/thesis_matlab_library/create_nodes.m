function out = create_nodes(num_satelites, num_stations, sat_inverse_vel, stat_inverse_vel, random_factor, theta_phi_ca, init_pos)
% Creates a list of node objects in the form of [satelite_objects .... station_objects]
%
% - num_satelites: number of satelites (set default as 9)
% - num_stations: number of stations (set default as 2)
% - sat_inverse_vel: vector, inverse velocity of satelites (smaller values imply high
% velocity) use values around [3 3 3 ... 3]
% - stat_inverse_vel: vector, inverse velocity of satelites (every station has the
% same value) (smaller values imply high velocity) use values around [80 80 80 ... 80]
% - random_factor: boolean, if true it adds some random noise from uniform distribution in [-10,10] on the satelite velocities
% - theta_phi: arg_theta and arg_phi values.

    earth_radius = 200;
    periods = 1;
    matr = [ ];
    inverse_velocities = [sat_inverse_vel, stat_inverse_vel]; % the bigger the value the slower the velocity
    
    if random_factor
        rng(79) %setting seed
        random_numbers_sats = (rand(1,num_satelites)*8+1)-(4-1);% (rand(1,num_satelites)*number*2)-number -> rand() returns values in (0,1) so
                                                                % I do this in order to turn it (the output range) into (-number,number)
        rands = [random_numbers_sats, zeros(1,num_stations)];% adding random factor only at satelite velocities
        
        inverse_velocities = inverse_velocities + rands;
    end
%     inverse_velocities = inverse_velocities * 1000;
       
    % Constructing satelites & stations:
    rand_int = 8;
    initializer_difs = (0:(num_satelites-1))*360*periods* rand_int ./ inverse_velocities(1:num_satelites); % 27 %used to initialiaze satelites at different positions in order to not overlap each other
    %     initializer_difs = init_pos; % UNCOMMENT THIS TO RUN THE SIMPLE COMPARISON EXAMPLE 
    for i = 1:(num_satelites + num_stations)
        if i <= num_satelites 
            matr = [matr, satelite3D(theta_phi_ca{i}(1), theta_phi_ca{i}(2), earth_radius+50, 20 + initializer_difs(i), periods, inverse_velocities(i), 'satellite')];
        else
            matr = [matr, satelite3D(theta_phi_ca{i}(1), theta_phi_ca{i}(2), earth_radius, initializer_difs(i-num_satelites)*2 + i, periods, inverse_velocities(i), 'station')];
        end
    end
    
%     disp('[~Report:] Satellites and stations successfully created')
    out = matr;
end
