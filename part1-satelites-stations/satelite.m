classdef satelite < handle %why do i have to use as a base class the handle one???
    %initialize as: satelite(arg_vel,arg_alt,arg_steps)
    %Contains all the information about a satelite-node.
    %Class describing a satelite object. Every satelite is also a source in graph
    %terms.
    
    % PROPERTIES
    properties(SetAccess = private) %Only satelite's methods must be able to access these
       longitude = 0 % current longitude of satelite (initialised as cartesian, after first call of net_pos it is transdormed to polar)
       latitude = 0 % current latitude of satelite (initialiased as cartesian, after first call of net_pos it is transdormed to polar)
       angular_velocity = 0 % measured in degrees/sec. Float number, it can take both positive and negative values. https://en.wikipedia.org/wiki/Orbital_speed
       satelite_id = 0  % CHANGE THIS. It must be an increasing number from 1 to N
       altitude = 0; % 408.000 meters (like the ISS)
    end
   
   properties %public properties
       neighbours = [] %binary list (of size N+M) containing the other sources (or sinks) that have a connection with the satelite.
   end
   % METHODS
   methods
       % --- Getters
       function out = get_lat(obj)%returns current latitude
           out = obj.latitude;            
       end
       
       function out = get_long(obj)%returns current longitude
           out = obj.longitude;
       end
       
       function sat_obj = satelite(arg_vel,arg_alt,arg_steps) %constructor, input arguments: longitude, latitude, angular velocity
            sat_obj.angular_velocity = arg_vel;
            sat_obj.altitude = arg_alt;
            sat_obj.longitude = arg_alt*cosd(arg_vel*arg_steps);
            sat_obj.latitude = arg_alt*sind(arg_vel*arg_steps);
       end
       
       function next_pos(obj, steps) %calculate and update the satelite's polar coordinates given it's angular velocity
           radius = obj.altitude;
           obj.longitude = radius*cosd(steps*obj.angular_velocity);%Long = X (polar)
           obj.latitude = radius*sind(steps*obj.angular_velocity); %Lat = Y (polar)
       end
   end
end

