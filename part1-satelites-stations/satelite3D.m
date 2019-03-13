classdef satelite3D < handle
    %satelite(arg_vel,arg_vel2, arg_alt, arg_init_pos, arg_periods,arg_bound)
    %3D version of the satellite object
   properties(SetAccess = private) %Only satelite's methods must be able to access these
       initial_position = 1
       angular_velocity = 0 %(theta) measured in degrees/sec. Float number, it can take both positive and negative values. https://en.wikipedia.org/wiki/Orbital_speed
       angular_velocity2 = 0 %(phi)
       satelite_id = 0  % CHANGE THIS. It must be an increasing number from 1 to N
       altitude = 0; %  = radius of sphere
    end
   
   properties %public properties
       neighbours = [] %binary list (of size N+M) containing the other sources (or sinks) that have a connection with the satelite.
       lifetime_coordinates = [];  %a matrix containing the coordinates of the satelite at any given time during it's lifetime
   end
   % METHODS------------------------
   methods
       function sat_obj = satelite3D(arg_vel,arg_vel2, arg_alt, arg_init_pos, arg_periods,arg_bound) %constructor, arg_init_pos=initial position(degrees), arg_periods=(e.x. 2*360)
            sat_obj.angular_velocity = arg_vel; %set angular velocity
            sat_obj.altitude = arg_alt;%set altitude
            sat_obj.initial_position = arg_init_pos;%set initial position
            steps = linspace(sat_obj.initial_position, sat_obj.initial_position + arg_periods*360, arg_bound);
%             sat_obj.lifetime_coordinates = [(arg_alt)*sind(steps*arg_vel).*cosd(steps*arg_vel2); (arg_alt)*sind(steps*arg_vel).*sind(steps*arg_vel2);(arg_alt)*cosd(steps.*arg_vel) ]; %[x;y;z]
            sat_obj.lifetime_coordinates = [(arg_alt)*sind(steps*arg_vel).*cosd(ones(1,length(steps))*arg_vel2); (arg_alt)*sind(ones(1,length(steps))*arg_vel).*sind(steps*arg_vel2);(arg_alt)*cosd(steps.*arg_vel) ]; %[x;y;z]
       end       
   end 
end
    
    