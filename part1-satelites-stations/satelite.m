classdef satelite < handle %why do i have to use as a base class the handle one???
    %initialize as: satelite(arg_vel,arg_alt,arg_init_pos,arg_periods,arg_bound)
    %Contains all the information about a satelite-node.
    %Class describing a satelite object. Every satelite is also a source in graph
    %terms.
    
    % PROPERTIES
    properties(SetAccess = private) %Only satelite's methods must be able to access these
       initial_position = 1
%        longitude = 0 % current longitude of satelite (initialised as cartesian, after first call of net_pos it is transdormed to polar)
%        latitude = 0 % current latitude of satelite (initialiased as cartesian, after first call of net_pos it is transdormed to polar)
       angular_velocity = 0 % measured in degrees/sec. Float number, it can take both positive and negative values. https://en.wikipedia.org/wiki/Orbital_speed
       satelite_id = 0  % CHANGE THIS. It must be an increasing number from 1 to N
       altitude = 0; % 408.000 meters (like the ISS)
    end
   
   properties %public properties
       neighbours = [] %binary list (of size N+M) containing the other sources (or sinks) that have a connection with the satelite.
       lifetime_coordinates = [];  %a matrix containing the coordinates of the satelite at any given time during it's lifetime
   end
   % METHODS
   methods
       
       function sat_obj = satelite(arg_vel, arg_alt, arg_init_pos, arg_periods,arg_bound) %constructor, arg_init_pos=initial position(degrees), arg_periods=(e.x. 2*360)
            sat_obj.angular_velocity = arg_vel; %set angular velocity
            sat_obj.altitude = arg_alt;%set altitude
            sat_obj.initial_position = arg_init_pos;%set initial position
            steps = linspace(sat_obj.initial_position, sat_obj.initial_position + arg_periods*360, arg_bound); %(initial angle, initial angle+-n circles, step_bound)
            %step creates a vector containing the range of a step done
            %every time the satelite moves. 
            sat_obj.lifetime_coordinates = [(arg_alt)*cosd(steps*arg_vel); (arg_alt)*sind(steps*arg_vel)]; %calculating life_time coordinates given the 'steps' vactor
       end
               
%        function next_pos(obj,step) %calculate and update the satelite's polar coordinates given it's angular velocity
%            radius = obj.altitude;
%            obj.longitude = radius*cosd(obj.initial_position*step*obj.angular_velocity);%Long = X (polar)
%            obj.latitude = radius*sind(obj.initial_position*step*obj.angular_velocity); %Lat = Y (polar)
%        end
       
       
   end
end

