classdef satelite < handle %why do i have to use as a base class the handle one???
    %{
    Contains all the information about a satelite-node.
    Each satelite has a constant communication range, therefore the satelite can be
    expressed as the center of a sphere and the communication as the radius of the sphere.
    If a sphere is 
    %}
    properties(SetAccess = private) %Only satelite's methods must be able to access these
       init_long = 0 % initial longitude, using this in order to have a point of reference
       init_lat = 0 % initial latitude, using this in order to have a point of reference
       longitude = 0 % current longitude of satelite
       latitude = 0 % current latitude of satelite
       satelite_id = 0  % change this. It must be an increasing number from 1 to N
       angular_velocity = 0 % measured in degrees/sec. Float number, it can take both positive and negative values. https://en.wikipedia.org/wiki/Orbital_speed
    end
   properties(Constant)
       altitude = 408; % 408.000 meters (like the ISS)
   end
   
   properties
       neighbours = [] %binary list (of size N+M) containing the other sources (or sinks) that have a connection with the satelite.
   end
   
   methods
       % --- Getters
       function out = get_lat(obj)
           out = obj.latitude;            
       end
       function out = get_long(obj)
           out = obj.longitude;
       end
       
       function sat_obj = satelite(arg_long,arg_lat,arg_vel) %constructor
          if nargin == 3 %if number of variables (...). If this statement is skipped, initialising the constructor with less than 3 errors will cause an error
            sat_obj.init_long = arg_long;
            sat_obj.init_lat = arg_lat;
            sat_obj.angular_velocity = arg_vel;
          end
       end
       
       function next_pos(obj, step)
           %step: next step (in degrees) at time t+1
           
           radius = obj.altitude;
           obj.longitude = radius*cosd(obj.init_long + step);
           obj.latitude = radius*sind(obj.init_lat + step);
       end
       
%        function out = is_neighbours(obj) % checks if a satelite has a
%        link with another entity
%            x = obj.latitude*obj.longitude;
%            out = x;
%        end
   end
end

