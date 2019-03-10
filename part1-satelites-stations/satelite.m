classdef satelite < handle %why do i have to use as a base class the handle one???
    %{
    Contains all the information about a satelite-node.
    Each satelite has a constant communication range, therefore the satelite can be
    expressed as the center of a sphere and the communication as the radius of the sphere.
    %}
    %% PROPERTIES
    properties(SetAccess = private) %Only satelite's methods must be able to access these
       init_long = 0 %(degrees), initial longitude, using this in order to have a point of reference
       init_lat = 0 %(degrees), initial latitude, using this in order to have a point of reference
       longitude = 0 % current longitude of satelite
       latitude = 0 % current latitude of satelite
       satelite_id = 0  % change this. It must be an increasing number from 1 to N
       angular_velocity = 0 % measured in degrees/sec. Float number, it can take both positive and negative values. https://en.wikipedia.org/wiki/Orbital_speed
    end
   properties(Constant)
       altitude = 408; % 408.000 meters (like the ISS)
   end
   
   properties %public properties
       neighbours = [] %binary list (of size N+M) containing the other sources (or sinks) that have a connection with the satelite.
   end
   %% METHODS
   methods
       % --- Getters
       function out = get_lat(obj)%returns current latitude
           out = obj.latitude;            
       end
       
       function out = get_long(obj)%returns current longitude
           out = obj.longitude;
       end
       
       function sat_obj = satelite(arg_long,arg_lat,arg_vel) %constructor, input arguments: longitude, latitude, angular velocity
          if nargin == 3 %if number of variables (...). If this statement is skipped, initialising the constructor with less than 3 errors will cause an error
            sat_obj.init_long = arg_long;
            sat_obj.init_lat = arg_lat;
            sat_obj.angular_velocity = arg_vel;
          end
       end
       
       function next_pos(obj, step) %calculate and update the satelite's polar coordinates given it's angular velocity
           %step: next step (in degrees) at time t+1
           radius = obj.altitude;
           obj.longitude = radius*cosd(obj.init_long + step*obj.angular_velocity);%Long = X
           obj.latitude = radius*sind(obj.init_lat + step*obj.angular_velocity); %Lat = Y 
       end
   end
end

