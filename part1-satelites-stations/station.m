classdef station < handle
    %initialized as: station(arg_vel, arg_alt, arg_steps)
    %Contains all the information about a station-node.
    %Class describing a station object. Every station is a sink in graph
    %terms.
    properties
        longitude = 0; 
        latitude = 0; 
        altitude = 0;
        angular_velocity = 0; % if earth is rotating it should be 0, else it should be the same number in every station
        neighbours = [] %binary list (of size N+M) containing the other sources (or sinks) that have a connection with the station
    end
    
    methods
        function station_obj = station(arg_vel, arg_alt, arg_steps) 
            station_obj.angular_velocity = arg_vel;
            station_obj.altitude = arg_alt;
            station_obj.longitude = arg_alt*cosd(arg_vel*arg_steps);
            station_obj.latitude = arg_alt*sind(arg_vel*arg_steps);
        end
        
        function next_pos(obj, steps)
            %warning: if angular_velocity = 0, the position will remain
            %the same
           radius = obj.altitude;
           obj.longitude = radius*cosd(steps*obj.angular_velocity);%Long = X (degrees)
           obj.latitude = radius*sind(steps*obj.angular_velocity); %Lat = Y  (degrees)
        end
       % Getters -------------------------------------------------------- 
       function out = get_lat(obj)%returns current latitude
           out = obj.latitude;            
       end
       
       function out = get_long(obj)%returns current longitude
           out = obj.longitude;
       end
    end %end methods
end