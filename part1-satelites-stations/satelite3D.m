classdef satelite3D < handle
    %satelite3D(arg_theta,arg_phi, arg_alt, arg_init_pos, arg_periods, arg_vel)
    %3D version of the satellite object

   properties %public properties
%        neighbours = [] %binary list (of size N+M) containing the other sources (or sinks) that have a connection with the satelite.
       lifetime_coordinates = [];  %a matrix containing the coordinates of the satelite at any given time during it's lifetime
   end
   
   methods
       function sat_obj = satelite3D(arg_theta,arg_phi, arg_alt, arg_init_pos, arg_periods, arg_vel) %constructor, arg_init_pos=initial position(degrees), arg_periods=(e.x. 2*360)
            steps = linspace(arg_init_pos, arg_init_pos + arg_periods*360, arg_vel);
%             sat_obj.lifetime_coordinates = [(arg_alt)*sind(steps*arg_theta).*cosd(ones(1,length(steps))*arg_phi);
%                                             (arg_alt)*sind(ones(1,length(steps))*arg_theta).*sind(steps*arg_phi);
%                                             (arg_alt)*cosd(steps.*arg_theta) ];
                                        
            sat_obj.lifetime_coordinates = [(arg_alt)*sind(steps*arg_theta).*cosd(ones(1,length(steps))*arg_phi);
                                            (arg_alt)*sind(ones(1,length(steps))*arg_theta).*sind(steps*arg_phi);
                                            (arg_alt)*cosd(steps.*arg_theta) ];
       
       end       
   end 
end
    
    