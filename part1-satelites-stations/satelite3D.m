classdef satelite3D < handle
    %satelite3D(arg_theta,arg_phi, arg_alt, arg_init_pos, arg_periods, arg_vel, arg_name)
    %3D version of the satellite object. arg_theta and arg_phi should be
    %equal.

   properties %public properties
       lifetime_coordinates = [];  %a matrix containing the coordinates of the satelite at any given time during it's lifetime
       name;
   end
   
   methods
       function sat_obj = satelite3D(arg_theta,arg_phi, arg_alt, arg_init_pos, arg_periods, arg_vel, arg_name) %constructor, arg_init_pos=initial position(degrees), arg_periods=(e.x. 2*360)
            if arg_vel>=0
                steps = linspace(arg_init_pos, arg_init_pos + arg_periods*360, arg_vel);
            else 
                if arg_vel<0
                    steps = -linspace(arg_init_pos, arg_init_pos + arg_periods*360, -arg_vel);
                end
            end
            sat_obj.name = arg_name;
            
            sat_obj.lifetime_coordinates = [(arg_alt)*sind(steps+arg_theta).*cosd(ones(1,length(steps))*arg_phi);
                                            (arg_alt)*sind(ones(1,length(steps))*arg_theta).*sind(steps+arg_phi);
                                            (arg_alt)*cosd(steps+arg_theta) ];
% EXPERIMENTAL:                                          
%             sat_obj.lifetime_coordinates = [(arg_alt)*sind(steps+arg_theta).*cosd(ones(1,length(steps))*arg_phi);
%                                             (arg_alt)*sind(ones(1,length(steps))*arg_theta).*sind(steps+arg_phi);
%                                             (arg_alt)*cosd(steps+arg_theta) ];
                                          
       end       
   end 
end
    
    