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
           sat_obj.name = arg_name;
 
           if arg_vel>=0
                steps = linspace(arg_init_pos, arg_init_pos + arg_periods*360, arg_vel);
           else 
                if arg_vel<0
                    steps = -linspace(arg_init_pos, arg_init_pos + arg_periods*360, -arg_vel);
                end
            end
            % Calculating orientation matrixes
            R1 = [cos(arg_theta) -sin(arg_theta) 0;...
                  sin(arg_theta) cos(arg_theta)  0;...
                  0               0                1];
            R2 = [1 0 0;...
                  0 cos(arg_phi) -sin(arg_phi);...
                  0 sin(arg_phi) cos(arg_phi)];
            R = R1 * R2;
            steps = deg2rad(steps);

            x = arg_alt * cos(steps); 
            y = arg_alt * sin(steps);
            z = zeros(1,length(steps));
            sat_obj.lifetime_coordinates = R * [x;y;z]; 
            % lifetime_coordinates: 3xn where 1st row is x, 2nd is y and
            % 3rd is z.
       end       
   end 
end
    
    