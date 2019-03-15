%% Main file, here everything is assembled.
%{
Usefull hotkeys:
Ctrl+/ auto-comments the selected lines
Ctrl+T uncomments
%}

%{ 
https://uk.mathworks.com/help/optim/ug/fmincon.html#busog7r-2

fmincon will be out tool for solving this optimization problem
Syntax of fmincon: ------------------------------------------------------
x = fmincon(fun,x0,A,b)
x = fmincon(fun,x0,A,b,Aeq,beq)
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub)
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
x = fmincon(problem)
[x,fval] = fmincon(___)
[x,fval,exitflag,output] = fmincon(___)
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(___)

[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
%}

%% Experimenting with fmincon
% fun = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
% x0 = [-1,2];
% A = [1,2];
% b = 1;
% %constraint is in the form Ax=b
% x = fmincon(fun,x0,A,b) %#ok<*NOPTS>


%% Space 2D:--------------------------------------------------------
hold on; % keep plotting on the existing figure
earth_radius = 300;
axis([-earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200]); %setting up figure size
title("Satelites orbiting the earth and some stations placed on the planet 's surface. -2D ");
step_bound = 50000;%should be the same as the 

sat1 = satelite(30,earth_radius+60,30,10,step_bound); %(arg_vel,arg_alt,arg_init_pos,arg_periods,arg_bound)
sat2 = satelite(60,earth_radius+80,65,10,step_bound);
sat3 = satelite(30,earth_radius+90,60,10,step_bound);
sat4 = satelite(80,earth_radius+50,55,10,step_bound);

station1 = satelite(10,earth_radius,90,10,step_bound);
station2 = satelite(10,earth_radius,10,10,step_bound);

for i=1:step_bound
    earth = viscircles([0 0],earth_radius,'Color','b'); %earth
    
    vis_sat1 = viscircles(sat1.lifetime_coordinates(:,i)',8,'Color','g');
    vis_sat2 = viscircles(sat2.lifetime_coordinates(:,i)',8,'Color','g');
    vis_sat3 = viscircles(sat3.lifetime_coordinates(:,i)',8,'Color','g');
    vis_sat4 = viscircles(sat4.lifetime_coordinates(:,i)',8,'Color','g');

    vis_station1 = viscircles(station1.lifetime_coordinates(:,i)',12,'Color','r');
    vis_station2 = viscircles(station2.lifetime_coordinates(:,i)',12,'Color','r');
    
    pause(0.05) %WARNING: all 'delete' (for graphic objects) functions must be after the 'pause' function 
    delete(vis_sat1) 
    delete(vis_sat2)
    delete(vis_sat3)
    delete(vis_sat4)
    
    delete(vis_station1)
    delete(vis_station2)
    
    delete(earth)
end
hold off;

%% Space 3D:--------------------------------------------------------
hold on; % keep plotting on the existing figure
earth_radius = 200;
axis([-earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200]); %setting up figure size
title("Satelites orbiting the earth and some stations placed on the planet 's surface. -3D ");
rounds = 10;

[x,y,z] = sphere();%get coordinates in order to create 3D spheres. A sphere can be the earth, a satelite or a station.
%%%% ---------[s s   s s s s s s st st]
% velocities = [3 3.2 6 4 3 5 2 3 80 80]*1000; % the smaller the faster (because of the use of linspace() in the satelite3D constructor)
velocities = [3 3 3 3 3 3 3 3 80 80]*10000; %debug
% velocities = [ 1 1 1 1 1 1 1 1 0.1 0.1]*10; %debug
sat1 = satelite3D(20, 20,earth_radius+50, 20, rounds, velocities(1)); % arg_vel is weirdly defined. The bigger it's value the slower the object moves (since it makes more steps to make a full circle)
sat2 = satelite3D(20, 20, earth_radius+50, 30, rounds, velocities(2));
sat3 = satelite3D(20, 20,earth_radius+50, 40, rounds, velocities(3));
sat4 = satelite3D(20, 20, earth_radius+50, 50, rounds, velocities(4));
sat5 = satelite3D(20, 20,earth_radius+50, 60, rounds, velocities(5)); 
sat6 = satelite3D(20, 20, earth_radius+50, 70, rounds, velocities(6));
sat7= satelite3D( 20, 20, earth_radius+50, 80, rounds, velocities(7));
sat8 = satelite3D(20, 20, earth_radius+50, 90, rounds, velocities(8));
sat9 = satelite3D(20, 20, earth_radius+50, 100, rounds, velocities(8));

station1 = satelite3D(20,20,earth_radius,20, rounds, velocities(9));
station2 = satelite3D(20,20,earth_radius,50, rounds, velocities(10));


%%%% WARNING! 't' must be always < 10000(from linspace) of all satelites and stations (avoiding index out of bounds error)
stop = 60;
counter = 1;
times = 10000;
for i = 1:times
   %Coordinates ---------------------------
   vis_sat1_coordinates = sat1.lifetime_coordinates(:,i)';
   vis_sat2_coordinates = sat2.lifetime_coordinates(:,i)';
   vis_sat3_coordinates = sat3.lifetime_coordinates(:,i)';
   vis_sat4_coordinates = sat4.lifetime_coordinates(:,i)';
   vis_sat5_coordinates = sat5.lifetime_coordinates(:,i)';
   vis_sat6_coordinates = sat6.lifetime_coordinates(:,i)';
   vis_sat7_coordinates = sat7.lifetime_coordinates(:,i)';
   vis_sat8_coordinates = sat8.lifetime_coordinates(:,i)';
   vis_sat9_coordinates = sat9.lifetime_coordinates(:,i)';
     
   vis_station1_coordinates = station1.lifetime_coordinates(:,i)';
   vis_station2_coordinates = station2.lifetime_coordinates(:,i)';
   
   %Visualize ---------------------------
   earth = surf( earth_radius*x, earth_radius*y, earth_radius*z );

   vis_sat1 = surf(vis_sat1_coordinates(1)+10*x,vis_sat1_coordinates(2)+10*y,vis_sat1_coordinates(3)+10*z);% create a sphere for the satelite
   vis_sat2 = surf(vis_sat2_coordinates(1)+10*x,vis_sat2_coordinates(2)+10*y,vis_sat2_coordinates(3)+10*z);
   vis_sat3 = surf(vis_sat3_coordinates(1)+10*x,vis_sat3_coordinates(2)+10*y,vis_sat3_coordinates(3)+10*z);
   vis_sat4 = surf(vis_sat4_coordinates(1)+10*x,vis_sat4_coordinates(2)+10*y,vis_sat4_coordinates(3)+10*z);
   vis_sat5 = surf(vis_sat5_coordinates(1)+10*x,vis_sat5_coordinates(2)+10*y,vis_sat5_coordinates(3)+10*z);
   vis_sat6 = surf(vis_sat6_coordinates(1)+10*x,vis_sat6_coordinates(2)+10*y,vis_sat6_coordinates(3)+10*z);
   vis_sat7 = surf(vis_sat7_coordinates(1)+10*x,vis_sat7_coordinates(2)+10*y,vis_sat7_coordinates(3)+10*z);
   vis_sat8 = surf(vis_sat8_coordinates(1)+10*x,vis_sat8_coordinates(2)+10*y,vis_sat8_coordinates(3)+10*z);
   vis_sat9 = surf(vis_sat9_coordinates(1)+10*x,vis_sat9_coordinates(2)+10*y,vis_sat9_coordinates(3)+10*z);
     
   vis_station1 = surf(vis_station1_coordinates(1)+20*x,vis_station1_coordinates(2)+20*y,vis_station1_coordinates(3)+20*z);
   vis_station2 = surf(vis_station2_coordinates(1)+20*x,vis_station2_coordinates(2)+20*y,vis_station2_coordinates(3)+20*z);
   
   if i == stop
       break
   end
   
   %Delete ---------------------------
   pause(0.05) %WARNING: all 'delete' (for graphic objects) functions must be after the 'pause' function 
   delete(vis_sat1)
   delete(vis_sat2)
   delete(vis_sat3)
   delete(vis_sat4)
   delete(vis_sat5)
   delete(vis_sat6)
   delete(vis_sat7)
   delete(vis_sat8)
   delete(vis_sat9)
   
   delete(vis_station1)
   delete(vis_station2)
   
   delete(earth)
end
% hold off %not necessary(?)

%%%% Since i = 
current_coordinates = [vis_sat1_coordinates;
                       vis_sat2_coordinates;
                       vis_sat3_coordinates;
                       vis_sat4_coordinates;
                       vis_sat5_coordinates;
                       vis_sat6_coordinates;
                       vis_sat7_coordinates;
                       vis_sat8_coordinates;
                       vis_sat9_coordinates;
                       vis_station1_coordinates;
                       vis_station2_coordinates;                       
                       ];

%% Creating matrix A:
n = length(current_coordinates(:,1));
DISTANCES = zeros(n,n); %(N+M)x(N+M)
LINKS = zeros(n,n); %(N+M)x(N+M)
communication_range = 200;
for i = 1:n
    for j = 1:n
        DISTANCES(i,j) = euclidean_dist(current_coordinates(i,:),current_coordinates(j,:)); 
        if euclidean_dist(current_coordinates(i,:),current_coordinates(j,:)) <= communication_range %200 is arbitrary
            LINKS(i,j) = 1;
        else
            LINKS(i,j) = 0;
        end
    end
end

%Each row concerns -> a satelite, last 2 rows -> stations
disp('Distance matrix:-----------------------------------------------------')
disp(DISTANCES)
disp('---------------------------------------------------------------------')
disp(LINKS)%should I include the diagonal elements (selfs)?
disp('---------------------------------------------------------------------')

%% Constructing the function
b = [5 ;6 ;7 ;5 ;1 ;2 ;5 ;3 ;4 ;3 ;2]; % capacities
% objective_function = @(s)sum(s);  % UNDER CONSTRUCTION
objective_function = @(x)sum(LINKS*x - LINKS*x);

s = zeros(1,n); %s: divergence, s is a 1x(N+M) vector
% for i = 1:n
%     s(i) = @(x)sum(LINKS(i,:).*x)-sum(LINKS(:,i).*x); %% This is WRONG!
%     IT WONT EVEN RUN (OBVIOUSLY)! If it run a zero sum is expected.
% end 
 



%{
--------------------------------------------------------------------
Function definitions in a script must appear at the end of the file.
Only the first fuction can be accessed by other files. The rest are
acounted as local functions. MATLAB does no support default values on 
function parameters.
--------------------------------------------------------------------
%}
function out = euclidean_dist(vec1, vec2) 
    out = sqrt(sum((vec1 - vec2) .^ 2));
end

