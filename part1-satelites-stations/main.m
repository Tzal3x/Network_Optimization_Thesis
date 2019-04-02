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
cd C:\Users\User\Documents\GitHub\Network_Optimization_Thesis\part1-satelites-stations
hold on; % keep plotting on the existing figure
earth_radius = 200;
axis([-earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200]); %setting up figure size
axis equal

title("Satelites orbiting the earth and some stations placed on the planet 's surface. -3D ");
rounds = 10;

[x,y,z] = sphere();%get coordinates in order to create 3D spheres. A sphere can be the earth, a satelite or a station.
%%%% ---------[s s   s s s s s s st st]
% velocities = [3 3.2 6 4 3 5 2 3 80 80]*1000; % the smaller the faster (because of the use of linspace() in the satelite3D constructor)
% velocities = [ 1 1 1 1 1 1 1 1 0.1 0.1]*10; %debug
velocities = [3 3 3 3 3 3 3 3 80 80]*1000; %easy mode
% velocities = [3 3.2 6 4 3 5 2 3 80 80]*1000; %hard mode

sat1 = satelite3D(20, 20, earth_radius+50, 20, rounds, velocities(1),'satelite'); % arg_vel is weirdly defined. The bigger it's value the slower the object moves (since it makes more steps to make a full circle)
sat2 = satelite3D(20, 20, earth_radius+50, 60, rounds, velocities(2),'satelite');
sat3 = satelite3D(20, 20, earth_radius+50, 100, rounds, velocities(3),'satelite');
sat4 = satelite3D(20, 20, earth_radius+50, 140, rounds, velocities(4),'satelite');
sat5 = satelite3D(20, 20, earth_radius+50, 180, rounds, velocities(5),'satelite'); 
sat6 = satelite3D(20, 20, earth_radius+50, 220, rounds, velocities(6),'satelite');
sat7= satelite3D( 20, 20, earth_radius+50, 260, rounds, velocities(7),'satelite');
sat8 = satelite3D(20, 20, earth_radius+50, 300, rounds, velocities(8),'satelite');
sat9 = satelite3D(20, 20, earth_radius+50, 340, rounds, velocities(8),'satelite');
station1 = satelite3D(20,20,earth_radius,21, rounds, velocities(9),'station');
station2 = satelite3D(20,20,earth_radius,260, rounds, velocities(10),'station');

%%%% Experimental orbits (harder, non-uniform case)
% sat1 = satelite3D(30, 30,earth_radius+50, 20, rounds, velocities(1),'satelite'); % arg_vel is weirdly defined. The bigger it's value the slower the object moves (since it makes more steps to make a full circle)
% sat2 = satelite3D(50, 50, earth_radius+50, 60, rounds, velocities(2),'satelite');
% sat3 = satelite3D(-20, -20,earth_radius+50, 100, rounds, velocities(3),'satelite');
% sat4 = satelite3D(100, 100, earth_radius+50, 140, rounds, velocities(4),'satelite');
% sat5 = satelite3D(35, 35,earth_radius+50, 180, rounds, velocities(5),'satelite'); 
% sat6 = satelite3D(90, 90, earth_radius+50, 220, rounds, velocities(6),'satelite');
% sat7= satelite3D(-90, -90, earth_radius+50, 260, rounds, velocities(7),'satelite');
% sat8 = satelite3D(-40, -40, earth_radius+50, 300, rounds, velocities(8),'satelite');
% sat9 = satelite3D(-120, -120, earth_radius+50, 340, rounds, velocities(8),'satelite');
% station1 = satelite3D(210,210,earth_radius,20, rounds, velocities(9),'station');
% station2 = satelite3D(230,230,earth_radius,260, rounds, velocities(10),'station');

%%%% SHOWING ORBITS!!! PROBLEM WITH VARIABLE ASSIGNMENT
% x1 = sat1.lifetime_coordinates(:,1)';
% x1 = sph2cart(x1(1),x1(2),x1(3));
% x2 = sat2.lifetime_coordinates(:,1)';
% x2 = sph2cart(x2(1),x2(2),x2(3));
% x3 = sat3.lifetime_coordinates(:,1)';
% x3 = sph2cart(x3(1),x3(2),x3(3));
% x4 = sat4.lifetime_coordinates(:,1)';
% x4 = sph2cart(x4(1),x4(2),x4(3));
% x5 = sat5.lifetime_coordinates(:,1)';
% x5 = sph2cart(x5(1),x5(2),x5(3));
% x6 = sat6.lifetime_coordinates(:,1)';
% x6 = sph2cart(x6(1),x6(2),x6(3));
% x7 = sat7.lifetime_coordinates(:,1)';
% x7 = sph2cart(x7(1),x7(2),x7(3));
% x8 = sat8.lifetime_coordinates(:,1)';
% x8 = sph2cart(x8(1),x8(2),x8(3));
% x9 = sat9.lifetime_coordinates(:,1)';
% x9 = sph2cart(x9(1),x9(2),x9(3));
% plotCircle3D([0,0,0],[x1(1),x1(2),x1(3)],earth_radius+50);
% plotCircle3D([0,0,0],[x2(1),x2(2),x2(3)],earth_radius+50);
% plotCircle3D([0,0,0],[x3(1),x3(2),x3(3)],earth_radius+50);
% plotCircle3D([0,0,0],[x4(1),x4(2),x4(3)],earth_radius+50);
% plotCircle3D([0,0,0],[x5(1),x5(2),x5(3)],earth_radius+50);
% plotCircle3D([0,0,0],[x6(1),x6(2),x6(3)],earth_radius+50);
% plotCircle3D([0,0,0],[x7(1),x7(2),x7(3)],earth_radius+50);
% plotCircle3D([0,0,0],[x8(1),x8(2),x8(3)],earth_radius+50);
% plotCircle3D([0,0,0],[x9(1),x9(2),x9(3)],earth_radius+50);

%%%% WARNING! 't' must be always < min(velocities)(lesser from linspace) of all satelites and stations (avoiding index out of bounds error)
stop = 600;
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
   pause(0.01) %WARNING: all 'delete' (for graphic objects) functions must be after the 'pause' function 
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

%% Creating objective function and constraints:
nodes = [sat1 sat2 sat3 sat4 sat5 sat6 sat7 sat8 sat9 station1 station2];
n = length(current_coordinates(:,1));% n = number of total nodes. (Remember N:#satelites, M:#stations)
M = 2; %number of stations in graph
DISTANCES = zeros(n,n); %(N+M)x(N+M)
LINKS = zeros(n,n); %(N+M)x(N+M)^2
communication_range = 200;
for i = 1:n
    for j = 1:n
        DISTANCES(i,j) = euclidean_dist(current_coordinates(i,:),current_coordinates(j,:)); 
        if euclidean_dist(current_coordinates(i,:),current_coordinates(j,:)) <= communication_range %200 is arbitrary
            if strcmp(nodes(i).name,'satelite') && strcmp(nodes(j).name,'satelite')
                LINKS(i,j) = 1; %satelite to satelite link
            else
                LINKS(i,j) = -1; %satelite to station link, '-1' since stations are sinks in graph terms
            end
        else
            LINKS(i,j) = 0;
        end
        if (i == j) % || ((i >= n-M) && j <= i ) % if self or (if station and lower trianglular matrix because "stations sink data only one way")
           LINKS(i,j)=0; % do not link the node to itself
        end
    end
end
num_of_links = length(find(LINKS~=0))/2;

%%%% Making negative the symmetric elements
% for i = 1:n
%     for j = (i+1):n
%         LINKS(i,j)=LINKS(i,j)*(-1);
%     end
% end

%Each row concerns -> a satelite, last 2 rows -> stations
disp('|=================================================================================|')
disp('Distance matrix:------------------------------------------------------------------')
disp(DISTANCES)
disp('Link matrix:----------------------------------------------------------------------')
disp(LINKS)%should I include the diagonal elements (selfs)?
disp('-- 1 = satelite to satelite connection, -1 = else --------------------------------')
disp('|=================================================================================|')


%% Constructing fmincon parameters (objective function, contraints: Aeq, beq, A, b)
%%%% Constructing the function
% Generate capacities:
objective_function = @(xs)[zeros(1,num_of_links*2),ones(1,n)]*xs'; %xs is the optimization vector. xs = [x1,x2,...,x(links_num),x(links_num+1),...,x(2*links_num),s1,s2,...,sn]

% Creating Aeq (Kirchhoff's law matrix)
xij = zeros(num_of_links,3); % xij[value][parent node 1][parent node 2]
counter = 1;
for i = 1:n
    for j = i:n
       if LINKS(i,j) ~= 0 
          xij(counter,1) = LINKS(i,j);
          xij(counter,2) = i;
          xij(counter,3) = j;
          counter = counter + 1;
       end
    end
end
disp('|--[value][parent node 1][parent node 2]|---|')
disp(xij)
disp('|-------------------------------------------|')

A_kirchhoff = zeros(0,num_of_links*2+n);
for i = 1:n %1:11
    temp = xijvec(i,xij,n);
    if isempty(A_kirchhoff)
        A_kirchhoff = [A_kirchhoff , temp]; %column bind
    else
        A_kirchhoff = [A_kirchhoff ; temp]; %row bind
    end
end

disp('|Aeq1:-------------------------------------------|')
disp(A_kirchhoff)
disp('|------------------------------------------------|')

% Constructing A, beq
A_ = diag(ones(1,length(A_kirchhoff(1,:)))); % each optimization variable must be equal or greater than zero.
disp('|A_:-----------------------------------------------|')
disp(A_)
disp('|------------------------------------------------|')

b_eq = zeros(1,length(A_kirchhoff(:,1)));
b_ = zeros(1,length(A_(:,1)))+20;
disp('|b_eq = b :------------------------------------------------|')
disp(b_eq)
disp('|------------------------------------------------|')


%% Using fmincon!
%%%%--------- fmincon(fun,x0,A,b,Aeq,beq) % x0 is the initial point used by the optimizer
opt_results = fmincon(objective_function, 1:(num_of_links*2+n), A_, b_,A_kirchhoff,b_eq,zeros(size(1:(num_of_links*2+n))));%,zeros(size(1:(num_of_links*2+n)))
for i = 1:length(opt_results)
    if i<=num_of_links*2
        sprintf("x%i= %e",[i,opt_results(i)])
    else
        sprintf("s%i= %e",[i-num_of_links*2,opt_results(i)])
    end
end

%% ---------------------- Function definitions ----------------------
%{
--------------------------------------------------------------------
Function definitions in a script must appear at the end of the file.
Only the first fuction can be accessed by other files. The rest are
acounted as local functions. MATLAB does no support default values on 
function parameters.
--------------------------------------------------------------------
%}
% plotCircle3D([0,0,0],[10,0,0],3)
% axis equal

function out = euclidean_dist(vec1, vec2)
%Definition: euclidean_dist(vec1, vec2)
%Calculates euclidean distance of two vectors
    out = sqrt(sum((vec1 - vec2) .^ 2));
end

function out = xijvec(i,x,num_nodes)
% Is used at the making of the xij matrix (kirchhoff matrix)
% i = current position of outter iterator, x = xij, num_nodes = number of nodes
% Output is a vector of A_kirchhoff matrix's row (fianl_Aeq).
    temp = zeros(1,length(x)*2);
    for it = 1:length(x)
        if i <= 9  % if inter-satelite communication
            if x(it,2)== i || x(it,3)== i % if it has a parent_node_1 or parent_node_2
                if x(it,2)== i % if xij, meaning if it is the parent_node_1
                    temp(it) = x(it,1); % because (Aeq1 = - Aeq2) and final_Aeq[] 
                    temp(it+length(x)) = -(-x(it,1));
                else % sidelestis_xij = - sidelestis_xji
                    temp(it) = -x(it,1);
                    temp(it+length(x)) = -(-(-x(it,1)));
                end
            end
            
        else % else if i > 9 a.k.a if it a station node...
            if x(it,2) == i || x(it,3)== i
%                temp(it) = x(it,1); 
                temp(it) = 0; %... ignore one direction (i.e. xi->xj but not xi<=>xj)
                temp(it+length(x)) = -x(it,1); %-(-x(it,1)); 
            end
        end
 
    end
    % Adding si (every single one of them = -1)
    temp2 = diag(ones(1,num_nodes));
    out = [temp, -temp2(i,:)]; % should it be(/2)? NO!
end

function plotCircle3D(center,normal,radius) %#ok<DEFNU>
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
plot3(points(1,:),points(2,:),points(3,:),'r-');
end