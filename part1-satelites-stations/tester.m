%% Space 3D:--------------------------------------------------------
%{
 Note that this is not a very robust implementation of the creation of the
 network.
- - - - - - - - - - - - - - - - -
 My Debug tools: 
1) disp('[~DEBUG:]'+string(object_of_interest))% debug
2) pause(1000)% debug 
%}

cd C:\Users\User\Documents\GitHub\Network_Optimization_Thesis\part1-satelites-stations

% Main parameters: - - - - -- - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% I should consider adding THETA, PHI and ALTITUDE vectors on the
% "create_nodes" function.
NUMBER_OF_SATELITES = 9; % integer, default 9
NUMBER_OF_STATIONS = 2; % integer, default 2
RANDOM_VELOCITIES = false; % boolean, default false
INVERSE_VELOCITIES_SATEL = ones(1,NUMBER_OF_SATELITES) * 3; % smaller value -> faster
INVERSE_VELOCITIES_STATIONS = ones(1,NUMBER_OF_STATIONS) * 80; % larger value -> slower
STOP_AT_TIME = 10; % integer
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

nodes = create_nodes(NUMBER_OF_SATELITES, NUMBER_OF_STATIONS, INVERSE_VELOCITIES_SATEL, ... 
                     INVERSE_VELOCITIES_STATIONS, RANDOM_VELOCITIES); 

% Setting up figure display options: - - - - - - - - - - - - - - - - -
figure('Name','3D Simulation of satelite orbits');
title("satelites:"+string(NUMBER_OF_SATELITES)+...
      ", stations:"+string(NUMBER_OF_STATIONS)+...
      ", random velocities:"+string(RANDOM_VELOCITIES));
  
hold on; % keep plotting on the existing figure
earth_radius = 200;
axis([-earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200]); %setting up figure size
axis equal

%%%% WARNING! 't' must be always < min(velocities)(lesser from linspace) of all satelites and stations (avoiding index out of bounds error)
rounds = 10;
stop = STOP_AT_TIME; % 30
times = 10000;
[x,y,z] = sphere();%get coordinates in order to create 3D spheres. A sphere can be the earth, a satelite or a station.
for i = 1:times

   %Coordinates ---------------------------
   coords = [];
   for j = 1:length(nodes)
      coords = [coords; nodes(j).lifetime_coordinates(:,i)' ];
   end

   %Visualize ---------------------------
   earth = surf( earth_radius*x, earth_radius*y, earth_radius*z );
   displays = [];
   for j = 1:length(nodes)
      if j <= NUMBER_OF_SATELITES
        displays = [displays, surf(coords(j,1)+10*x,coords(j,2)+10*y,coords(j,3)+10*z)]; 
      else
        displays = [displays, surf(coords(j,1)+20*x,coords(j,2)+20*y,coords(j,3)+20*z)]; 
      end
   end
   
   % Placing here the 'stop' break for the figure to stay active with all
   % objects:
   if i == stop
       break
   end
   
   %Delete (to create the animation illusion)
   pause(0.01) %WARNING: all 'delete' (for graphic objects) functions must be after the 'pause' function 
   for j = 1:length(displays)
      delete(displays(j))
   end
   delete(earth)
end
% hold off %not necessary(?)
disp('--- TIME STOPPED ---')


%% Creating objective function and constraints:
disp('Creating objective function and constraints...')
n = length(nodes);% number of total nodes
DISTANCES = zeros(n,n); %(N+M)x(N+M) , N: #sats, M: #stats
LINKS = zeros(n,n); %(N+M)x(N+M)^2
communication_range = 200;
for i = 1:n
    for j = 1:n
        DISTANCES(i,j) = euclidean_dist(coords(i,:),coords(j,:)); 
        if euclidean_dist(coords(i,:),coords(j,:)) <= communication_range %200 is arbitrary
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
disp('---------- OPT_RESULTS ----------')
disp(opt_results)
disp('---------------------------------')


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% - - - - - - - - - - - - - - - - - - - - F U N C T I O N S - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

function out = create_nodes(num_satelites, num_stations, sat_inverse_vel, stat_inverse_vel, random_factor)
% Creates a list of nodes in the form of [satelite_objects .... station_objects]
%
% - num_satelites: number of satelites (set default as 9)
% - num_stations: number of stations (set default as 2)
% - sat_inverse_vel: vector, inverse velocity of satelites (smaller values imply high
% velocity) use values around [3 3 3 ... 3]
% - stat_inverse_vel: vector, inverse velocity of satelites (every station has the
% same value) (smaller values imply high velocity) use values around [80 80 80 ... 80]
% - random_factor: boolean, if true it adds some random noise from uniform distribution in [-10,10] on the satelite velocities

    earth_radius = 200;
    rounds = 10;
    matr = [ ];
    inverse_velocities = [sat_inverse_vel, stat_inverse_vel]; % the bigger the value the slower the velocity
    
    if random_factor
        rng(79) %setting seed
        random_numbers_sats = (rand(1,num_satelites)*8+1)-(4-1);%(rand(1,num_satelites)*number*2)-number -> rand() returns values in (0,1) so
                                                                % I do this in order to turn it into (-number,number)
        rands = [random_numbers_sats, zeros(1,num_stations)];
        
        inverse_velocities = inverse_velocities + rands;
    end
    inverse_velocities = inverse_velocities * 1000;
    
    % Remember: satelite3D(arg_theta,arg_phi, arg_alt, arg_init_pos, arg_periods, arg_vel, arg_name)
    % Constructing satelites & stations:
    initializer_difs = (0:(num_satelites+num_stations-1))*27; %used to initialiaze satelites at different positions in order to not overlap each other
    for ii = 1:(num_satelites + num_stations)
        if ii <= num_satelites
            matr = [matr, satelite3D(20, 20, earth_radius+50, 20 + initializer_difs(ii), rounds, inverse_velocities(ii), 'satelite')];
        else
            matr = [matr, satelite3D(20, 20, earth_radius, 200 + initializer_difs(ii)*4, rounds, (num_satelites+num_stations)*1000, 'station')];
        end
    end
    
    
    disp('[~Report:] Satelites and stations successfully created')
    out = matr;
end

function out = euclidean_dist(vec1, vec2)
% Definition: euclidean_dist(vec1, vec2)
% Calculates euclidean distance of two vectors
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


