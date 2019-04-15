%% Space 3D:--------------------------------------------------------
%{
 Note that this is not a very robust implementation of the creation of the
 network. Only a simple version (everyones orbit exists in the same circle) of the problem is being tested.
 I should consider adding THETA, PHI and ALTITUDE vectors on the
 "create_nodes" function.
%}

cd C:\Users\User\Documents\GitHub\Network_Optimization_Thesis\part1-satelites-stations

% Main parameters: - - - - -- - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
NUMBER_OF_SATELLITES = 9; % integer, default 9
NUMBER_OF_STATIONS = 2; % integer, default 2
RANDOM_VELOCITIES = false; % boolean, default false
INVERSE_VELOCITIES_SATEL = ones(1,NUMBER_OF_SATELLITES) * 3; % smaller value -> faster, it can be a vector of the wanted speeds [v1 v2 ... vn], where n == NUMBER_OF_SATELLITES
INVERSE_VELOCITIES_STATIONS = ones(1,NUMBER_OF_STATIONS) * 80; % larger value -> slower, >> >> >> >> >> >> >> >> >> >> >> >> 
STOP_AT_TIME = 300; % integer, declares when the time should be stopped
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

nodes = create_nodes(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, INVERSE_VELOCITIES_SATEL, ... 
                     INVERSE_VELOCITIES_STATIONS, RANDOM_VELOCITIES); 

% Setting up figure display options: - - - - - - - - - - - - - - - - -
figure('Name','3D Simulation of satelite orbits');
title("satelites:"+string(NUMBER_OF_SATELLITES)+...
      ", stations:"+string(NUMBER_OF_STATIONS)+...
      ", random velocities:"+string(RANDOM_VELOCITIES)+...
      ", stop:"+string(STOP_AT_TIME));
  
hold on; % keep plotting on the existing figure
earth_radius = 200;
axis([-earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200]); %setting up figure size
axis equal
view(30,0)% setting azemuth and elevation angles of camera for nicer visualization 

%%%% WARNING! 'times' must be always < lesser from linspace length of all satelites and stations (avoiding index out of bounds error)
rounds = 10;
stop = STOP_AT_TIME; % 30
times = 10000; % upper bound of iterations
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
      if j <= NUMBER_OF_SATELLITES
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
n = NUMBER_OF_SATELLITES + NUMBER_OF_STATIONS;% number of total nodes
DISTANCES = zeros(n,n); %(N+M)x(N+M) , N: #sats, M: #stats
LINKS = zeros(n,n); %(N+M)x(N+M)^2
communication_range = 150; % links exist at a distance smaller or equal of communication_range
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
% xij = zeros(num_of_links,3); % xij[value][parent node 1][parent node 2]
counter = 1;
xij={}; % cell array (like lists in the R programming language)
for i = 1:n
    for j = i:n
       if LINKS(i,j) ~= 0 
          xij{counter} = 1:3;
          xij{counter}(1) = LINKS(i,j);
          xij{counter}(2) = i;
          xij{counter}(3) = j;
          counter = counter + 1;
       end
    end
end
disp('|[value][parent node 1][parent node 2]------|')
for i = 1:(counter-1)
%     disp(string(xij{i}(1)) + "|" + string(xij{i}(2)) + "|" + string(xij{i}(3)) + "|")
   disp(xij{i})
end
disp('|-------------------------------------------|')

A_kirchhoff = zeros(0,num_of_links*2+n);
for i = 1:n
    temp = xijvec(i,xij,n,NUMBER_OF_SATELLITES,NUMBER_OF_STATIONS);
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
A = diag(ones(1,length(A_kirchhoff(1,:)))); % each optimization variable must be equal or greater than zero.
% disp('|A:-----------------------------------------------|')
% disp(A)
% disp('|------------------------------------------------|')

b_eq = zeros(1,length(A_kirchhoff(:,1)));
b = zeros(1,length(A(:,1)))+20;
% disp('|b_eq = b :------------------------------------------------|')
% disp(b_eq)
% disp('|------------------------------------------------|')



%% Using fmincon!
%%%%--------- fmincon(fun,x0,A,b,Aeq,beq) % x0 is the initial point used by the optimizer
opt_results = fmincon(objective_function, 1:(num_of_links*2+n), A, b, A_kirchhoff, b_eq, zeros(size(1:(num_of_links*2+n)))); % zeros(size(1:(num_of_links*2+n))) -> every rate should be positive
% for i = 1:length(opt_results) % Print in a 'nicer' format. Only prints non zero values
%     if opt_results(i)>1
%         if i<=num_of_links*2
%             sprintf("x%i= %e",[i,opt_results(i)])
%         else
%             sprintf("s%i= %e",[i-num_of_links*2,opt_results(i)])
%         end    
%     end
% end
disp('---------- OPT_RESULTS ----------------------------------------------------------------------------')
disp(opt_results)
disp('---------------------------------------------------------------------------------------------------')

%% Plot 2D world map and links
%{
 THIS SCRIPT IS FOR THE VISUALIZATION OF THE SATELITE ORBITS ON A 2D MAP
 OF THE EARTH. BEFORE CONTINUING IT SHOULD BE CONSIDERED THAT THERE IS NO 
 COMPLETE MATCHING OF THE 3D ORBITS TO THE 2D MAPPING BECAUSE STATIONS 
 SHOULD BE REGARDED AS STATIONARY (ROTATION OF THE EARTH WILL BE IGNORED) 
 TO KEEP THE VELOCITY CALCULATIONS OF EACH SATELITE SIMPLE.
%}




%% FUNCTIONS
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
        random_numbers_sats = (rand(1,num_satelites)*8+1)-(4-1);% (rand(1,num_satelites)*number*2)-number -> rand() returns values in (0,1) so
                                                                % I do this in order to turn it (the output range) into (-number,number)
        rands = [random_numbers_sats, zeros(1,num_stations)];% adding random factor only at satelite velocities
        
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
            matr = [matr, satelite3D(20, 20, earth_radius, 200 + initializer_difs(ii)*4, rounds, inverse_velocities(ii), 'station')];
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


%Another version of xijvec that was working correctly, returning results
%of 13 instead of 17 on the same optimization variables.

function out = xijvec(i,x,num_nodes,num_sats,num_stats)
% This function is used for the creation of kirchoff matrix.
% However, each time it is called, one row of kirchoff matrix is created
% (regarding link i).
% Therefore, the output is a vector of A_kirchhoff matrix's row (fianl_Aeq).
% Parameters: - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% i = current position of outter iterator
% x = xij where xij is a cell array. xij{i} contains value, parent_node_1, parent_node2 of link i 
% num_nodes = number of nodes
% num_sats = number of satelites

    temp = zeros(1,length(x)*2);
    for it = 1:length(x)
        if  x{it}(2)<(num_sats+num_stats) && x{it}(3)<(num_sats+num_stats)% if inter-satelite communication
            if x{it}(2)==i || x{it}(3)==i % if i is parent_node_1 or parent_node_2
                if x{it}(2)==i % if it is the parent_node_1:
                    temp(it) = x{it}(1); % get it's value
                    temp(it+length(x)) = -(x{it}(1));% (Aeq1 = - Aeq2) <=> 
                else % if it is the parent_node_2:
                    temp(it) = x{it}(1);
                    temp(it+length(x)) = -(-x{it}(1));
                end
            end
            
        else % a.k.a if its station-station link...
            if x{it}(2)==i || x{it}(3)==i
               temp(it) = x{it}(1); 
               temp(it+length(x)) = 0;%... ignore one direction (i.e. xi->xj but not xi<=>xj)
            else
                temp(it) = 0; 
                temp(it+length(x)) = 0;
            end
        end
 
    end
    % Adding si: (where every si = -1)
    temp2 = diag(ones(1,num_nodes));
    out = [temp, -temp2(i,:)];
end


% function out = xijvec(i,x,num_nodes,num_sats,num_stats)
% % This function is used for the creation of kirchoff matrix.
% % However, each time it is called, one row of kirchoff matrix is created
% % (regarding link i).
% % Therefore, the output is a vector of A_kirchhoff matrix's row (fianl_Aeq).
% % Parameters: - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% % i = current position of outter iterator
% % x = xij where xij is a cell array. xij{i} contains value, parent_node_1, parent_node2 of link i 
% % num_nodes = number of nodes
% % num_sats = number of satelites
% 
%     temp = zeros(1,length(x)*2);
%     for it = 1:length(x)
%         if  x{it}(2)<(num_sats+num_stats) && x{it}(3)<(num_sats+num_stats)% if inter-satelite communication
%             if x{it}(2)==i || x{it}(3)==i % if i is parent_node_1 or parent_node_2
%                 if x{it}(2)==i % if it is the parent_node_1:
%                     temp(it) = x{it}(1); % get it's value
%                     temp(it+length(x)) = -(-(x{it}(1)));% (Aeq1 = - Aeq2) <=> 
%                 else % if it is the parent_node_2:
%                     temp(it) = -x{it}(1);
%                     temp(it+length(x)) = -(-(-x{it}(1)));
%                 end
%             end
%             
%         else % if its a station to satelite link...
%             if x{it}(2)==i || x{it}(3)==i
%                temp(it) = x{it}(1); 
%                temp(it+length(x)) = 0;%... ignore one direction (i.e. xi->xj but not xi<=>xj)
%             else
%                 temp(it) = -x{it}(1); 
%                 temp(it+length(x)) = 0;
%             end
%         end
%  
%     end
%     % Adding si: (where every si = -1)
%     temp2 = diag(ones(1,num_nodes));
%     out = [temp, -temp2(i,:)];
% end