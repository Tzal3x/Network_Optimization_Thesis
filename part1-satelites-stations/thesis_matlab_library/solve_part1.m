function out = solve_part1(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, RANDOM_VELOCITIES, INVERSE_VELOCITIES_SATEL,...
                    INVERSE_VELOCITIES_STATIONS, STOP_AT_TIME, THETA_PHI, LINK_CAPACITY, PRINT_DETAILS, PRINT_MAIN_PARAMETERS,SHOW_TOPOLOGY)
    % Solves the first (and simplest) version of a network flow utility maximization problem 
    %
    % Input parameters:
    % NUMBER_OF_SATELLITES = integer, number of satellites 
    % NUMBER_OF_STATIONS = integer, number of stations 
    % RANDOM_VELOCITIES = boolean, true if a random factor should be added at the satellite velocities
    % INVERSE_VELOCITIES_SATEL = vector, each element is the inverse velocity of a satellite
    % INVERSE_VELOCITIES_STATIONS = vector, each element is the inverse
    % velocity of a station. Shouldl be the same for all, since it coincide
    % with earths angular speed around itself
    % STOP_AT_TIME = integer, time when the optimization should take place
    % THETA_PHI = vector, consisting of 2 elements [azimuth elevation]
    % LINK_CAPACITY = float, maximum link capacity
    % PRINT_DETAILS = boolean
    %
    % Output:
    % A vector consisting of optimal values
    % ---------------------------------------------------------------------
   
    disp("|========================================= INSIDE PART 1 =============================================|")
    nodes = create_nodes(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, INVERSE_VELOCITIES_SATEL, ... 
                         INVERSE_VELOCITIES_STATIONS, RANDOM_VELOCITIES, THETA_PHI); 

    if PRINT_MAIN_PARAMETERS
       disp("NUMBER_OF_SATELLITES: " + string(NUMBER_OF_SATELLITES))
       disp("NUMBER_OF_STATIONS: " + string(NUMBER_OF_STATIONS))
       disp("RANDOM_VELOCITIES: " + string(RANDOM_VELOCITIES))
       disp("INVERSE_VELOCITIES_SATEL: " + print_list(INVERSE_VELOCITIES_SATEL))
       disp("INVERSE_VELOCITIES_STATIONS: " + print_list(INVERSE_VELOCITIES_STATIONS))
       disp("STOP_AT_TIME: " + string(STOP_AT_TIME))
       disp("THETA_PHI: " + print_list(THETA_PHI))
       disp("LINK_CAPACITY: " + string(LINK_CAPACITY))
       disp("PRINT_DETAILS: " + string(PRINT_DETAILS))
    end

    % Setting up figure display options: - - - - - - - - - - - - - - - - -
    earth_radius = 200;
    if SHOW_TOPOLOGY
    figure('Name','3D Simulation of satelite orbits');
    title("satelites:"+string(NUMBER_OF_SATELLITES)+...
          ", stations:"+string(NUMBER_OF_STATIONS)+...
          ", random velocities:"+string(RANDOM_VELOCITIES)+...
          ", stop:"+string(STOP_AT_TIME));
        hold on; % keep plotting on the existing figure
        axis equal
        view(30,0)% setting azemuth and elevation angles of camera for nicer visualization 
    axis([-earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200]); %setting up figure size
    end    
    
    %%%% WARNING! 'times' must be always < lesser from linspace length of all satelites and stations (avoiding index out of bounds error)
    stop = STOP_AT_TIME; % 30
    times = 10000; % upper bound of iterations
    [x,y,z] = sphere();%get coordinates in order to create 3D spheres. A sphere can be the earth, a satelite or a station.   
    for i = 1:times
       %Coordinates ---------------------------
       coords = [];
       for j = 1:length(nodes)
          coords = [coords; nodes(j).lifetime_coordinates(:,i)' ];
       end
       
       if SHOW_TOPOLOGY
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
       end
           % Placing here the 'stop' break for the figure to stay active with all
           % objects:
           if i == stop
               break
           end
       if SHOW_TOPOLOGY
           %Delete (to create the animation illusion)
           pause(0.01) %WARNING: all 'delete' (for graphic objects) functions must be after the 'pause' function 
           for j = 1:length(displays)
              delete(displays(j))
           end
           delete(earth)
       end
    end
    % hold off %not necessary(?)
    disp('[~Report:] Time stopped! A network has been created.')

    % Creating objective function and constraints: ========================================================================================================================================================================================================================================================
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
    disp('Done!')
    %Each row concerns -> a satelite, last 2 rows -> stations

    if PRINT_DETAILS
        disp('|=================================================================================================|')
        disp('Distance matrix:------------------------------------------------------------------')
        disp(DISTANCES)
        disp('Link matrix:----------------------------------------------------------------------')
        disp(LINKS)%should I include the diagonal elements (selfs)?
        disp('-- [1] = satelite to satelite connection, [-1] = satelite to station -------------')
        disp('-- Number of links: '+string(num_of_links))
        disp('|=================================================================================================|')
    end

    % Constructing fmincon parameters (objective function, contraints: Aeq, beq, A, b) ========================================================================================================================================================================================================================================================
    %%%% Constructing the function
    disp('Constructing optimization parameters...')

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
    if PRINT_DETAILS
        disp('|[value][parent node 1][parent node 2]------|')
        for i = 1:(counter-1)
        %     disp(string(xij{i}(1)) + "|" + string(xij{i}(2)) + "|" + string(xij{i}(3)) + "|")
           disp(xij{i})
        end
        disp('') % break line
    end

    num_station_links = 0; % number of links that have a station as a parent
    for i = 1:length(xij)
        if xij{i}(2) > NUMBER_OF_SATELLITES || xij{i}(3) > NUMBER_OF_SATELLITES
            num_station_links = num_station_links + 1;
        end
    end

    A_kirchhoff = null(1,1); % Each row is a node. Each column until column[num of satellites] expresses which links belong to the node. Multiplying A_kirchhoff with optimization vector creates all the kirchhoff equalities that must hold.
    for i = 1:n
        temp = xijvec(i,xij,n,NUMBER_OF_SATELLITES);%,NUMBER_OF_STATIONS); % DEBUG , I WAS PREVIOUSLY USING xijvec()
        A_kirchhoff = [A_kirchhoff ; temp]; %row bind
    end
    %A_kirchhoff = [A_kirchhoff ; [zeros(1,num_of_links*2), ones(1,NUMBER_OF_SATELLITES + NUMBER_OF_STATIONS)]]; % Adding sum(divergencies)==0 (NOT NECESSARY, IT SHOULD HOLD NEVERTHELESS)

    beq = zeros(1,length(A_kirchhoff(:,1))); % OTHER KIRCHOFF EQUALITY PART
    A = [];% A = diag(ones(1,length(A_kirchhoff(1,:))));
    b = [];% b = [zeros(1,num_of_links*2-num_station_links)+LINK_CAPACITY,zeros(1,NUMBER_OF_SATELLITES)+inf, zeros(1,NUMBER_OF_STATIONS)]; %zerosNUMBER_OF_SATELLITES+100 % CAPACITIES (UPPER BOUNDS)
    upper_bounds = [zeros(1,num_of_links*2-num_station_links)+LINK_CAPACITY, zeros(1,NUMBER_OF_SATELLITES+NUMBER_OF_STATIONS)+inf];
    % lower_bounds = [zeros(1,(num_of_links*2-num_station_links)) zeros(1,NUMBER_OF_SATELLITES)+lbsi (1:NUMBER_OF_STATIONS)*(-inf)]; % zeros for rates (and sources(?)), -inf for sinks
    lower_bounds = zeros(1,length(A_kirchhoff(1,:)));
    disp('Done!')
    % Prints: ------------------------------------------------------------
    if PRINT_DETAILS
        disp('|- - -------\ (Aeq) K I R C H H O F F /-------------- - -|')
        disp(A_kirchhoff)
        disp('|- - -------\ (beq) EQUALITY VECTOR /---------------- - -|')
        disp(beq)
    %     disp('|- - -------\ (A) INEQUALITY MATRIX /---------------- - -|')
    %     disp(A)
        disp('|- - -------\ (b) CAPACITIES (UPPER BOUNDS)/--------- - -|')
    %     disp(b)
        disp(upper_bounds)
        disp('|- - -------\ (lower_bounds) LOWER BOUNDS /---------- - -|')
        disp(lower_bounds)
    end
    % --------------------------------------------------------------------
    % Using fmincon! ========================================================================================================================================================================================================================================================
    %%%%--------- fmincon(fun,x0,A,b,Aeq,beq) % x0 is the initial point used by the optimizer
    objective_function = @(xs)[(10^(-3))*ones(1,num_of_links*2 - num_station_links), zeros(1,NUMBER_OF_SATELLITES),-ones(1,NUMBER_OF_STATIONS)]*xs'; %xs is the optimization vector. xs = [x1,x2,...,x(links_num),x(links_num+1),...,x(2*links_num),s1,s2,...,sn]
    x0 = zeros(1,(num_of_links*2-num_station_links+n));
    % x0 = [ones(1,(num_of_links*2-num_station_links + n - 3)) 1 2 -100]; % debug
    opt_results = fmincon(objective_function, x0, A, b', A_kirchhoff, beq', lower_bounds, upper_bounds); % zeros(size(1:(num_of_links*2+n))) -> every rate should be positive
    % for i = 1:length(opt_results) % Print in a 'nicer' format. Only prints non zero values
    %     if opt_results(i)>1
    %         if i<=num_of_links*2
    %             sprintf("x%i= %e",[i,opt_results(i)])
    %         else
    %             sprintf("s%i= %e",[i-num_of_links*2,opt_results(i)])
    %         end    
    %     end
    % end
    disp('---------- OPTIMIZATION RESULTS ----------------------------------------------------------------')
    disp(opt_results)
    disp('------------------------------------------------------------------------------------------------')
    out = opt_results;
    
    %%% Create Simple Graph (figure): _____________________________________
    createGraph(NUMBER_OF_SATELLITES, xij, opt_results, nodes, 0)
    
    % Create Graph on Earth's map (figure):________________________________
    
end% end of main function