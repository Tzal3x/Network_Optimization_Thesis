%{
!> PART 2: Formulating a problem like the "problem 1" of the paper:
" Delay-Tolerant Network Utility Maximization, Toumpis et al."
*** [Self-note:] A 'cd C:/...' may be needed for the function to run smoothly *** 
%}

function out = solve_part2(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, RANDOM_VELOCITIES, INVERSE_VELOCITIES_SATEL,...
                    INVERSE_VELOCITIES_STATIONS, STOP_AT_TIME, THETA_PHI, LINK_CAPACITY, PRINT_DETAILS, PRINT_MAIN_PARAMETERS,SHOW_TOPOLOGY)
    % Solves the second version of the problem Delay-Tolerant Network
    % Utility Maximization Problem [Problem I: DTNUM]
    % Having added data buffers
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
    % (number of epochs)
    % THETA_PHI = vector, consisting of 2 elements [azimuth elevation]
    % LINK_CAPACITY = float, maximum link capacity
    % PRINT_DETAILS = boolean
    %
    % Output:
    % A vector consisting of optimal values
    % ---------------------------------------------------------------------
   
    disp("|========================================= INSIDE PART 2 =============================================|")
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
    
    n = NUMBER_OF_SATELLITES + NUMBER_OF_STATIONS; % number of total nodes (#sats + #stats)

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
    stop = STOP_AT_TIME; % epochs, keep it low like 30
    times = 10000; % upper bound of iterations
    [x,y,z] = sphere();%get coordinates in order to create 3D spheres. A sphere can be the earth, a satelite or a station.       
%     epoch_link_list = get_list_num_links(STOP_AT_TIME, nodes); % number of links per epoch 
    
%===================================================================================================================================================
%===================================================================================================================================================
%===================================================================================================================================================
%===================================================================================================================================================
%===================================================================================================================================================
    for epoch = 1:times
       %Coordinates ---------------------------
       coords = [];
       for j = 1:length(nodes)
          coords = [coords; nodes(j).lifetime_coordinates(:,epoch)' ];
       end
       
       LINKS = abs(create_LINKS(coords, nodes, n));
%        xij = create_flow_info(LINKS, n); 
%        num_station_links = get_nstat_links(xij, NUMBER_OF_SATELLITES);    
       
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
           if epoch == stop
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
%===================================================================================================================================================
%===================================================================================================================================================
%===================================================================================================================================================
%===================================================================================================================================================
%===================================================================================================================================================
%===================================================================================================================================================
    disp('[~Report:] Final epoch reached!')   
    if PRINT_DETAILS
        disp('|=================================================================================================|')
        disp('Link matrix:----------------------------------------------------------------------')
        disp(abs(LINKS))%should I include the diagonal elements (selfs)?
        disp('|=================================================================================================|')
    end
    xij = create_flow_info(LINKS, n); % get xij (parent1 parent2 vector)
    if PRINT_DETAILS
        disp('|[value][parent node 1][parent node 2]------|')
        for epoch = 1:length(xij)
        %     disp(string(xij{i}(1)) + "|" + string(xij{i}(2)) + "|" + string(xij{i}(3)) + "|")
           disp(xij{epoch})
        end
        disp('') % break line
    end
    disp('[~Report:] Constructing optimization parameters...')
    aeq_epochs={}; %contains the rows of Aeq, each element is an xijvec vector % MIGHT REMOVE LATER
    counter = 1;
    for epoch = 1:STOP_AT_TIME % construct sub-matrix of aeq for each epoch
        for node = 1:n
           aeq_epochs{counter} = xijvec2(epoch, node, STOP_AT_TIME, xij, nodes, NUMBER_OF_SATELLITES);
           counter = counter + 1;
        end
    end
    Aeq = assemble_Aeq(aeq_epochs);
    if PRINT_DETAILS
       disp('|- - -------------\ (Aeq)  /-------------- - -|')
       disp(Aeq)
    end
%===================================================================================================================================================
%===================================================================================================================================================
    % Constructing fmincon parameters (objective function, contraints: Aeq, beq, A, b) 
    %%%% Constructing the function
    beq = zeros(1,length(Aeq(:,1))); % OTHER KIRCHOFF EQUALITY PART
    A = [];
    b = [];
    
    num_links = length(xij); % number of unidirectional links
    num_station_links = 0; % number of links that have a station as a parent
    for it = 1:length(xij)
        if xij{it}(2) > NUMBER_OF_SATELLITES || xij{it}(3) > NUMBER_OF_SATELLITES
            num_station_links = num_station_links + 1;
        end
    end
    % At upper_bounds vector, flows & buffers should have a bound. Divergencies not
    upper_bounds = create_upper_bounds(num_links*2 - num_station_links, n, STOP_AT_TIME);
    lower_bounds = zeros(1,length(Aeq(1,:))); % remains the same as in PART_1
    
    disp('[~Report:] Done!')
    if PRINT_DETAILS
        disp('|- - -------\ (beq) EQUALITY VECTOR /---------------- - -|')
        disp(beq)
        disp('|- - -------\ UPPER BOUNDS /--------- - -|')
        disp(upper_bounds)
        disp('|- - -------\ LOWER BOUNDS /---------- - -|')
        disp(lower_bounds)
    end

    % Using fmincon! ============================================================================================================================================================================================================================
    %%%%--------- fmincon(fun,x0,A,b,Aeq,beq) % x0 is the initial point used by the optimizer
    objective_function = @(xs)objective_function_coefs(num_links*2 - num_station_links, n, STOP_AT_TIME, NUMBER_OF_SATELLITES)*xs'; %xs is the optimization vector. xs = [x1,x2,...,x(links_num),x(links_num+1),...,x(2*links_num),s1,s2,...,sn]
    x0 = rand(1,length(Aeq(1,:)));
    opt_results = fmincon(objective_function, x0, A, b', Aeq, beq', lower_bounds, upper_bounds); % zeros(size(1:(num_of_links*2+n))) -> every rate should be positive
    disp('---------- OPTIMIZATION RESULTS ----------------------------------------------------------------')
    disp(opt_results)
    disp('------------------------------------------------------------------------------------------------')
    out = opt_results;
end% end of main function