function out = create_Aeq_row(epoch, epoch_link_list, node, xij, num_nodes, num_station_links, NUMBER_OF_SATELLITES)
% - Used in solve_part2.m -
% This function is used for the creation of A_kirchoff matrix (see solve_part1 function).
% Each time it is called, one row of kirchoff matrix is created
% (regarding link i).
% The output is a vector of A_kirchhoff matrix's row (fianl_Aeq).
%
% Input parameters: 
% epoch_link_list = a list contating the number of links per epoch.
% node = node i, current position of outter iterator (used in main script
% at Aeq creation/assemblying).
% xij = a cell array. xij{i} contains: value, parent_node_1, parent_node2 of link i. The value is 1 if it is an intersatellite connection else -1 
% num_nodes = number of nodes.
% num_sats = number of satelites.
% ----------------------------------------------------------------------------------------------- 
    tnumflows = sum(epoch_link_list); % total number of flows
    tnumepochs = length(epoch_link_list); % total number of epochs
    final_result = zeros(1,tnumflows + tnumepochs*num_nodes*2); % *2 epeidh length_of([divergencies, buffers] = 2*n)
    %______________________________________________________________________________________________________________
    flows1 = zeros(1,length(xij));
    for i =1:length(xij) % for every link, check if any edge is node i
        parent1 = xij{i}(2);
        parent2 = xij{i}(3);
        if parent1 == node  % if node i is parent 1 of the link
           flows1(i) = 1;
        elseif parent2 == node % if node i is parent 2 of the link
           flows1(i) = -1;
        end
    end
    divergencies = diag(ones(1,num_nodes));
    divergencies = -divergencies(node,:);  
    if node > NUMBER_OF_SATELLITES
       divergencies = -divergencies; 
    end
    flows2 = -flows1(1:(length(flows1)-num_station_links)); % se periptwsh pou den valoume tis sthles
    %______________________________________________________________________________________________________________
%     disp('flows1')% debug
%     disp(flows1)% debug
%     disp('flows1')% debug
%     disp(flows2)% debug
%     disp('divergencies')% debug
%     disp(divergencies)% debug
    
    if epoch == 1
        vima = 0; 
    else
        vima = epoch_link_list(epoch-1); % arithmos links prohgoumenis epoxis
    end
    
    disp('epoch:'+string(epoch)) % debug
    
    for i = (vima + 1):(vima + length(flows1) + length(flows2) + length(divergencies) + 2)
        if i <= vima + length(flows1) % adding flows 1
            final_result(i) = flows1(i);
        elseif i <= vima + length(flows1) + length(flows2) % adding flows 2 
            final_result(i) = flows2(i);
        elseif i <= vima + length(flows1) + length(flows2) + length(divergencies)
            final_result(i) = divergencies(i);
        else % adding buffers
            if epoch == 1 % I do notn add si(t-1) since s(0)==0
                final_result(vima + length(flows1) + length(flows2) + length(divergencies) + node) = -1; % si(t)            
            else
                final_result(vima + length(flows1) + length(flows2) + length(divergencies) + node) = -1; % si(t)
                final_result( vima + 1 - 1 - num_nodes + node) = 1; % si(t-1)
            end
        end
    end
    out = final_result;
end