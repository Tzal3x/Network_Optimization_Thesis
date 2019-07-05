function out = xijvec3(node, xij, num_nodes, NUMBER_OF_SATELLITES) 
% - Used in solve_part2.m -
% This function is used for the creation of Aeq matrix (see solve_part1 function).
% Each time it is called, one row of Aeq matrix is created (regarding node i).
%
% The output is a vector of Aeq matrix's row.
% IMPORTANT NOTE: there is some crucial information missing.
%
% Input parameters: 
% i = node == current position of outter iterator (used in main script at assemblying)
% x = xij where xij is a cell array. xij{i} contains value, parent_node_1, parent_node2 of link i 
% num_nodes = total number of nodes
% num_sats = total number of satellites
% ----------------------------------------------------------------------------------------------- 
    flows = zeros(1,length(xij)*2);
    sat_to_stat_positions = null(1,1); % hafies!
    for it =1:length(xij) % for every link, check if any edge is node i
        value = xij{it}(1);
        parent1 = xij{it}(2);
        parent2 = xij{it}(3);
        if value == -1 % if it is inter-satellite connection
            sat_to_stat_positions = [sat_to_stat_positions, it];
        end
        if parent1 == node  % if node i is parent 1 of the link
           flows(it) = 1;
           flows(length(xij)+it) = -1; 
        elseif parent2 == node % if node i is parent 2 of the link 
           flows(it) = -1;
           flows(length(xij)+it) = 1; 
        end
    end
    disp("inside xijvec3")
    disp(sat_to_stat_positions)
    flows(length(xij) + sat_to_stat_positions) = [];
    % Final step: out == [Aeq1 Aeq2 si]
    divergencies = diag(ones(1,num_nodes));
    if node > NUMBER_OF_SATELLITES
       divergencies = -divergencies; 
    end
    buffers = -diag(ones(1,num_nodes)); % current buffer
    out = [flows, -divergencies(node,:), buffers(node,:)]; % [flows_1st_half, flows_2nd_half, divergencies, current_buffer]
end
