function out = xijvec(node,x,num_nodes,NUMBER_OF_SATELLITES) 
% - Used in solve_part1.m -
% This function is used for the creation of A_kirchoff matrix (see solve_part1 function).
% Each time it is called, one row of kirchoff matrix is created
% (regarding link i).
% The output is a vector of A_kirchhoff matrix's row (fianl_Aeq).
%
% Input parameters: 
% i = link == current position of outter iterator (used in main script at assemblying)
% x = xij where xij is a cell array. xij{i} contains value, parent_node_1, parent_node2 of link i 
% num_nodes = number of nodes
% num_sats = number of satelites
% ----------------------------------------------------------------------------------------------- 
    temp = zeros(1,length(x)*2);
    sat_to_stat_positions = null(1,1); % hafies!
    for it =1:length(x) % for every link, check if any edge is node i
        value = x{it}(1);
        parent1 = x{it}(2);
        parent2 = x{it}(3);
        if value == -1 % if it is inter-satellite connection
            sat_to_stat_positions = [sat_to_stat_positions, it];
        end
        if parent1 == node  % if node i is parent 1 of the link
           temp(it) = 1;
           temp(length(x)+it) = -1; 
        elseif parent2 == node % if node i is parent 2 of the link
           temp(it) = -1;
           temp(length(x)+it) = 1; 
        end
    end
    temp(length(x) + sat_to_stat_positions) = [];
    temp3 = diag(ones(1,num_nodes)); % divergencies
    if node > NUMBER_OF_SATELLITES
       temp3 = -temp3; 
    end
    out = [temp, -temp3(node,:)];
end
