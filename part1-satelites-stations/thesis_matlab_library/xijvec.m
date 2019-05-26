function out = xijvec(i,x,num_nodes,num_station_links,NUMBER_OF_SATELLITES) 
% - Used in solve_part1.m -
% This function is used for the creation of A_kirchoff matrix (see solve_part1 function).
% Each time it is called, one row of kirchoff matrix is created
% (regarding link i).
% The output is a vector of A_kirchhoff matrix's row (fianl_Aeq).
%
% Input parameters: 
% i = current position of outter iterator (used in main script at assemblying)
% x = xij where xij is a cell array. xij{i} contains value, parent_node_1, parent_node2 of link i 
% num_nodes = number of nodes
% num_sats = number of satelites
% ----------------------------------------------------------------------------------------------- 
    temp = zeros(1,length(x));
    for it =1:length(x) % for every link, check if any edge is node i
        parent1 = x{it}(2);
        parent2 = x{it}(3);
        if parent1 == i  % if node i is parent 1 of the link
           temp(it) = 1;
        else %elif:
            if parent2 == i % if node i is parent 2 of the link
                temp(it) = -1;
            end
        end
    end
    % Final step: out == [Aeq1 Aeq2 si]
    temp3 = diag(ones(1,num_nodes));
    if i > NUMBER_OF_SATELLITES
       temp3 = -temp3; 
    end
    temp2 = temp(1:(length(temp)-num_station_links)); % se periptwsh pou den valoume tis sthles
    out = [temp, -temp2, -temp3(i,:)];
end
