% START: TESTING create_Aeq_row
% epochs = 1:3;
% nodes = 1:3;
% LINKS = [[ 0 1 1 ];[1 0 1];[1 1 0]];
% for i = 1:length(epochs)
%     for j = 1:length(nodes)
%         if i == 1 && j == 1
%             Aeq = create_Aeq_ro(i,length(epochs),j,LINKS, 1);
%         else
%             Aeq = [Aeq ; create_Aeq_ro(i,length(epochs),j,LINKS, 1)]; 
%         end
%     end
% end
% disp(Aeq)
% END: TESTING create_Aeq_row
disp(create_Aeq_ro(1,3,2,LINKS, 1))
function out = create_Aeq_ro(epoch, total_epochs, node, LINKS, NUMBER_OF_SATELLITES) % ALLAKSE THS TO ONOMA
% IT IS WRONG! DO NOT USE IT! - Was about to be used in solve_part2.m -  
    num_nodes = length(LINKS(1, :));
    final_result = (1:(((num_nodes-1)*2 + num_nodes*2 ) * total_epochs)) *0;
    % Estw o LINKS matrix:
    % [ 0 1 1 | Ston anw trigwniko vazw thetikous sydelestes twn rates
    % | 1 0 1 | Ston katw, arnhtikous. (bidirectional graph).
    % | 1 1 0 ] H ylopoihsh einai h apo katw:
    % FLOWS ---------------------------------------------------------------
    i = (epoch-1) * num_nodes ;
    % Vazodas sydelestes twn rates:
    for j = 1:num_nodes % iterating LINKS matrix
        if node <= num_nodes - NUMBER_OF_SATELLITES %(an einai satellite)
            if j < node
               final_result(i+j) = -LINKS(node, j); % arnhtikoi sydelestes (an yparxei link)
               final_result(i+j+num_nodes) = -final_result(i+j);
            elseif j > node
               final_result(i+j) = LINKS(node, j); % thetikoi sydelestes (an yparxei link)
               final_result(i+j+num_nodes) = -final_result(i+j);
            else
                continue; % den ginetai na exw link me ton eauto mou ara continue
            end
        elseif node >  num_nodes - NUMBER_OF_SATELLITES %(an einai station)
            if j < node
               final_result(i+j) = -LINKS(node, j); % arnhtikoi sydelestes (an yparxei link)
            else
               continue;
            end
        end % end if satellite/station
    end % end LINKS loop
    % DIVERGENCIES --------------------------------------------------------
    if node < num_nodes - NUMBER_OF_SATELLITES %(an einai satellite)
        final_result(epoch * num_nodes + num_nodes + node) = -1; 
    else
        final_result(epoch * num_nodes + num_nodes + node) = 1;
    end
    % BUFFERS -------------------------------------------------------------
    if epoch == 1
        final_result(epoch * num_nodes + num_nodes + num_nodes + node) = -1;
    elseif epoch >= 1
        final_result(epoch * num_nodes + num_nodes + num_nodes + node) = -1; % s(t)
        final_result((epoch-1) * num_nodes + num_nodes + num_nodes + node) = 1; % s(t-1)
    end
    out = final_result;
end