function out = create_R_capacity(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, xij_ca, Aeqs_cell_array, distances_ca,LINK_CAPACITY)
% Creates A matrix. (Remember: Ax<=b)
% Adding variable constraints relative to the distances between nodes.
% Link capacities should have the form of x1/R1 + ... + xn/Rn <= 1 
% Where x1 and xn are the *OUTFLOWS* of a node, and Ri the link capacities
% that are calculated as a decreasing function with distance as output
% (i.e. R(d) = -a*d + C, C <- max capacity, and 0 < a < 1)
% e.g. R(d) = -0.5*log(d) + 20, we use log(d) because distances have large
% values.
   total_epochs = length(xij_ca);
   num_nodes = NUMBER_OF_SATELLITES + NUMBER_OF_STATIONS;
   
   % Get number of columns of Aeq per epoch:
   lenAeqs = null(1,1);
   for i = 1:total_epochs
       lenAeqs = [lenAeqs, length(Aeqs_cell_array{i}(1,:))];
   end
   
   % Get all links per epoch:
   mat_xij = {};
   for epoch = 1:total_epochs
       xij = xij_ca{epoch};
       mat_xij{epoch} = zeros(length(xij),2);
       unwanted = null(1,1); 
       for i = 1:length(xij)
           mat_xij{epoch}(i,1) = xij{i}(2);
           mat_xij{epoch}(i,2) = xij{i}(3);
           if xij{i}(1) == -1
               unwanted = [unwanted, i];
           end
       end

       temp = [mat_xij{epoch}(:,2), mat_xij{epoch}(:,1)];
       
       % Add bidirectional links:
       if ~isempty(unwanted)
           temp(unwanted,:) = [];
       end
       mat_xij{epoch} = [mat_xij{epoch};...  
                         temp]; 
   end
   % Construct A matrix:
   A = null(1,sum(lenAeqs));
   for epoch = 1:total_epochs
%         disp("EPOCH ______ "+ string(epoch))

       % Epoch's characteristics
       links_parents = mat_xij{epoch}; % matrix
       distances = distances_ca{epoch};
              
       Atemp = zeros(num_nodes,lenAeqs(epoch)); % (!)
       for node = 1:num_nodes
%            disp("NODE: "+string(node)); % DEBUG
                    
           links_in_use = links_parents(links_parents(:,1) == node,:); % select the links that the given node is a parent of.        
           which_links = find(links_parents(:,1) == node); % which links are in use
%            disp("links_parents")% debug
%            disp(links_parents)% debug
%            disp("links_in_use") % debug
%            disp(links_in_use)% debug
%            disp("which_links")% debug
%            disp(which_links)% debug
           
           if ~isempty(links_in_use)
               for i = 1:length(which_links)
                   Atemp(node, which_links(i)) = (LINK_CAPACITY - 0.5 * log(distances(links_in_use(i,1),links_in_use(i,2))))^(-1); % R(d)^(-1)
               end % end for
           end % end if
       end % end for each node  
       
       % Add previous and future variable zero coefs
%        disp("lenAeqs:") % DEBUG
%        disp(lenAeqs) % DEBUG
       if epoch == 1
           right_pad = zeros(num_nodes,sum(lenAeqs(epoch+1:length(lenAeqs))));
%            disp("right_pad size:")
%            disp(size(right_pad))
%            disp("Atemp size:")
%            disp(size(Atemp))
           Atemp = [Atemp, right_pad];
       elseif epoch < total_epochs
           left_pad = zeros(num_nodes, sum(lenAeqs(1:(epoch-1))));
           right_pad = zeros(num_nodes,sum(lenAeqs((epoch+1):length(lenAeqs))));
%            disp("right_pad size:")
%            disp(size(right_pad))
%            disp("left_pad size:")
%            disp(size(left_pad))
%            disp("Atemp size:")
%            disp(size(Atemp))
           Atemp = [left_pad, Atemp, right_pad];
       elseif epoch == total_epochs
           left_pad = zeros(num_nodes, sum(lenAeqs(1:(epoch-1))));
%            disp("left_pad size:")
%            disp(size(left_pad))
%            disp("Atemp size:")
%            disp(size(Atemp))
           Atemp = [left_pad, Atemp];
       end

       % Row bind every A row:
        A = [A ; 
            Atemp];
   end % end for each epoch
   out = A;
end % end function