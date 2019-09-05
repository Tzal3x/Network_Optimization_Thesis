function out = create_R_capacity(NUMBER_OF_SATELLITES, NUMBER_OF_STATIONS, xij_ca, Aeqs_cell_array, distances_ca, LINK_CAPACITY)
% Creates A matrix v.2.1 (Containing interference factor)

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
       links = xij_ca{epoch};
       mat_xij{epoch} = zeros(length(links),2);
       unwanted = null(1,1); 
       for i = 1:length(links)
           mat_xij{epoch}(i,1) = links{i}(2);
           mat_xij{epoch}(i,2) = links{i}(3);
           if links{i}(1) == -1
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
   A = null(1,1);
   for epoch = 1:total_epochs
       % Epoch's characteristics
       links = mat_xij{epoch}; % matrix [parent 1, parent2 ; parent 1 parent2; ...]
       if isempty(links)
          continue;
       end          
           
       distances = distances_ca{epoch};
       num_links = length(links(:,1));
       
       % Construct R(d) vector for every link:
       Atemp = null(1); % (!)
       for link1 = 1:num_links % for every link
           par_1 = links(link1,1);
           par_2 = links(link1,2);
           temp_Rd_vector = zeros(1, num_links); % each row corresponds to a part of A matrix row.
           for link2 = 1:num_links % iterate every other link and pick it if...
              temp_par_1 = links(link2,1);
              temp_par_2 = links(link2,2);
              if (temp_par_1 == par_1) || (temp_par_1 == par_2) || (temp_par_2 == par_1) || (temp_par_2 == par_2)
                 temp_Rd_vector(link2) = (LINK_CAPACITY - 0.5 * log( distances(temp_par_1, temp_par_2) ))^(-1); % R(d)^(-1)
              end
           end
           Atemp = [Atemp; temp_Rd_vector];
       end % end for each link 
       disp("DEBUG PAUSE RD")
       disp(Atemp);pause;% DEBUG
       
       % Add previous and future variable zero coefs
       disp("lenAeqs:") % DEBUG
       disp(lenAeqs) % DEBUG
       if epoch == 1
           right_pad = zeros(num_links, num_nodes*2 + sum(lenAeqs(epoch+1:length(lenAeqs))));
%            disp("right_pad size:")
%            disp(size(right_pad))
%            disp("Atemp size:")
%            disp(size(Atemp))
%            pause;
           Atemp = [Atemp, right_pad];
       elseif epoch < total_epochs
           left_pad = zeros(num_links, sum(lenAeqs(1:(epoch-1))));
           right_pad = zeros(num_links, num_nodes*2 + sum(lenAeqs((epoch+1):length(lenAeqs))));
%            disp("right_pad size:")
%            disp(size(right_pad))
%            disp("left_pad size:")
%            disp(size(left_pad))
%            disp("Atemp size:")
%            disp(size(Atemp))
           Atemp = [left_pad, Atemp, right_pad];
       elseif epoch == total_epochs
           left_pad = zeros(num_links, sum(lenAeqs(1:(epoch-1))));
           right_pad = zeros(num_links, num_nodes*2);
%            disp("left_pad size:")
%            disp(size(left_pad))
%            disp("Atemp size:")
%            disp(size(Atemp))
           Atemp = [left_pad, Atemp, right_pad];
       end

       % Row bind every A row:
        A = [A ; 
            Atemp];
   end % end for each epoch
   disp(A)
   out = A;
end % end function