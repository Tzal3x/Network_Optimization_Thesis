function out = create_flow_info(LINKS,n)
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
    out = xij;
end