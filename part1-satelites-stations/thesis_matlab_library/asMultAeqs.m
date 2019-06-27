function out = asMultAeqs(Aeqs_cell_array)
% Assembles the Aeqs_cell_array. 
% Each element of Aeqs_cell_array corresponds to a Aeq matrix containingf
% the Aeq matrix of epoch i. Each Aeq contains the current buffer but not
% the previous.
% The final result is the generic Aeq used by the optimizer.
    disp("[~Report:] Assembling Aeq matrixes...")
    total_epochs = length(Aeqs_cell_array);
    ncols = null(1,1); % number of columns of Aeq at epoch i
    for i = 1:total_epochs
       ncols = [ncols, length(Aeqs_cell_array{i}(1,:))];
    end
    Aeq = null(1,1);
    n = length(Aeqs_cell_array{1}(:,1));    
    for i = 1:total_epochs
        % step 1
        Aeq_at_i = Aeqs_cell_array{i};
        % step 2
        right_bind = zeros(n, sum(ncols((i+1):length(ncols))));
        Aeq_at_i = [Aeq_at_i, right_bind];
        % step 3
        if i>1
           left_bind = zeros(n, sum(ncols(1:(i-1))));
           left_bind(:,(length(left_bind(1,:))-n+1):length(left_bind(1,:))) = diag(ones(1,n)); % previous buffers;
           Aeq_at_i = [left_bind, Aeq_at_i];
        end
        % step 4
        Aeq = [Aeq ; Aeq_at_i];
    end
    disp("[~Report:] Done!")
    out = Aeq;
end