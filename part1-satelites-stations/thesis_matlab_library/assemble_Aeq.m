function out = assemble_Aeq(aeq_rows)
% Row binding the vectors of aeq_rows vector.
    disp("[~Report:] Assembling Aeq matrix...")
    Aeq = [];
    for i = 1:length(aeq_rows)
        Aeq = [Aeq; aeq_rows{i}];
    end
    disp("[Report:] Done!")
    out = Aeq;
end