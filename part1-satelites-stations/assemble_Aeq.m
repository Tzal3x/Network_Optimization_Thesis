function out = assemble_Aeq(aeq_rows)
    disp("[~Report:] Assembling Aeq matrix...")
    Aeq = [];
    for i = 1:length(aeq_rows)
        Aeq = [Aeq; aeq_rows{i}];
    end
    disp("[Report:] Done!")
    out = Aeq;
end