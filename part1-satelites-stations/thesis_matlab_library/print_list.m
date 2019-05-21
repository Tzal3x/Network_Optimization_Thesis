function out = print_list(lista)
    temp = '';
    for i = 1:length(lista)
        temp = [temp  ' '  convertStringsToChars(string(lista(i)))];
    end
    out = temp;
end
