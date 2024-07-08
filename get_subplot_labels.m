function charlbl = get_subplot_labels(strchain,nIDs)
    
    alphabet = (strchain).';
    chars = num2cell(alphabet(1:nIDs));
    chars = chars.';
    charlbl = strcat(chars,')');

end