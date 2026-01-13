function boolValue = parseBoolean(str)
    str = lower(str);  % Convert to lowercase for consistency
    if strcmp(str, 'true') || strcmp(str, '1')
        boolValue = true;
    elseif strcmp(str, 'false') || strcmp(str, '0')
        boolValue = false;
    else
        error('Invalid boolean string: %s', str);
    end
end