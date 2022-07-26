function [ res ] = CORE_add_HT( head, body, tail )

if ( (numel(head) > 0)  &&  (numel(tail) > 0) )
    res = [head, body, tail];
elseif (numel(head) > 0)
    res = [head, body];
elseif (numel(tail) > 0)
    res = [tail, body];
else
    res = body;
end

end

