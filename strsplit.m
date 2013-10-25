function strings = strsplit(str0)

str = str0;
for s = 1:length(str0)
    if isspace(str0(s))
        str(s) = ' ';
    end
end
str = regexprep(str, ' +', ' ');

s = 1;
while str(s) == ' '
    s = s + 1;
end

if s > length(str)
    strings = cell(0);
    return;
end

t = length(str);
while str(t) == ' '
    t = t - 1;
end

A = str(s:t);
strings = cell(sum(isspace(A)) + 1, 1);

s = 1;
l = 1;
t = s;
while t <= length(A)
    if A(s) == ' ' && A(t) == ' '
        t = t + 1;
    elseif A(s) == ' ' && A(t) ~= ' '
        s = t;
        t = t + 1;
    elseif A(s) ~= ' ' && (A(t) == ' ' || t == length(A))
        if t == length(A)
            strings{l} = A(s:end);
        else
            strings{l} = A(s:t-1);
        end
        l = l + 1;
        s = t;
        t = t + 1;
    else
        t = t + 1;
    end
end
