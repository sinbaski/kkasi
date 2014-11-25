function central = central_moments(noncentral)
central = NaN(1, length(noncentral));

for a = 1:length(noncentral)
    s = 0;
    for l = 1:a
        s = s + nchoosek(a, l) * (-noncentral(1))^(a-l) * noncentral(l);
    end
    s = s + (-noncentral(1))^a;
    central(a) = s;
end

