function index = find_index(c)
l=length(c);
s = sum(c);
index = 1;
temp = rand(1);
frac = c(1) / s;
while frac < temp && index <= l
    index = index + 1;
    frac = frac + c(index) / s;
end
end