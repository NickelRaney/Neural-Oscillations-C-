function a = real2int(b)
f = floor(b);
if rand(1) < b - f
    a = f + 1;
else
    a = f;
end
end