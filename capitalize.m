function s = capitalize(s)
mask = s == ' ';
mask = [true mask(1:end-1)];
s = lower(s);
s(mask) = upper(s(mask));