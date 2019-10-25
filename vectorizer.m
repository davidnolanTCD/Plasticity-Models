function [vec] = vectorizer(tens)
vec = zeros(6,1);
for i = 1:3
    vec(i) = tens(i,i);
end

vec(4) = tens(1,2); vec(5) = tens(1,3); vec(6) = tens(2,3);


end