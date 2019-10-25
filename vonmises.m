function [out] = vonmises(s)

p1 = (s(1)-s(2))^2 + (s(2)-s(3))^2 + (s(3)-s(1))^2;
p2 = 6 * ( s(4)^2 + s(5)^2 + s(6)^2 );

out = sqrt( (p1 + p2)/2 );

end