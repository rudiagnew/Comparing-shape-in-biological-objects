function [u,v] = TPS_NO_LINE_ELEMENTS(xi,yi,Xi,Yi,fineness,img)
% Creates a TPS interpolation (with no landmark crosses) given landmark data 
% ----GENERATING MATRIX EQAUTIONS----
m = length(xi); % number of points mapped
A = zeros(m);

for i = 1:m 
    for j = 1:m
        A(i,j) = gk(xi(i,1),yi(i,1),xi(j,1),yi(j,1));
    end
end

B = ones(m,3);

for i = 1:m
    B(i,2) = xi(i,1);
    B(i,3) = yi(i,1);
end

C = [A,B;transpose(B),zeros(3)];

F = linsolve(C,[Xi(1:m,1);0;0;0]); 
H = linsolve(C,[Yi(1:m,1);0;0;0]);

[coord1,coord2] = meshgrid(1:fineness:length(img),1:fineness:height(img));

% ----CALCULATING MAPPINGS----
u = 0;
for k = 1:m
    u = u + F(k)*gk(coord1,coord2,xi(k),yi(k));
end

u = u + F(m+1) + F(m+2)*coord1 + F(m+3)*coord2;

v = 0;
for k = 1:m
    v = v + H(k)*gk(coord1,coord2,xi(k),yi(k));
end

v = v + H(m+1) + H(m+2)*coord1 + H(m+3)*coord2;

end

