function output = gk(x,y,xk,yk)
% Solution for the biharmonic
% Parameters are coordinates (x,y) and point of impulse (xk,yk)
output = (((x - xk).^2 + (y - yk).^2) .* 1/2 .* log((x - xk).^2 + (y - yk).^2));
output(isnan(output)) = 0; % accounts for variable type erros
end

