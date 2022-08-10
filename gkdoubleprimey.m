function output = gkdoubleprimey(x,y,xk,yk)
% Solution for the biharmonic
% Parameters are coordinates (x,y) and point of impulse (xk,yk)
output =  (2 .* ((y - yk).^2)./((x - xk).^2 + (y - yk).^2)) + (log((x - xk).^2 + (y - yk).^2));
output(isnan(output)) = 0; % accounts for variable type erros
end

