function [M] = LaplaceEquation(fx,gx,fy,gy,a,b,c,d,n,m)
syms x y

% Consider the Laplace Equation 
% Uxx+Uyy=0 for a < x < b, c < y < d
% Subject to the following boundary conditions: 
% U(a,y) = f(y), U(b,y) = g(y)
% U(x,c) = f(x), U(x,d) = g(x)
%
% According to the finite difference method a mesh (n x m) should be
% created as a representation of the region, each boundary point holding
% a known value and each internal point being the average from its for 
% surrounding points. The former situation generates a linear system of
% equations of (n-1) x (m-1) in which the coefficients give an approximate
% value for the Laplace equation in (x,y)
% Consider no value for the corners.

% Generate Mesh
M = zeros(n,m);
% Fill edges with boundary condition
for i = 1:n
    for j = 1:m
        if i == 1 && j ~= 1 && j ~= m
            M(i,j) = eval(subs(fy,y,a));
        elseif i == n && j ~= 1 && j ~= m
            M(i,j) = eval(subs(gy,y,b));
        end
        if j == 1 && i ~= 1 && i ~= n
            M(i,j) = eval(subs(fx,x,c));
        elseif j == m && i ~= 1 && i ~= n
            M(i,j) = eval(subs(gx,x,d));
        end
    end
end

% old
MO = M;

flag = 1;
% Smart iteration stop missing
while flag < 100
    % The core idea is mapping the matrix into a vector with index phi
    % this way you can treat each interior point as consecutive, which
    % easily lets you to generalize the average formula for every interior
    % point.
    %
    % Given that the first matrix position is addressed by (1,1):
    % i(phi) = fix(phi/m) + 1
    % j(phi) = mod(phi/m)
    % phi(i,j) = m*(i-1) +j
    %
    for phi = m+2:n*m
        % Unless the phi index adresses a boundary it is set equal to the
        % average of the new upper and left points and the old down and
        % right points.
        if mod(phi,m) ~= 0 && mod(phi,m) ~= 1 && fix(phi/m) + 1 ~= n
            % M(i,j) = [ M(i-1,j) + M(i,j-1) + M(i+1,j) + M(i,j+1) ] / 4
            % Once mapped and changing M for MO:
            M(fix(phi/m)+1,mod(phi,m)) = ( MO(fix(phi/m)+1+1,mod(phi,m)) + ... 
                MO(fix(phi/m)+1,mod(phi,m)+1) + M(fix(phi/m)+1-1,mod(phi,m)) + ...
                M(fix(phi/m)+1,mod(phi,m)-1) ) / 4;
        end
    end
    flag = flag +1;
    MO = M;
end

mesh(M)

return
end