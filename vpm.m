% Perform a simple numerical calculation of the lift coefficient
% using a first order lumped vortex panel method.  This is for thin
% airfoils (i.e. thickness is neglected and a single sheet of lumped
% vortex panels is constructed along the mean camber line).  It uses
% the Aerodynamic Influence Coefficient technique.
%
% Arguments:
%   - alpha:   Angle of attack of the freestream, in radians
%   - n:       Number of vortex panels to use
%   - mcl:     Function handle, mcl(x_bar) = z_bar | (x_bar, z_bar) is on
%              the mean camber line of the airfoil.  Note that both x_bar
%              and y_bar must be chord normalized/non-dimensionalized.
%
% Return Value:
%   - c_l:     Section lift 2-D Lift Coefficient
%   - D_cp:    Two column, pressure jump coeff. vs. xbar.  This is a
%              returned in a cell in order to keep e.g. arrayfun() from
%              complaining in the common use case.
%   - c_m_ac:  Pitching moment coefficient about the aerodynamic center,
%              which we assume to be the quarter chord.

function[c_l, D_cp, c_m_ac] = vpm(alpha, n, mcl)
    % get endpoints of the vortex panels
    points = linspace(0, 1, n+1);
    points = cat(1, points, mcl(points))';

    % initialize variables
    theta = zeros(n, 1);
    vortex_coords = zeros(n, 2);
    cp_coords = zeros(n, 2);

    % get points and theta angles
    for i = 1:n
        delta = points(i + 1,:) - points(i,:);
        % place vortex at quarter chord of panel
        vortex_coords(i,:) = points(i,:) + delta*0.25;
        % place CP at three quarter chord of panel
        cp_coords(i,:) = points(i,:) + delta*0.75;
        theta(i) = atan2(delta(2), delta(1));
    end

    % generate AICs
    a = zeros(n);
    for i = 1:n
        for j = 1:n
            r_vec = (cp_coords(i,:) - vortex_coords(j,:));
            n_vec = [-cos(theta(i)), sin(theta(i))];
            a(i,j) = dot(r_vec, n_vec)/(norm(r_vec)^2);
        end
    end

    % calculate circulation and lift, using matrix math
    gamma_bar = a\(theta - alpha);
    c_l = 4*pi*sum(gamma_bar);
    
    % calculate pressure jump coefficient, clamping to zero at LE and TE
    D_cp = cat(2, linspace(1/(2*n), 1-1/(2*n), n)', 4*pi*n*gamma_bar);
    D_cp = {cat(1, [0 0], D_cp, [1 0])};
    
    % calculate c_m_ac, assuming ac falls at x = c/4 -> xbar = 0.25
    c_m_ac = -4*pi*dot(vortex_coords(:,1) - 0.25, gamma_bar);
end
