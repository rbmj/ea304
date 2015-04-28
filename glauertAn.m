% Note, we will assume symmetric loading.
%
% function[] = glauertAn(alpha, a0, alpha0, c, b, {N,y})
%
% Return Values:
%   - A:  N fourier coefficients of the lift distribution.  These are all
%         odd numbered coefficients due to symmetric loading.
%
% Arguments:
%   - alpha: wing angle of attack
%   - a0:  function handle that returns the section lift curve slope
%               (per radian) as a function of spanwise position, y.
%   - alpha0: function handle that returns the zero lift angle of attack
%               (in radians) as a function of spanwise position, y.
%   - c: function handle that returns the chord as a function of spanwise
%               position, y.
%   - b: the span of the wing
%   - N: the number of stations to use
%   - y:  The spanwise locations of the stations
%
% Notes:
%   - All function handles must support vectorization
%   - Both geometric and aerodynamic twist are controlled through alpha0
%   - Units for length are arbitrary but must be consistent
%   - Give either N or y - not both. If given as N, stations will be
%     linearly spaced in theta-space 

function[A, cn] = glauertAn(alpha, a0, alpha0, c, b, N)
    if isscalar(N)
        % The tip of the wing having zero circulation causes some issues,
        % so ignore the last 1% of the span.
        theta = linspace(pi/2, pi*0.99, N);
        y = cos(theta)*(b/-2);
    else
        y = N;
        N = numel(y);
        theta = acos(-2*y/b);
    end
    % Calculate planform parameter mu for all stations
    mu = a0(y) .* c(y) ./ (4*b);
    % Calculate I/B matrices
    I = zeros(N);
    B = zeros(N, 1);
    for i = 1:N
        B(i) = (alpha - alpha0(y(i)))*mu(i)*sin(theta(i));
        for j = 1:N
            n = 2*j - 1; % Only need odd fourier coefficients
            I(i, j) = (n*mu(i) + sin(theta(i))) * sin(theta(i)*n);
        end
    end
    A = (I\B)';
    cn = rcond(I);
end
