function[] = project2()
    close all;
    %% Setup
    % Flight parameters, from the project description
    V_inf = 322.7; %fps
    rho_inf = atmos_rho(8000);
    alpha = 0;
    % Wing parameters, from the project description
    % Note:  There are some simplifications here.  We assume a linearly
    % tapered wing with 0 sweep at the quarter chord.  This is somewhat
    % different than the actual wing geometry.
    b = 24 + (4/12); %feet
    root_twist = 1.5; % degrees
    washout = 2.5; %degrees
    taper = 0.5246;
    c_r = 5 + (1/12);
    % These numbers are from A&VD for NACA 64A210, NACA 64A212
    a0_t = 1.824; % times pi, per radian
    a0_r = 1.954;
    alpha0_t = -2; % degrees
    alpha0_r = -1.5; % degrees
    C_m = -0.04; % roughly constant along span
    % Functions to determine the wing properties at any spanwise location
    a0 = linear_a0(a0_r, a0_t, b);
    alpha0 = linear_alpha0(alpha0_r, alpha0_t, root_twist, washout, b);
    c = linear_c(c_r, taper, b);
    % Dependent Wing Parameters
    cbar = c_r*(1+taper)/2;
    S = b*cbar;
    AR = b*b/S;
    MAC = (2/3)*c_r*(1+taper+taper^2)/(1+taper);
    %% Convergence Study
    function[a1, cn] = conv_study(num)
        [coeffs, cn] = glauertAn(alpha, a0, alpha0, c, b, num);
        a1 = coeffs(1);
    end
    N = 2:100;
    [A, cn] = arrayfun(@conv_study, N);
    plot(N, A);
    title('A_1 vs Number of Stations');
    figure();
    plot(N, cn);
    title('Reciprocal Condition Number vs. Number of Stations');
    N = find(abs(A - A(end)) < 1e-5);
    N = N(1);
    display(N);
    %% Study of initial configuration
    % Compute Lift Distribution
    n = (1:N)*2 - 1;
    A = glauertAn(alpha, a0, alpha0, c, b, N);
    % Get some sample points in y-space for graphing, and convert to
    % theta-space.
    y_pos = linspace(0, (b/2)*0.995, 64);
    y = cat(2, fliplr(-y_pos(2:end)), y_pos);
    theta = acos(-2*y/b);
    % Compute coefficients
    [C_L, C_Di, e] = getCoeffs(AR, A, n);
    C_M = C_m*MAC/cbar; % C_m is roughly constant so this works
    % Get spanwise distributions
    alpha_i = arrayfun(@(t) sum(n .* A .* sin(n*t)) / sin(t), theta);
    gamma = arrayfun(@(t) 2*b*V_inf*sum(A.*sin(n*t)), theta);
    C_l = gamma .* (2./(V_inf*c(y)));
    dL = gamma * V_inf * rho_inf; % dy
    pos_i = (y >= 0);
    i_begin = find(y == 0, 1);
    MP_spanwise = arrayfun(@(i) trapz(y(i:end), dL(i:end) .* ...
        (y(i:end) - y(i)), 2), (i_begin:numel(y)));
    MY_spanwise = arrayfun(@(i) trapz(y(i:end), dL(i:end) .* ...
        alpha_i(i:end) .* (y(i:end) - y(i)), 2), (i_begin:numel(y)));
    % Graphs & Output
    figure();
    plot(y, alpha_i*(180/pi));
    title('Induced AOA Distribution');
    figure();
    alpha_eff = alpha - alpha_i - alpha0(y);
    plot(y, alpha_eff*(180/pi));
    title('Effective AOA Distribution');
    figure();
    plot(y, gamma);
    title('Circulation Distribution');
    figure();
    plot(y, -alpha_i*V_inf);
    title('Downwash Distribution');
    figure();
    plot(y, C_l./C_L);
    title('Normalized Section Lift Coeff. Distribution');
    figure();
    plot(y, (C_l .* alpha_i)./C_Di);
    title('Normalized Section Induced Drag Coeff. Distribution');
    figure();
    plot(y(pos_i), MP_spanwise);
    title('Spanwise Pitching Moment Distribution');
    figure();
    plot(y(pos_i), MY_spanwise);
    title('Spanwise Yawing Moment Distribution');
    display(C_L);
    display(C_Di);
    display(e);
    MP_root = MP_spanwise(1);
    MY_root = MY_spanwise(1);
    display(MP_root);
    display(MY_root);
    display(C_M);
    q = rho_inf*V_inf*V_inf/2;
    % Total aerodynamic force
    F = q*S*sqrt(C_L^2 + C_Di^2);
    load_per_bolt = F/4; % assume uniform loading
    display(load_per_bolt);
    %% Existing Configuration Variable \alpha performance:
    alphav = linspace(-5, 5, 50)*pi/180;
    [clv, cdiv, ev] = arrayfun(@(a) ...
        getCoeffs(AR, glauertAn(a, a0, alpha0, c, b, N), n), alphav);
    figure();
    plot(alphav*180/pi, clv);
    title('A/C Lift Curve');
    figure();
    plot(clv, cdiv);
    title('A/C Drag Polar');
    figure();
    plot(alphav*180/pi, ev);
    title('Effect of \alpha on e');
    figure();
    plot(alphav*180/pi, clv./cdiv);
    title('Effect of \alpha on L/D');
    %% Taper effects study
    a0_fn = linear_a0(a0_r, a0_t, b);
    alpha0_fn = linear_alpha0(alpha0_r, alpha0_t, root_twist, washout, b);
    lambda_vec = [(2/3) (1/2) (1/3) (1/6)];
    c_fn = arrayfun(@(lambda) linear_c(c_r, lambda, b), lambda_vec, ...
        'UniformOutput', false);
    [~, ~, ev] = cellfun(@(c_fn_i) getCoeffs(AR, ...
        glauertAn(alpha, a0_fn, alpha0_fn, c_fn_i, b, N), n), c_fn);
    figure();
    plot(lambda_vec, ev);
    title('e vs \lambda');
    [~, i] = max(ev);
    taper = lambda_vec(i);
    c_fn = c_fn{i};
    display(taper);
    %% Washout effects study
    a0_fn = linear_a0(a0_r, a0_t, b);
    wash_vec = [1 2.5 5];
    alpha0_fn = arrayfun(@(wash) ...
        linear_alpha0(alpha0_r, alpha0_t, root_twist, wash, b), wash_vec, ...
        'UniformOutput', false);
    [~, ~, ev] = cellfun(@(alpha0_fn_i) getCoeffs(AR, ...
        glauertAn(alpha, a0_fn, alpha0_fn_i, c_fn, b, N), n), alpha0_fn);
    figure();
    plot(wash_vec, ev);
    title('e vs washout');
    [~, i] = max(ev);
    wash = wash_vec(i);
    display(wash);
    %% Combined taper/washout
    lambda = linspace(0, 1, 50);
    washout = linspace(1, 5, 10);
    values = zeros(numel(washout), numel(lambda));
    alpha0_fn = arrayfun(@(wash) linear_alpha0(alpha0_r, alpha0_t, ...
        root_twist, wash, b), washout, 'UniformOutput', false);
    c_fn = arrayfun(@(tap) linear_c(c_r, tap, b), lambda, ...
        'UniformOutput', false);
    for i = 1:numel(washout)
        alpha0_h = alpha0_fn{i};
        for j = 1:numel(lambda)
            c_h = c_fn{j};
            % It's OK not to recompute AR here as it's not actually used
            % to compute e, which is all we care about.
            [~, ~, e] = getCoeffs(AR, glauertAn(alpha, a0_fn, alpha0_h, ...
                c_h, b, N), n);
            values(i, j) = e;
        end
    end
    figure();
    surf(lambda, washout, values, 'EdgeColor', 'None', 'FaceColor', 'interp');
    view(2);
    colorbar('EastOutside');
    title('e vs. \lambda and washout');
    xlabel('\lambda');
    ylabel('washout (degrees)');
    %% Study of Proposed Design
    a0 = linear_a0(a0_r, a0_t, b);
    alpha0 = linear_alpha0(alpha0_r, alpha0_t, root_twist, wash, b);
    c = linear_c(c_r, taper, b);
    A = glauertAn(alpha, a0, alpha0, c, b, N);
    % Compute coefficients
    cbar = c_r*(1+taper)/2;
    S = b*cbar;
    AR = b*b/S;
    MAC = (2/3)*c_r*(1+taper+taper^2)/(1+taper);
    [C_L, C_Di, e] = getCoeffs(AR, A, n);
    C_M = C_m*MAC/cbar; % C_m is roughly constant so this works
    % Get spanwise distributions
    C_l = gamma .* (2./(V_inf*c(y)));
    figure();
    plot(y, C_l./C_L);
    title('Normalized Section Lift Coeff. Distribution');
    display(C_L);
    display(C_Di);
    display(C_M);
    display(e);
    display(taper);
    display(wash);
end

function[rho] = atmos_rho(h)
    [~, ~, ~, rho] = atmosisa(h*0.3048);
    rho = rho * 0.00194; % conv to slug/ft^3
end

function[C_L, C_Di, e] = getCoeffs(AR, A, n)
    C_L = pi*AR*A(1);
    delta = sum(n(2:end) .* (A(2:end)./A(1)).^2);
    e = 1/(1+delta);
    C_Di = C_L*C_L/(pi*e*AR);
end

function[a0h] = linear_a0(a0_root, a0_tip, b)
    function[x] = a0(y)
        x = pi*(a0_root + ((a0_tip-a0_root)*(y./(b/2))));
    end
    a0h = @a0;
end
function[alpha0h] = linear_alpha0(alpha0_root, alpha0_tip, geom_twist_root, washout, b)
    function[x] = alpha0(y)
        ybar = abs(y)./(b/2);
        % Account for aerodynamic twist
        alpha0_aero = alpha0_root + (alpha0_tip-alpha0_root)*ybar;
        alpha0_aero = alpha0_aero * (pi/180);
        % Account for geometric washout
        alpha0_geom = geom_twist_root - washout*ybar;
        alpha0_geom = -alpha0_geom * (pi/180);
        x = alpha0_aero + alpha0_geom;
    end
    alpha0h = @alpha0;
end
function[ch] = linear_c(c_root, taper, b)
    function[x] = c(y)
        x = c_root*(1 - (1-taper)*(abs(y)./(b/2)));
    end
    ch = @c;
end
