% This function, using the equations from Wikipedia, returnes two closures
% which can compute the airfoil shape and mean camber line given the
% NACA airfoil number.
% 
% As an example, to graph the NACA 2412 airfoil:
% >> xbar = linspace(0, 1);
% >> airfoil_coords = naca_coords(2412);
% >> [upper, lower] = airfoil_coords(xbar);
% >> plot(upper(1,:), upper(2,:), lower(1,:), lower(2,:));
% >> axis([0 1 -0.5 0.5]);
%
% The MCL can be graphed in a similar manner, with e.g.:
% >> [~, mcl] = naca_coords(2412);

function[coord_fn, mcl_fn] = naca_coords(airfoil_number)
    [max_camber, camber_loc, percent_thick] = naca4digit(airfoil_number);
    function[yc] = mcl(xbar)
        if xbar < camber_loc
            yc = max_camber.*xbar.*(2*camber_loc - xbar)/(camber_loc^2);
        else
            yc = max_camber.*(1-xbar).*(1 + xbar - 2*camber_loc)/(1-camber_loc)^2;
        end
    end
    function[yc_prime] = mcl_slope(xbar)
        if xbar < camber_loc
            yc_prime = 2*max_camber*(camber_loc - xbar)/camber_loc^2;
        else
            yc_prime = 2*max_camber*(camber_loc - xbar)/(1 - camber_loc)^2;
        end
    end
    function[yt] = thickness(xbar)
        yt = 5*percent_thick*(0.2969*sqrt(xbar) - 0.1260*xbar ...
            - 0.3516*xbar.^2 + 0.2843*xbar.^3 - 0.1015*xbar.^4);
    end
    function [upper, lower] = get_coords(xbar)
        theta = atan(mcl_slope(xbar));
        yt = thickness(xbar);
        yc = mcl(xbar);
        upper = [(xbar - yt.*sin(theta)); (yc + yt.*cos(theta))];
        lower = [(xbar + yt.*sin(theta)); (yc - yt.*cos(theta))];
    end
    coord_fn = @get_coords;
    mcl_fn = @mcl;
end

function[max_camber, camber_loc, thickness] = naca4digit(airfoil_number)
    max_camber = round(airfoil_number/1000)/100;
    camber_loc = round(mod(airfoil_number, 1000)/100)/10;
    thickness = mod(airfoil_number, 100)/100;
end