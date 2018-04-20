function [h, i, a, e, Omega, omega, theta] = rv2oe(r, v, mu)
    xhat = [1, 0, 0];
    yhat = [0, 1, 0];
    zhat = [0, 0, 1];

    % mu = 398600; % km^3/s^2, earth
    h = cross(r, v);
    epsilon = dot(v, v)/2 - mu/norm(r);
    a = -mu/(2*epsilon);
    i = acosd(h(3)/norm(h));
    e = sqrt(1 + 2*dot(h, h)*epsilon/(mu^2));
    n = cross(zhat, h);

    % RAAN with quadrant check
    Omega = abs(acosd(dot(n, xhat)/norm(n)));
    if dot(n, yhat) < 0
        Omega = -Omega;
    end

    e_vec = cross(v, h)/mu - r/norm(r);

    % argument of periapsis with quadrant check
    omega = abs(acosd(dot(n, e_vec)/(e*norm(n))));
    if dot(e_vec, zhat) < 0
        omega = -omega;
    end

    % true anomaly with quadrant check
    theta = abs(acosd(dot(r, e_vec)/(e*norm(r))));
    if dot(r, v) < 0
        theta = -theta;
    end
end
