clear all; clc;
xhat = [1; 0; 0];
yhat = [0; 1; 0];
zhat = [0; 0; 1];

mu = 42828; % mars, km^3/s^2

r  = [-3424.7; -47.5; 1172]; % km (x, y, z)
vy = -3.333; % km/s
vx = (r(1)*vy - 11394.338)/r(2);
vz = (3948.694 + r(3)*vy)/r(2);
v  = [vx; vy; vz];

% == Keplerian Orbital Parameters == %
% angular momentum
h = cross(r, v);
% specific energy
epsilon = (norm(v)^2)/2 - mu/norm(r);
% eccentricity
e = sqrt(1 + (2*epsilon*norm(h)^2)/mu^2);
% eccentricity vector
e_vec = cross(v, h)/mu - r/norm(r);
% inclination
i = acosd(h(3)/norm(h));
% semimajor axis
a = -mu/(2*epsilon);
% line of nodes
n_vec = cross(zhat, h); % zhat x h
% Right ascension of ascending node
Omega = acosd(dot(n_vec, xhat)/norm(n_vec));
% argument of periapsis
omega = acosd(dot(n_vec, e_vec)/(norm(n_vec)*e));
% true anomaly
true_anom = acosd(dot(r, e_vec)/(norm(r)*e));
% orbital period, seconds
T = 2*pi*sqrt((a^3)/mu);
