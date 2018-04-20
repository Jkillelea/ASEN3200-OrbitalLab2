% Part 2 of Orbital Lab 2
clear all; clc; close all;

mu_sun = 132712e6; % km^3/s^2
AU = 1.496e8; % km

r1 = [0.5887; -0.2206; 0.0239]*AU;
r2 = [0.5027;  0.2289; 0.0436]*AU;
r3 = [0.3243;  0.4560; 0.0453]*AU;

v2 = gibbs(r1, r2, r3, mu_sun, 2);

r2
v2

[h, i, a, e, Omega, omega, theta] = rv2oe(r2, v2, mu_sun)

T = 2*pi*sqrt(a^3/mu_sun)
