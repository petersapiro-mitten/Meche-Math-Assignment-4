%Rate function describing Newton's law of gravitation
%INPUTS:
%t: the current time
%V: the vector of position and velocity of the planet
% V = [x_p; y_p; dxdt_p; dydt_p]
%orbit_params: a struct describing the system parameters
% orbit_params.m_sun is the mass of the sun
% orbit_params.m_planet is the mass of the planet
% orbit_params.G is the gravitational constant
%OUTPUTS:
%dVdt: a column vector describing the time derivative of V:
% dVdt = [dxdt_p; dydt_p; d2xdt2_p; d2ydt2_p]
function dVdt = gravity_rate_func(t,V,orbit_params)

% Extract position and velocity
x = V(1);
y = V(2);
vx = V(3);
vy = V(4);

% Extract system parameters
G = orbit_params.G;
m_sun = orbit_params.m_sun;

% Distance from the planet to the sun
r = sqrt(x^2 + y^2);

% Calculate acceleration due to gravity
ax = -G*m_sun*x/r^3;
ay = -G*m_sun*y/r^3;

% Construct derivative of the state vector
dVdt = [vx;
        vy;
        ax;
        ay];

end