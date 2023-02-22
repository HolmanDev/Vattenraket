% Konstanter
g = 9.82; % [m/s^2]
density_water = 997.13; % [kg/m^3]
kelvin_celsius_offset = 273.15;
p_0 = 101300; % [N/m^2]

% Variabler
angle = 45; % Uppskjutningsvinkel [deg]
p = 100; % Lufttryck i raketen [N/m^2]
T_C = 20; % [C]
T_K = T_C + kelvin_celsius_offset; % [K]
m_raket = 0.1; % Raketmassa [kg]
r_nozzle = 0.01025; % Mynningsradie [m]
d_nozzle = 2*r_nozzle; % Mynningsdiameter [m]
A_nozzle = r_nozzle * r_nozzle * pi; % Mynningsarea [m^2]
cd = 1; % "Discharge coefficient [dimensionslös]
p_air = 300000; %! MÅSTE RÄKNAS UT [N/m^2]

% Räkna ut massflöde genom mynningen
m_flow = A*Cd*sqrt(2*density_water*(p_air-p_0));

% Uträkning
speed0 = 10; % räkna ut sen
v0 = speed0 * [cosd(angle) sind(angle)]; 
% v_y*t - (g*t^2)/2 = 0
% t^2 - (2/g)*v_y*t = 0
% t = (2/g)*v_y
t = v0(2) * 2/g;
d = v0(1) * t; % Nu räknar jag bara första kicken
tVec = linspace(0, t, 20);
a = [0 -g];
P = get_position(v0, a, tVec);
plot(P(1, :), P(2, :));
disp("Max höjd: " + max(P(2, :)) + "m")

% d = Integral(v_x * dx)

function P = get_position(v, a, tVec)
    n = size(tVec, 2);
    P = zeros(2, n);
    P(1, :) = v(1).*tVec + a(1)*0.5.*tVec.*tVec;
    P(2, :) = v(2).*tVec + a(2)*0.5.*tVec.*tVec;
end
