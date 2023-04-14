clear;

% Konstanter
g = 9.82; % [m/s^2]
density_w = 1000; % [kg/m^3] Saltvatten!
kelvin_celsius_offset = 273.15;
V_bottle = 0.00157; % [m^3]
bottle_height = 0.34; % [m]
RH = 0.85;
R_spec_water_vapour = 461.4;
R_spec_dry_air = 286.9;
R_spec_ambient = R_spec_dry_air + (R_spec_water_vapour-R_spec_dry_air);
n = 1.135;
adiabatic_index_air = 1.4;
adiabatic_index_w = 1.333;
r_nozzle = 0.01025; % Mynningsradie [m]
d_nozzle = 2*r_nozzle; % Mynningsdiameter [m]
A_nozzle = r_nozzle * r_nozzle * pi; % Mynningsarea [m^2]
r_bottle = 0.04405; % [m]
A_cross_section = r_bottle^2 * pi; % [m^2]

figure('Position', [100 100 900 600]);

n_tests = 4;
volumes = linspace(0.0004, 0.0010, n_tests); % 0.8L verkar vara bäst
for v_i = 1:n_tests  
    % Variabler
    p_air_0 = 800000;
    V_water = volumes(v_i); % [m^3]
    V_air_0 = V_bottle - V_water;
    angle = 45; % Uppskjutningsvinkel [deg]
    m_body = 0.2; % Raketmassa [kg]
    m_fuel = V_water * density_w; %! Add air mass?
    T_C = -1; % [C]
    T_K = T_C + kelvin_celsius_offset; % [K]
    C_discharge = 0.98; % "Discharge coefficient [dimensionslös]
    C_drag = 0.345;
    wind = [0; 0]; % m/s
    p_atm = 102400; % [N/m^2]
    % https://www.omnicalculator.com/physics/air-density
    density_amb_air = 1.30844; %p_atm / (R_spec_ambient * T_K);

    disp(density_amb_air);

    tspan = [0 9];
    y0 = [0; 0; V_air_0; 0; 0]; % [m_e_water_0; m_e_air_0; V_air_0; v_x_0; v_y_0; height]
    ODEOptions = odeset('RelTol', 1e-8);
    [t_y, y] = ode15s(@(t, y) ODESystem(t, y, p_air_0, R_spec_water_vapour, T_K, V_air_0, adiabatic_index_air, adiabatic_index_w, m_body, m_fuel, density_w, density_amb_air, p_atm, A_nozzle, A_cross_section, g, bottle_height, angle, wind, C_discharge, C_drag, n), tspan, y0, ODEOptions);
    t_total = t_y;
    
    velocity = [y(:, 4)'; y(:, 5)'];
    velocity_magnitudes = sqrt(velocity(1, :).^2 + velocity(2, :).^2);
    [t_x_pos, x_pos] = ode15s(@(t, x_pos) interp1(t_total, velocity(1, :), t), tspan, 0, ODEOptions);
    [t_y_pos, y_pos] = ode15s(@(t, y_pos) interp1(t_total, velocity(2, :), t), tspan, 0, ODEOptions);
    y_pos_sampled = interp1(t_y_pos, y_pos, t_x_pos);
    position = [x_pos, y_pos_sampled];
    acceleration = [diff(velocity(1, :) ); diff(velocity(2, :))] ./ diff(t_total');
    acceleration_magnitudes = sqrt(acceleration(1, :).^2 + acceleration(2, :).^2);

    subplot(2, 3, 1);
    hold on;
    plot(t_total, velocity(2, :));
    title("Y-Velocity")
    xlabel("t [s]")
    ylabel("v [m/s]")
    yline(0, 'HandleVisibility', 'off');
    hold off
    
    % First seconds velocity
    subplot(2, 3, 2);
    hold on;
    T = 0.15;
    dt = 0.001;
    Tdt = round(T/dt);
    plot(0:dt:T, interp1(t_total, velocity_magnitudes, 0:dt:T)); % v_y over time
    title("Velocity", "first " + T + " seconds")
    xlabel("t [s]")
    ylabel("v [m/s]")
    yline(0, 'HandleVisibility', 'off');
    hold off
    
    % Position
    subplot(2, 3, 3);
    hold on
    plot(position(:, 1), position(:, 2));
    title("Position")
    xlabel("x [m]")
    ylabel("y [m]")
    yline(0, 'HandleVisibility', 'off');
    hold off
    
    subplot(2, 3, 4);
    hold on
    plot(0:dt:T, interp1(t_total(1:(end-1)), acceleration_magnitudes, 0:dt:T)); % v_y over time
    title("Acceleration", "first " + T + " seconds")
    xlabel("t [s]")
    ylabel("a [m/s^2]")
    yline(0, 'HandleVisibility','off');
    hold off
    
    % Mass flow
    subplot(2, 3, 5);
    hold on
    m_flow = diff(y(:, 1) + y(:, 2)) ./ diff(t_total);
    plot(0:dt:T, interp1(t_total(1:(end-1)), m_flow, 0:dt:T));
    title("Mass flow", "first " + T + " seconds")
    xlabel("t [s]")
    ylabel("dm/dt [kg/s]")
    hold off
end

% Visa legends
subplot(2, 3, 6);
plot(nan(n_tests));
set(gca, 'Visible', 'off');
lgd_values = round(volumes'*1000, 3);
lgd = legend(strtrim(cellstr(num2str(lgd_values))'), 'Location', 'east');
title(lgd, "Volumes [L]");