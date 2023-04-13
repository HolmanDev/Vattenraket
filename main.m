% Konstanter
g = 9.82; % [m/s^2]
density_water = 1000; % [kg/m^3]
kelvin_celsius_offset = 273.15;
p_atm = 102400; % [N/m^2]
V_bottle = 0.00157; % [m^3]

figure('Position', [100 100 900 600]);

n_tests = 4;
volumes = linspace(0.0004, 0.0010, n_tests); % 0.8L verkar vara bäst
for v_i = 1:n_tests  
    % Variabler
    p_air_0 = 800000;
    V_water = volumes(v_i); % [m^3]
    V_air_0 = V_bottle - V_water;
    angle = 45; % Uppskjutningsvinkel [deg]
    m_rocket = 0.2; % Raketmassa [kg]
    m_fuel = V_water * density_water; %! Add air mass?
    T_C = -1; % [C]
    T_K = T_C + kelvin_celsius_offset; % [K]
    r_nozzle = 0.01025; % Mynningsradie [m]
    d_nozzle = 2*r_nozzle; % Mynningsdiameter [m]
    A_nozzle = r_nozzle * r_nozzle * pi; % Mynningsarea [m^2]
    C_discharge = 0.98; % "Discharge coefficient [dimensionslös]
    
    % Hastighet
    N = 90000;
    dt = 0.0001;
    % Create a new method, from the ground up...

    [a_vec, v_vec, s_vec, m_flow_vec] = FlightIntegral(N, dt, V_air_0, m_rocket, m_fuel, p_atm, p_air_0, density_water, A_nozzle, C_discharge, angle);
    N = length(s_vec) - 1;
    t_vec = 0:dt:(N*dt);
    subplot(2, 3, 1);
    hold on;
    plot(t_vec, v_vec(2, :)); % v_y over time
    title("Velocity")
    xlabel("t [s]")
    ylabel("v_y [m/s]")
    yline(0, 'HandleVisibility', 'off');
    hold off
    
    % First seconds velocity
    subplot(2, 3, 2);
    hold on;
    v_y_vec = v_vec(2, :);
    T = 0.15;
    Tdt = round(T/dt);
    plot(t_vec(1:Tdt), v_y_vec(1:Tdt)); % v_y over time
    title("Velocity", "first " + T + " seconds")
    xlabel("t [s]")
    ylabel("v_y [m/s]")
    yline(0, 'HandleVisibility', 'off');
    hold off
    
    % Position
    subplot(2, 3, 3);
    hold on
    plot(s_vec(1, :), s_vec(2, :));
    title("Position")
    xlabel("x [m]")
    ylabel("y [m]")
    yline(0, 'HandleVisibility', 'off');
    hold off
    
    subplot(2, 3, 4);
    hold on
    a_t_vec = t_vec(1:N-1);
    a_y_vec = a_vec(2, :);
    plot(a_t_vec(1:Tdt), a_y_vec(1:Tdt)); % v_y over time
    title("Acceleration", "first " + T + " seconds")
    xlabel("t [s]")
    ylabel("a_y [m/s^2]")
    yline(0, 'HandleVisibility','off');
    hold off
    
    % Mass flow
    subplot(2, 3, 5);
    hold on
    plot(t_vec(1:Tdt), m_flow_vec(1:Tdt));
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