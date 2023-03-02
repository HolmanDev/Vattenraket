volumes = [0.0005, 0.00075, 0.001, 0.00125, 0.0015];

for v_i = 1:5
    % Konstanter
    g = 9.82; % [m/s^2]
    density_water = 997.13; % [kg/m^3]
    kelvin_celsius_offset = 273.15;
    p_0 = 101300; % [N/m^2]
    V_bottle = 0.0016; % [m^3]
    
    % Variabler
    p_air_0 = 800000;
    V_water = volumes(v_i); % [m^3]
    V_air_0 = V_bottle - V_water;
    angle = 45; % Uppskjutningsvinkel [deg]
    m_rocket = 0.2; % Raketmassa [kg]
    m_fuel = V_water * density_water;
    T_C = 20; % [C]
    T_K = T_C + kelvin_celsius_offset; % [K]
    r_nozzle = 0.01025; % Mynningsradie [m]
    d_nozzle = 2*r_nozzle; % Mynningsdiameter [m]
    A_nozzle = r_nozzle * r_nozzle * pi; % Mynningsarea [m^2]
    C_discharge = 0.5; % "Discharge coefficient [dimensionslös]
    
    % Hastighet
    N = 90000;
    dt = 0.0001;
    [v_vec, m_flow_vec] = Velocity(N, dt, V_air_0, m_rocket, m_fuel, g, p_0, p_air_0, density_water, A_nozzle, C_discharge);
    t_vec = 0:dt:(N*dt);
    figure(1)
    hold on;
    plot(t_vec, v_vec(2, :)); % v_y over time
    title("Velocity")
    xlabel("t [s]")
    ylabel("v_y [m/s]")
    yline(0, 'HandleVisibility', 'off');
    hold off
    
    % First seconds velocity
    figure(2)
    hold on;
    v_y_vec = v_vec(2, :);
    T = 0.5;
    Tdt = round(T/dt);
    plot(t_vec(1:Tdt), v_y_vec(1:Tdt)); % v_y over time
    title("Velocity", "first " + T + " seconds")
    xlabel("t [s]")
    ylabel("v_y [m/s]")
    yline(0, 'HandleVisibility', 'off');
    hold off
    
    % Position
    s_vec = [0;0];
    for i=1:N
        s_vec(:, i+1) = s_vec(:, i) + v_vec(:, i)*dt;
    end
    figure(3)
    hold on
    plot(s_vec(1, :), s_vec(2, :));
    title("Position")
    xlabel("x [m]")
    ylabel("y [m]")
    yline(0, 'HandleVisibility', 'off');
    hold off
    
    % Acceleration
    a_vec = zeros(2, N-1);
    for i=1:(N-1)
        a_vec(:, i) = (v_y_vec(:, i+1) - v_y_vec(:, i))/dt;
    end
    figure(4)
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
    figure(5)
    hold on
    plot(t_vec(1:Tdt), m_flow_vec(1:Tdt));
    title("Mass flow", "first " + T + " seconds")
    xlabel("t [s]")
    ylabel("dm/dt [kg/s]")
    hold off
end

for i=1:5
    figure(i)
    legend(strtrim(cellstr(num2str(volumes'))')) % Color mismatch?
end