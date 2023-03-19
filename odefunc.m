function dsdt = odefunc(t, vel, dt, V_air_0, m_rocket, m_fuel, p_0, p_air, density_water, A_nozzle, C_discharge)
    % vel = y
    g = 9.82; % [m/s^2]
    m_0 = m_rocket + m_fuel;
    m_e = 0;
    local_up = 1/sqrt(2)*[1; 1]; % Borde "pitcha" över tid

    V_air = V_air_0 + m_e / density_water;
    p = p_air * V_air_0/V_air;
    delta_p = p - p_0; % p_air borde minska
    
    % Exhaust-hastighet. OBS: Relativ till raketen!  
    beta = 20.5/88.1; % flaskhalsdiameter / flaskdiameter
    v_e = C_discharge * sqrt(2*delta_p / (density_water * (1 - beta^4))) * -local_up;

    % Massflöde       
    m_flow = norm(v_e) * density_water * A_nozzle;
    m_e = m_e + m_flow*dt;
    is_empty = m_e >= m_fuel;
    if is_empty
        m_e = m_fuel;
        m_flow = 0;
    end

    % Krafter
    F_prop = -m_flow*v_e;
    F_g = [0; -g] * (m_0 - m_e);
    drag_coeff = 0.4;
    F_air = -drag_coeff * 1.225 * 0.5 * 0.04405^2 * pi * norm(vel) * vel; % Ändra sen

    % Euler framåt
    a = (F_prop + F_g + F_air) ./ (m_0 - m_e);
    vel = vel + dt*a; % Velocity (saves)
end