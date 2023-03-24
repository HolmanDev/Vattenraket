function [a_vec, v_vec, s_vec, m_flow_vec] = FlightIntegral(N, dt, V_air_0, m_rocket, m_fuel, p_0, p_air_0, density_water, A_nozzle, C_discharge, angle)
    g = 9.82; % [m/s^2]
    a_vec(:, 1) = [0; 0];
    v_vec(:, 1) = [0; 0];
    s_vec(:, 1) = [0; 0];
    m_0 = m_rocket + m_fuel;
    m_e = 0;
    local_up = [cosd(angle); sind(angle)]; % Borde "pitcha" över tid
    m_flow_vec = zeros(1, N);
    p_air = p_air_0;
    v_e_air = [0; 0];
    m_e_air = 0;
    V_air = V_air_0;
    T = 10 + 273.15;
    R_spec = 350;
    p_atm = 103000;
    m_flow_air = 0;
    n = 1.135;
    r = (2 / (n + 1)) ^ (n/(n-1));
    for i=1:N
        if m_e < m_fuel
            V_air = V_air_0 + m_e / density_water;
            p_air = p_air_0 * V_air_0/V_air;
            delta_p = p_air - p_0; % p_air borde minska
        end
        
        % Exhaust-hastighet. OBS: Relativ till raketen!  
        beta = 20.5/88.1; % flaskhalsdiameter / flaskdiameter
        v_e = C_discharge * sqrt(2*delta_p / (density_water * (1 - beta^4))) * -local_up;
        % Borde inte flaskans acceleration driva ut vätskan snabbare?

        % Massflöde       
        m_flow = norm(v_e) * density_water * A_nozzle;
        m_e = m_e + m_flow*dt;
        is_empty = m_e >= m_fuel;
        density_air = p_air / (R_spec * T);
        if is_empty
            m_e = m_fuel;
            m_flow = 0;
            
            % https://en.wikipedia.org/wiki/Density_of_air
            % https://www.engineersedge.com/fluid_flow/convergent_nozzle_flow_velocity_14032.htm
            if p_air > p_atm
                v_e_air = sqrt(2*n/(n-1) * (p_air-p_atm)/density_air * (1-2/(n+1))) * -local_up; %sqrt(2*(p_air-p_atm)/density_air) * -local_up; % Lite osäker på v_e_air...
                m_flow_air = norm(v_e_air) * density_air * A_nozzle;
                m_e_air = m_e_air + m_flow_air * dt;
                p_air = p_air_0 * V_air_0/V_air * (1 - m_e_air/air_mass_at_cutoff);
                if p_air < p_atm
                    v_e_air = 0;
                    m_flow_air = 0;
                end
            else
                v_e_air = 0;
                m_flow_air = 0;
            end
        else
            air_mass_at_cutoff = V_air * density_air;
        end
        m_flow_vec(i) = m_flow;

        % Krafter
        F_prop_water = -m_flow*v_e;
        F_prop_air = -m_flow_air*v_e_air; % + delta_p*A_nozzle*(v_e_air/norm(v_e_air));
        F_g = [0; -g] * (m_0 - m_e);
        drag_coeff = 0.4;
        F_air = -drag_coeff * 1.225 * 0.5 * 0.04405^2 * pi * norm(v_vec(:, i)) * v_vec(:, i); % Ändra sen

        % Euler framåt
        a = (F_prop_water + F_prop_air + F_g + F_air) ./ (m_0 - m_e);
        a_vec(:, i+1) = a;
        v_vec(:, i+1) = v_vec(:, i) + dt*a;
        s_vec(:, i+1) = s_vec(:, i) + v_vec(:, i)*dt + 0.5*a_vec(:, i)*dt*dt;

        if(s_vec(2, i+1) < 0)
            break
        end
        
        local_up = v_vec(:, i)/norm(v_vec(:, i+1));
    end
end