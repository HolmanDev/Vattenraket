function dy = ODESystem(t, y, p_air_0, R_spec_air, T, V_air_0, adiabatic_index_air, adiabatic_index_w, m_body, m_fuel_0, density_w, density_amb_air, p_atm, A_nozzle, A_cross_section, g, bottle_height, angle, wind, C_discharge, C_drag, n)
    m_e_water = y(1);
    m_e_air = y(2);
    V_air = y(3);
    v_vec = [y(4); y(5)];

    % Misc
    local_up = DirectionOf(v_vec);
    if t < 0.05
        local_up = [cosd(angle); sind(angle)];
    end
    m_tot = m_body + m_fuel_0 - m_e_water; % Air isn't counted (fix!)
    v_e_water_vec = [0; 0];
    v_e_air_vec = [0; 0];
    m_flow_air = 0;

    % Properties of air inside the bottle
    density_air_0 = p_air_0 / (R_spec_air * T);
    density_air = density_air_0 * (V_air_0 / V_air)^adiabatic_index_air;
    m_air_0 = V_air_0 * density_air_0;
    m_e_air = max(0, min(m_air_0, m_e_air));
    air_left_coeff = (1 - m_e_air/m_air_0); % n ratio in ideal gas law (?)
    p_air = p_air_0 * (V_air_0 / V_air)^adiabatic_index_air * air_left_coeff;

    % Mass flow in different stages of flight
    is_empty = m_e_water >= m_fuel_0;
    if is_empty
        if p_air > p_atm
            delta_p = p_air - p_atm;
            if delta_p < 0 % Prevent backwards flow weirdness
                delta_p = 0;
                p_air = p_atm;
            end
            v_e_air_vec = -local_up * sqrt(2*n/(n+1) * delta_p/density_air);
            m_flow_air = density_air * norm(v_e_air_vec) * A_nozzle;
        end
        m_flow_water = 0;
    else
        % Water pressure
        air_height = V_air / A_cross_section;
        h = max(0, bottle_height - air_height);
        p_water_depth = density_w * g * h;
        p_w = p_air + p_water_depth;
        p_e = p_atm;
    
        density_e = density_w; % Assume laminar flow
        density_ratio = density_e / density_w;
        area_ratio = A_nozzle / A_cross_section;
        % Exhaust velocity relative to rocket
        v_e_water_vec = -local_up * C_discharge * sqrt( ...
            2 * adiabatic_index_w / (adiabatic_index_w - 1) * ...
            (p_e/density_e - p_w/density_w) / ...
            (density_ratio^2 * area_ratio^2 - 1) ...
        );

        m_flow_water = density_w * norm(v_e_water_vec) * A_nozzle;
        m_flow_air = 0;
    end
    delta_V_air = m_flow_water / density_w;

    % Aerodynamics
    v_rel_vec = v_vec - wind; % Velocity relative to wind
    p_dynamic = 0.5 * density_amb_air * norm(v_rel_vec)^2;

    % Forces
    F_prop_w = -m_flow_water * v_e_water_vec;
    F_prop_air = -m_flow_air * v_e_air_vec; % + delta_p*A_nozzle*local_up;
    F_g = [0; -g] * m_tot;
    F_d = -C_drag * density_amb_air * p_dynamic * A_cross_section * DirectionOf(v_rel_vec);

    % Total
    F_var_mass = F_prop_w + F_prop_air; % Force from mass change
    F_ext = F_g + F_d; % Sum of external forces
    F_tot = F_var_mass + F_ext; % Net force

    a_vec = F_tot ./ m_tot;

    dy(1) = m_flow_water;
    dy(2) = m_flow_air;
    dy(3) = delta_V_air;
    dy(4) = a_vec(1);
    dy(5) = a_vec(2);
    dy = dy';
    %disp(dy);
end

