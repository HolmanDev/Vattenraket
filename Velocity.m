function [v_vec] = Velocity(N, dt, m_rocket, m_fuel, g, p_0, p_air, density_water, A_nozzle, cd)
    v_vec(:, 1) = [0; 0];
    m_0 = m_rocket + m_fuel;
    a_g = [0; -g];
    m_e = 0;
    v_e = [0; 0];
    local_up = 1/sqrt(2)*[1; 1]; % Should pitch over time
    for i=1:N
        delta_p = p_air-p_0; % p_air should decrease

        v_e_last = v_e;
        v_e = cd * sqrt(2*delta_p/density_water) * local_up;
        a_e = (v_e - v_e_last) ./ dt;

        m_flow = norm(v_e) * density_water * A_nozzle;
        m_e = m_e + m_flow*dt;
        is_empty = m_e >= m_fuel;
        if is_empty
            m_e = m_fuel;
            m_flow = 0;
        end

        a_luft = [0; 0]; % Ändra sen

        % Euler framåt
        v_vec(:, i+1) = v_vec(:, i) + dt*(m_flow*(v_e+v_vec(:, i)) - m_e*a_e)...
            ./ (m_0 - m_e) + a_g*dt + a_luft*dt;
    end
end