function [a_vec, v_vec, s_vec, m_flow_vec] = FlightIntegral(N, dt, V_air_0, m_rocket, m_fuel, p_atm, p_air_0, density_water, A_nozzle, C_discharge, angle)
    g = 9.82; % [m/s^2]
    a_vec(:, 1) = [0; 0];
    v_vec(:, 1) = [0; 0];
    s_vec(:, 1) = [0; 0];
    m_0 = m_rocket + m_fuel;
    m_e = 0;
    local_up = [cosd(angle); sind(angle)];
    m_flow_vec = zeros(1, N);
    p_air = p_air_0;
    v_e_air = [0; 0];
    m_e_air = 0;
    V_air = V_air_0;
    T = -1 + 273.15;
    RH = 0.85;
    R_spec_water_vapor = 461.4;
    R_spec_dry_air = 286.9;
    R_spec = R_spec_dry_air + (R_spec_water_vapor-R_spec_dry_air) * RH;
    m_flow_air = 0;
    n = 1.135;
    tailwind = 0; % m/s
    water_vapour_gamma = 1.37; % Adiabatic index
    %r = (2 / (n + 1)) ^ (n/(n-1));
    v_e = 0;
    a_e = 0;
    volume_air_loss = 0;
    for i=1:N
        if m_e < m_fuel
            V_air = V_air_0 + m_e / density_water - volume_air_loss;
            p_air = p_air_0 * (V_air_0/V_air)^water_vapour_gamma;
            delta_p = p_air - p_atm; % p_air borde minska
        end
        density_air = p_air / (R_spec * T);
        
        % Exhaust-hastighet. OBS: Relativ till raketen!  
        beta = 20.5/88.1; % flaskhalsdiameter / flaskdiameter
        % Kom ihåg vattentrycket också!
        A_front = 0.04405^2 * pi;
        delta_z = max(0, 0.34 - V_air / A_front); % Ska vi ha med detta?
        extra_water_pressure = density_water * 9.82 * delta_z;
        nozzle_delta_p = delta_p + extra_water_pressure;
        turbulance_ratio = 0.8 + 0.2 * p_atm/p_air; % Hur turbulent är exhausten? 1 = vatten, 0.1 = luft
        density_exhaust = density_water * turbulance_ratio; % Är exhausten vatten eller vattenånga?        
        %v_e = C_discharge * sqrt(2*nozzle_delta_p / (density_exhaust * (1 - beta^4))) * -local_up;
        p_water = p_air + extra_water_pressure;
        p_exhaust = p_atm;
        k = water_vapour_gamma / (water_vapour_gamma - 1);
        v_e = -local_up * C_discharge * sqrt(2*k*(p_exhaust/density_exhaust - p_water/density_water) / ((density_exhaust/density_water)^2 * beta^4 - 1)); % I think this is the best equation
        %v_e = -local_up * C_discharge * sqrt(2 * k * (p_water/density_water - p_exhaust/density_exhaust));
        % Using (scuffed) thermodynamics:
        %old_v_e = v_e;
        %v_e = -local_up * sqrt(abs(2/(3*(1 - beta^4)) * (k * (p_water/density_water - p_exhaust/density_exhaust) - norm(a_e)*delta_z))); % cursed absolute value
        %a_e = v_e - old_v_e;
        %disp(v_e)
        %v_e = sqrt(2 * delta_p / density_water) * -local_up; <- Bernoulli
        %disp(v_e)
        % Borde inte flaskans acceleration driva ut vätskan snabbare?

        exhaust_volume = norm(v_e) * A_nozzle * dt;
        volume_air_loss = (1 - turbulance_ratio) * exhaust_volume; % Volume of ejected air 

        % Massflöde       
        m_flow = norm(v_e) * A_nozzle * density_exhaust;
        m_e = m_e + m_flow*dt;
        is_empty = m_e >= m_fuel;
        %disp(p_air)
        if is_empty
            m_e = m_fuel;
            m_flow = 0;
            
            % https://en.wikipedia.org/wiki/Density_of_air
            % https://www.engineersedge.com/fluid_flow/convergent_nozzle_flow_velocity_14032.htm
            if p_air > p_atm
                %r = p_air/p_atm;
                v_e_air = sqrt(2*n/(n+1) * delta_p/density_air) * -local_up; %sqrt(2*(p_air-p_atm)/density_air) * -local_up; % Lite osäker på v_e_air...
                m_flow_air = norm(v_e_air) * density_air * A_nozzle;
                m_e_air = m_e_air + m_flow_air * dt;
                p_air = p_air_0 * V_air_0/V_air * (1 - m_e_air/air_mass_at_cutoff);
                delta_p = p_air - p_atm;
                if delta_p < 0
                    p_air = p_atm;
                    v_e_air = [0; 0];
                    m_flow_air = 0;
                    delta_p = 0;
                end
            else
                v_e_air = [0; 0];
                m_flow_air = 0;
            end
        else
            air_mass_at_cutoff = V_air * density_air;
        end
        m_flow_vec(i) = m_flow + m_flow_air;

        % Krafter
        % RÄKNA MED KRAFTEN FRÅN HÖGTRYCKSOMRÅDET UNDER RAKETEN
        %disp(delta_p)
        F_prop_water = -m_flow*v_e;
        F_prop_air = -m_flow_air*v_e_air; % + delta_p*A_nozzle*local_up;
        F_g = [0; -g] * (m_0 - m_e);
        drag_coeff = 0.4;
        v_tailwind = v_vec(:, i) - [sign(v_vec(1, i)) * tailwind; 0];

        F_air = -drag_coeff * 1.3 * 0.5 * A_front * norm(v_tailwind) * v_tailwind;

        % Euler framåt
        m_fuel_left = m_fuel - m_e;
        % (F - v*dm/dt) / m = a
        F_var_mass = F_prop_water + F_prop_air; % v * dm/dt
        F_ext = F_g + F_air;
        a = (F_ext + F_var_mass) ./ (m_rocket + m_fuel_left);
        a_vec(:, i+1) = a;
        v_vec(:, i+1) = v_vec(:, i) + dt*a;
        s_vec(:, i+1) = s_vec(:, i) + v_vec(:, i)*dt + 0.5*a_vec(:, i)*dt*dt;

        %disp(s_vec(2, i+1))
        if(s_vec(2, i+1) < 0)
            break
        end
        
        local_up = v_vec(:, i+1)/norm(v_vec(:, i+1));
    end
end