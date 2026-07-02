function main_hydrotreating_simulation_CS()
    clear;
    clc;
    close all;
    fprintf('Starting CHE251 (Group 1) Process Simulation...\n');
    fprintf('--- Using Chao-Seader (VLE) + Cp(T) Enthalpy Model ---\n');

    [constants, specs] = populate_constants();
    
    diesel_feed.F_comp = zeros(1, constants.num_comps); 
    diesel_feed.T = specs.T_inlet;
    diesel_feed.P = specs.P_inlet;
    
    F_total_diesel = 2.226; 
    z_CH3SH = 0.058;        
    
    diesel_feed.F_comp(1) = F_total_diesel * z_CH3SH;
    
    num_inerts = 8;
    z_inert_each = (1 - z_CH3SH) / num_inerts;
    diesel_feed.F_comp(5:11) = F_total_diesel * z_inert_each; 
    diesel_feed.F_comp(5) = diesel_feed.F_comp(5) + (F_total_diesel * z_inert_each); 
    
    diesel_feed.H = get_stream_enthalpy(diesel_feed, 'L', constants, specs);
    
    recycle_guess.F_comp = zeros(1, constants.num_comps);
    recycle_guess.F_comp(2) = 150; 
    recycle_guess.T = specs.T_cooler_out; 
    recycle_guess.P = specs.P_flash; 
    recycle_guess.H = get_stream_enthalpy(recycle_guess, 'V', constants, specs);
    
    makeup_H2_flow_guess = 120; 
    max_iter = 100;
    tol = 1e-4;
    damping_factor = 0.5; 
    for iter = 1:max_iter
        fprintf('--- Iteration %d ---\n', iter);
        
        recycle_old = recycle_guess;
        makeup_H2_old = makeup_H2_flow_guess;
        [pump_out, W_pump] = unit_pump(diesel_feed, specs.P_system, ...
                                        specs.eta_pump, constants, specs);
        makeup_H2.F_comp = zeros(1, constants.num_comps);
        makeup_H2.F_comp(2) = makeup_H2_flow_guess;
        makeup_H2.T = specs.T_H2_makeup;
        makeup_H2.P = specs.P_system; 
        makeup_H2.H = get_stream_enthalpy(makeup_H2, 'V', constants, specs);
        mixer_out = unit_mixer(pump_out, makeup_H2, recycle_guess, constants, specs);
        F_S_feed = pump_out.F_comp(1); 
        F_H2_target = F_S_feed * specs.H2_Sulfur_Ratio;
        F_H2_recycle = recycle_guess.F_comp(2);
        F_H2_makeup_new = F_H2_target - F_H2_recycle;
        if F_H2_makeup_new < 0
            F_H2_makeup_new = 0; 
        end
        makeup_H2_flow_guess = (1 - damping_factor) * makeup_H2_old + ...
                               damping_factor * F_H2_makeup_new;
        [furnace_out, Q_furnace] = unit_furnace(mixer_out, specs.T_reactor, ...
                                               'V', constants, specs);
        [reactor_out, W_cat, Q_reactor, chi_S] = ...
            unit_pfr(furnace_out, specs, constants);
        cooler_out = reactor_out;
        cooler_out.T = specs.T_cooler_out;
        
        [vapor_out, liquid_out] = unit_flash_sep(cooler_out, constants, specs);
        
        Q_cooler = unit_cooler(reactor_out, vapor_out, liquid_out, constants, specs);
        [purge_stream, recycle_pre_comp] = ...
            unit_splitter(vapor_out, specs.split_ratio, constants);
        [recycle_calculated, W_comp] = ...
            unit_compressor(recycle_pre_comp, specs.P_system, ...
                            specs.eta_comp, constants, specs);
        err_H2 = abs(makeup_H2_flow_guess - makeup_H2_old) / ...
                 (makeup_H2_old + tol);
        err_recycle = sum(abs(recycle_calculated.F_comp - recycle_old.F_comp)) / ...
                      (sum(recycle_old.F_comp) + tol);
        
        fprintf('  Makeup H2 Flow: %.2f kmol/h (Err: %.2e)\n', ...
                makeup_H2_flow_guess, err_H2);
        fprintf('  Recycle H2 Flow: %.2f kmol/h (Err: %.2e)\n', ...
                recycle_calculated.F_comp(2), err_recycle);
        
        if err_H2 < tol && err_recycle < tol
            fprintf('\n*** CONVERGENCE ACHIEVED IN %d ITERATIONS ***\n', iter);
            break;
        end
        
        recycle_guess.F_comp = (1 - damping_factor) * recycle_old.F_comp + ...
                               damping_factor * recycle_calculated.F_comp;
        recycle_guess.T = recycle_calculated.T;
        recycle_guess.P = recycle_calculated.P;
        recycle_guess.H = get_stream_enthalpy(recycle_guess, 'V', constants, specs);
        
        if iter == max_iter
            fprintf('\n*** FAILED TO CONVERGE AFTER %d ITERATIONS ***\n', iter);
        end
    end
    
    fprintf('\n--- Base Case Results ---\n');
    fprintf('Core Performance:\n');
    fprintf('  Single-Pass Conversion (chi_S): %.2f %%\n', chi_S * 100);
    
    % --- ADDED OVERALL CONVERSION CALCULATION ---
    sulfur_in = diesel_feed.F_comp(1);
    sulfur_out_liq = liquid_out.F_comp(1);
    sulfur_out_vap = purge_stream.F_comp(1);
    overall_conversion = (sulfur_in - (sulfur_out_liq + sulfur_out_vap)) / sulfur_in;
    fprintf('  Overall Conversion:           %.2f %%\n', overall_conversion * 100);
    % --- END ADDITION ---
    
    fprintf('\nDesign Parameters:\n');
    fprintf('  Catalyst Weight (W_cat): %.2f kg\n', W_cat / 1000); 
    
    fprintf('\nEnergy Duties (kW):\n');
    fprintf('  Pump Work (W_pump):       %.2f\n', W_pump);
    fprintf('  Furnace Duty (Q_furnace):   %.2f\n', Q_furnace);
    fprintf('  Reactor Duty (Q_reactor):   %.2f\n', Q_reactor);
    fprintf('  Cooler Duty (Q_cooler):     %.2f\n', Q_cooler);
    fprintf('  Compressor Work (W_comp):   %.2f\n', W_comp);
    fprintf('\nStream Properties (T, H, F_total):\n');
    fprintf('  Stream              T (C)   H (kJ/kmol)  Flow (kmol/h)\n');
    fprintf('  --------------------------------------------------------\n');
    fprintf('  1. Diesel Feed      %.2f    %.2f         %.2f\n', diesel_feed.T, diesel_feed.H, sum(diesel_feed.F_comp));
    fprintf('  2. Pump Out         %.2f    %.2f         %.2f\n', pump_out.T, pump_out.H, sum(pump_out.F_comp));
    fprintf('  3. H2 Makeup        %.2f    %.2f         %.2f\n', makeup_H2.T, makeup_H2.H, sum(makeup_H2.F_comp));
    fprintf('  4. Recycle In       %.2f    %.2f         %.2f\n', recycle_calculated.T, recycle_calculated.H, sum(recycle_calculated.F_comp));
    fprintf('  5. Furnace Out      %.2f    %.2f         %.2f\n', furnace_out.T, furnace_out.H, sum(furnace_out.F_comp));
    fprintf('  6. Reactor Out      %.2f    %.2f         %.2f\n', reactor_out.T, reactor_out.H, sum(reactor_out.F_comp));
    fprintf('  7. Product (Liquid) %.2f    %.2f         %.2f\n', liquid_out.T, liquid_out.H, sum(liquid_out.F_comp));
    fprintf('  8. Purge (Vapor)    %.2f    %.2f         %.2f\n', purge_stream.T, purge_stream.H, sum(purge_stream.F_comp));
    
    fprintf('\nStream Flows & Compositions (kmol/h):\n');
    fprintf('  Component       Feed     Product (Liquid)   Purge (Vapor)\n');
    fprintf('  ---------------------------------------------------------\n');
    for i = 1:constants.num_comps
        fprintf('  %-15s %-8.3f %-18.3f %-12.3f\n', ...
                constants.comp_names{i}, ...
                diesel_feed.F_comp(i), ...
                liquid_out.F_comp(i), ...
                purge_stream.F_comp(i));
    end
    
end
% =========================================================================
% == DATA AND PROPERTY FUNCTIONS (Chao-Seader + Cp(T))
% =========================================================================
function [constants, specs] = populate_constants()
    
    constants.comp_names = {
        'CH3SH', 'H2', 'CH4', 'H2S', 'n-Hexane', 'n-Decane', ...
        'n-Hexadecane', 'Benzene', 'Toluene', 'P-xylene', 'Cyclohexane'
    };
    constants.num_comps = length(constants.comp_names);
    
    % --- Cp(T) Data ---
    % Cp_vap = A + B*T + C*T^2 + D*T^3 
    % T is in Kelvin, Cp will be in J/(mol*K)
    % Cp_liq is constant kJ/(kmol*K)
    i = 1; % CH3SH
    comps(i).Name = 'CH3SH';
    comps(i).MW = 48.11;
    comps(i).rho_liq = 867;
    comps(i).Cp_liq = 100.0;
    comps(i).Hvap = 24.0e3;
    comps(i).T_boil = 6.2; % degC
    comps(i).Antoine = [4.195, 1075.1, 231.9]; 
    comps(i).Cp_vap_coeffs = [30.40, 9.615e-2, -1.91e-5, -1.79e-9]; % A, B, C, D
    
    i = 2; % H2
    comps(i).Name = 'H2';
    comps(i).MW = 2.016;
    comps(i).rho_liq = 70.8; 
    comps(i).Cp_liq = 28.0;
    comps(i).Hvap = 0.9e3;
    comps(i).T_boil = -252.9;
    comps(i).Antoine = []; 
    comps(i).Cp_vap_coeffs = [27.14, 9.274e-3, -1.381e-5, 7.645e-9];
    
    i = 3; % CH4
    comps(i).Name = 'CH4';
    comps(i).MW = 16.04;
    comps(i).rho_liq = 422;
    comps(i).Cp_liq = 52.7;
    comps(i).Hvap = 8.2e3;
    comps(i).T_boil = -161.5;
    comps(i).Antoine = []; 
    comps(i).Cp_vap_coeffs = [34.31, 5.469e-2, 0.366e-5, -1.10e-8];
    
    i = 4; % H2S
    comps(i).Name = 'H2S';
    comps(i).MW = 34.08;
    comps(i).rho_liq = 993;
    comps(i).Cp_liq = 65.0;
    comps(i).Hvap = 18.7e3;
    comps(i).T_boil = -60.0;
    comps(i).Antoine = [4.162, 946.8, 246.5];
    comps(i).Cp_vap_coeffs = [32.69, 1.246e-2, 1.913e-5, -1.41e-8];
    
    i = 5; % n-Hexane
    comps(i).Name = 'n-Hexane';
    comps(i).MW = 86.18;
    comps(i).rho_liq = 655;
    comps(i).Cp_liq = 195.0;
    comps(i).Hvap = 28.9e3;
    comps(i).T_boil = 69.0;
    comps(i).Antoine = [4.002, 1171.5, 224.4];
    comps(i).Cp_vap_coeffs = [-4.43, 3.95e-1, -2.23e-4, 4.85e-8];
    i = 6; % n-Decane
    comps(i).Name = 'n-Decane';
    comps(i).MW = 142.28;
    comps(i).rho_liq = 730;
    comps(i).Cp_liq = 313.0;
    comps(i).Hvap = 51.4e3;
    comps(i).T_boil = 174.0;
    comps(i).Antoine = [4.058, 1504.1, 209.7];
    comps(i).Cp_vap_coeffs = [-12.33, 6.79e-1, -4.20e-4, 1.05e-7];
    
    i = 7; % n-Hexadecane
    comps(i).Name = 'n-Hexadecane';
    comps(i).MW = 226.44;
    comps(i).rho_liq = 773;
    comps(i).Cp_liq = 470.0;
    comps(i).Hvap = 82.0e3;
    comps(i).T_boil = 287.0;
    comps(i).Antoine = [4.135, 1858.1, 175.7];
    comps(i).Cp_vap_coeffs = [-24.10, 1.06, -6.64e-4, 1.70e-7];
    
    i = 8; % Benzene
    comps(i).Name = 'Benzene';
    comps(i).MW = 78.11;
    comps(i).rho_liq = 874;
    comps(i).Cp_liq = 136.0;
    comps(i).Hvap = 30.7e3;
    comps(i).T_boil = 80.1;
    comps(i).Antoine = [4.018, 1211.0, 220.8];
    comps(i).Cp_vap_coeffs = [-3.39, 3.23e-1, -2.48e-4, 7.63e-8];
    
    i = 9; % Toluene
    comps(i).Name = 'Toluene';
    comps(i).MW = 92.14;
    comps(i).rho_liq = 867;
    comps(i).Cp_liq = 156.0;
    comps(i).Hvap = 33.5e3;
    comps(i).T_boil = 111.0;
    comps(i).Antoine = [4.080, 1344.8, 219.5];
    comps(i).Cp_vap_coeffs = [-1.75, 4.04e-1, -3.00e-4, 9.38e-8];
    
    i = 10; % P-xylene
    comps(i).Name = 'P-xylene';
    comps(i).MW = 106.16;
    comps(i).rho_liq = 861;
    comps(i).Cp_liq = 175.0;
    comps(i).Hvap = 35.8e3;
    comps(i).T_boil = 138.0;
    comps(i).Antoine = [4.106, 1453.4, 215.3];
    comps(i).Cp_vap_coeffs = [-3.31, 4.75e-1, -3.51e-4, 1.09e-7];
    
    i = 11; % Cyclohexane
    comps(i).Name = 'Cyclohexane';
    comps(i).MW = 84.16;
    comps(i).rho_liq = 779;
    comps(i).Cp_liq = 156.0;
    comps(i).Hvap = 29.9e3;
    comps(i).T_boil = 80.7;
    comps(i).Antoine = [3.966, 1182.2, 222.9];
    comps(i).Cp_vap_coeffs = [-43.47, 5.04e-1, -3.33e-4, 8.87e-8];
    
    constants.comps = comps;
    constants = populate_thermo_data(constants);
    
    constants.R_gas_kcal = 1.987e-3; 
    constants.R_gas_L_MPa = 8.314e-3; 
    constants.R_gas_J_mol_K = 8.314;  
    constants.R_gas_L_bar = 0.08314;  
    
    constants.Ea = 30.0;             
    constants.ln_k0 = 25.0;          
    constants.k0 = exp(constants.ln_k0);
    
    constants.H_rxn = -72.5e3; 
    constants.stoich_vec = [-1, -1, 1, 1, 0, 0, 0, 0, 0, 0, 0];
    
    specs.P_system = 12.5;         
    
    % --- MODIFICATIONS FOR 85%+ CONVERSION ---
    specs.T_reactor = 400.0;       % Changed from 375.0
    specs.LHSV = 2.0;              % Changed from 2.9
    specs.split_ratio = 0.85;      % Changed from 0.75
    % --- END MODIFICATIONS ---
    
    specs.H2_Sulfur_Ratio = 1000;  % THIS IS THE LINE THAT WAS MISSING
    specs.eta_pump = 0.75;         
    specs.eta_comp = 0.75;         
    
    specs.T_ref = 25.0;            
    specs.T_inlet = 25.0;          
    specs.P_inlet = 0.5;           
    specs.T_H2_makeup = 25.0;      
    specs.T_cooler_out = 40.0;     
    
    specs.P_flash = 12.0;          
end
function constants = populate_thermo_data(constants)
    
    i = 1; 
    constants.comps(i).Tc = 469.0;
    constants.comps(i).Pc = 72.4;
    constants.comps(i).omega = 0.198;
    constants.comps(i).V_L = 0.061;
    constants.comps(i).delta = 18.2;
    
    i = 2; 
    constants.comps(i).Tc = 33.2;
    constants.comps(i).Pc = 13.1;
    constants.comps(i).omega = -0.21;
    constants.comps(i).V_L = 0.028;
    constants.comps(i).delta = 3.7;
    
    i = 3; 
    constants.comps(i).Tc = 190.6;
    constants.comps(i).Pc = 46.0;
    constants.comps(i).omega = 0.011;
    constants.comps(i).V_L = 0.038;
    constants.comps(i).delta = 10.6;
    
    i = 4; 
    constants.comps(i).Tc = 373.2;
    constants.comps(i).Pc = 89.4;
    constants.comps(i).omega = 0.100;
    constants.comps(i).V_L = 0.046;
    constants.comps(i).delta = 19.0;
    
    i = 5; 
    constants.comps(i).Tc = 507.6;
    constants.comps(i).Pc = 30.2;
    constants.comps(i).omega = 0.301;
    constants.comps(i).V_L = 0.132;
    constants.comps(i).delta = 14.9;
    i = 6; 
    constants.comps(i).Tc = 617.7;
    constants.comps(i).Pc = 21.0;
    constants.comps(i).omega = 0.489;
    constants.comps(i).V_L = 0.230;
    constants.comps(i).delta = 15.8;
    
    i = 7; 
    constants.comps(i).Tc = 722.0;
    constants.comps(i).Pc = 15.3;
    constants.comps(i).omega = 0.741;
    constants.comps(i).V_L = 0.347;
    constants.comps(i).delta = 16.3;
    
    i = 8; 
    constants.comps(i).Tc = 562.2;
    constants.comps(i).Pc = 48.9;
    constants.comps(i).omega = 0.210;
    constants.comps(i).V_L = 0.089;
    constants.comps(i).delta = 18.7;
    
    i = 9; 
    constants.comps(i).Tc = 591.8;
    constants.comps(i).Pc = 41.1;
    constants.comps(i).omega = 0.263;
    constants.comps(i).V_L = 0.107;
    constants.comps(i).delta = 18.2;
    
    i = 10; 
    constants.comps(i).Tc = 616.2;
    constants.comps(i).Pc = 35.1;
    constants.comps(i).omega = 0.321;
    constants.comps(i).V_L = 0.124;
    constants.comps(i).delta = 18.0;
    
    i = 11; 
    constants.comps(i).Tc = 553.6;
    constants.comps(i).Pc = 40.7;
    constants.comps(i).omega = 0.210;
    constants.comps(i).V_L = 0.109;
    constants.comps(i).delta = 16.8;
end
function H_integral = integral_Cp_vap(T_K, coeffs)
    % Calculates the integral of Cp_vap(T) dT from 0 to T_K
    % Cp_vap = A + B*T + C*T^2 + D*T^3
    % Integral = A*T + (B/2)*T^2 + (C/3)*T^3 + (D/4)*T^4
    A = coeffs(1);
    B = coeffs(2);
    C = coeffs(3);
    D = coeffs(4);
    
    H_integral = A*T_K + (B/2)*(T_K^2) + (C/3)*(T_K^3) + (D/4)*(T_K^4);
    % Result is in J/mol. Convert to kJ/kmol (J/mol * 1000 mol/kmol / 1000 J/kJ)
    % This is a 1:1 conversion, so J/mol = kJ/kmol
end
function H_comp = get_component_enthalpy(T_C, phase, comp_data, T_ref)
    
    T_K = T_C + 273.15;
    T_ref_K = T_ref + 273.15;
    
    % Reference State: Liquid at T_ref
    % H_ref = 0
    
    if strcmpi(phase, 'L')
        % Enthalpy of Liquid at T_C
        % H(T,L) = integral(Cp_liq, T_ref, T)
        % We use constant Cp_liq (in kJ/kmol*K)
        H_comp = comp_data.Cp_liq * (T_K - T_ref_K);
    else
        % Enthalpy of Vapor at T_C
        % H(T,V) = H(T_boil,L) + H_vap + H(T,V)
        % H(T,V) = [integral(Cp_liq, T_ref, T_boil)] + [H_vap] + [integral(Cp_vap, T_boil, T)]
        
        T_boil_K = comp_data.T_boil + 273.15;
        
        % 1. Enthalpy of liquid from T_ref to T_boil
        H_liq_sensible = comp_data.Cp_liq * (T_boil_K - T_ref_K);
        
        % 2. Heat of Vaporization (kJ/kmol)
        H_vap = comp_data.Hvap;
        
        % 3. Enthalpy of vapor from T_boil to T
        % Get integral(Cp_vap) from T_boil to T_K
        coeffs = comp_data.Cp_vap_coeffs;
        H_vap_sensible = integral_Cp_vap(T_K, coeffs) - integral_Cp_vap(T_boil_K, coeffs);
        
        % Total enthalpy (kJ/kmol)
        H_comp = H_liq_sensible + H_vap + H_vap_sensible;
        
        % Handle light gases (H2, CH4) where T_boil < T_ref
        if T_boil_K < T_ref_K
            % H(T,V) = H(T_ref,V) + integral(Cp_vap, T_ref, T)
            % H(T_ref,V) = H(T_boil,L) + H_vap + integral(Cp_vap, T_boil, T_ref)
            H_vap_at_ref = H_liq_sensible + H_vap + ...
                           (integral_Cp_vap(T_ref_K, coeffs) - integral_Cp_vap(T_boil_K, coeffs));
            
            H_vap_sensible_from_ref = integral_Cp_vap(T_K, coeffs) - integral_Cp_vap(T_ref_K, coeffs);
            
            H_comp = H_vap_at_ref + H_vap_sensible_from_ref;
        end
    end
end
function H_stream = get_stream_enthalpy(stream, phase, constants, specs)
    F_total = sum(stream.F_comp);
    if F_total == 0
        H_stream = 0;
        return;
    end
    
    H_total_kJ_h = 0;
    for i = 1:constants.num_comps
        if stream.F_comp(i) > 0
            H_i = get_component_enthalpy(stream.T, phase, constants.comps(i), specs.T_ref);
            H_total_kJ_h = H_total_kJ_h + stream.F_comp(i) * H_i;
        end
    end
    H_stream = H_total_kJ_h / F_total; 
end
function T_out = get_temp_from_enthalpy(stream, H_in, phase, constants, specs)
    
    F_total = sum(stream.F_comp);
    if F_total == 0
        T_out = specs.T_ref;
        return;
    end
    
    % This is now the inverse of a complex integral function
    err_func = @(T_C) (get_stream_enthalpy(struct('T', T_C, 'F_comp', stream.F_comp), phase, constants, specs) - H_in);
    
    % Use fzero to find the temperature T_C that matches the enthalpy H_in
    options = optimset('Display','off');
    try
        T_out = fzero(err_func, stream.T, options); % Use last T as guess
    catch
        T_out = fzero(err_func, specs.T_ref, options); % Fallback guess
    end
end
function V_liq_L_h = get_liquid_vol_flow(stream, constants)
    V_liq_m3_h = 0;
    for i = 1:constants.num_comps
        F_kmol_h = stream.F_comp(i);
        MW_kg_kmol = constants.comps(i).MW;
        rho_liq_kg_m3 = constants.comps(i).rho_liq;
        
        Mass_flow_kg_h = F_kmol_h * MW_kg_kmol;
        V_liq_m3_h = V_liq_m3_h + (Mass_flow_kg_h / rho_liq_kg_m3);
    end
    V_liq_L_h = V_liq_m3_h * 1000;
end
function [V_F, x, y] = solve_rachford_rice(z, K)
    g_func = @(psi) sum(z .* (K - 1) ./ (1 + psi .* (K - 1)));
    
    psi_min = 1 / (1 - max(K));
    psi_max = 1 / (1 - min(K));
    
    psi_low = max(0, psi_min + 1e-5);
    psi_high = min(1, psi_max - 1e-5);
    
    if psi_low >= psi_high
       if g_func(0.01) > 0 
           V_F = 0;
       else 
           V_F = 1;
       end
    else
        options = optimset('Display','off');
        V_F = fzero(g_func, [psi_low, psi_high], options);
    end
    x = z ./ (1 + V_F .* (K - 1));
    y = K .* x;
    
    x = x / sum(x);
    y = y / sum(y);
end
function [stream_out, W_pump] = unit_pump(stream_in, P_out, eta, constants, specs)
    stream_out = stream_in;
    stream_out.P = P_out;
    
    V_liq_L_h = get_liquid_vol_flow(stream_in, constants);
    V_liq_m3_s = V_liq_L_h / 1000 / 3600;
    
    P_in_Pa = stream_in.P * 1e6;
    P_out_Pa = P_out * 1e6;
    
    W_ideal_W = V_liq_m3_s * (P_out_Pa - P_in_Pa);
    W_pump_W = W_ideal_W / eta;
    W_pump = W_pump_W / 1000; 
    
    W_pump_kJ_h = W_pump * 3600;
    H_in_kJ_kmol = get_stream_enthalpy(stream_in, 'L', constants, specs);
    F_total = sum(stream_in.F_comp);
    
    H_out_kJ_kmol = H_in_kJ_kmol + (W_pump_kJ_h / F_total);
    stream_out.H = H_out_kJ_kmol;
    
    stream_out.T = get_temp_from_enthalpy(stream_out, H_out_kJ_kmol, 'L', constants, specs);
end
function stream_out = unit_mixer(stream_1, stream_2, stream_3, constants, specs)
    
    stream_out.F_comp = stream_1.F_comp + stream_2.F_comp + stream_3.F_comp;
    F_total_out = sum(stream_out.F_comp);
    stream_out.P = stream_1.P; 
    
    H1_kJ_h = get_stream_enthalpy(stream_1, 'L', constants, specs) * sum(stream_1.F_comp);
    H2_kJ_h = get_stream_enthalpy(stream_2, 'V', constants, specs) * sum(stream_2.F_comp);
    H3_kJ_h = get_stream_enthalpy(stream_3, 'V', constants, specs) * sum(stream_3.F_comp);
    
    H_total_in = H1_kJ_h + H2_kJ_h + H3_kJ_h;
    stream_out.H = H_total_in / F_total_out;
    
    stream_out.T = (stream_1.T + stream_2.T + stream_3.T) / 3; 
end
function [stream_out, Q_furnace] = unit_furnace(stream_in, T_out, phase_out, constants, specs)
    stream_out = stream_in;
    stream_out.T = T_out;
    
    H_in_kJ_kmol = stream_in.H;
    H_in_kJ_h = H_in_kJ_kmol * sum(stream_in.F_comp);
    
    H_out_kJ_kmol = get_stream_enthalpy(stream_out, phase_out, constants, specs);
    H_out_kJ_h = H_out_kJ_kmol * sum(stream_out.F_comp);
    
    Q_kJ_h = H_out_kJ_h - H_in_kJ_h;
    Q_furnace = Q_kJ_h / 3600; 
end
function dF_dW = pfr_odes(W, F_kmol_h, T_K, P_MPa, constants)
    
    F_kmol_h(F_kmol_h < 0) = 0; 
    F_mol_h = F_kmol_h * 1000;
    F_total_mol_h = sum(F_mol_h);
    
    if F_total_mol_h == 0
        dF_dW = zeros(constants.num_comps, 1);
        return;
    end
    
    y = F_mol_h / F_total_mol_h; 
    
    K_rate_const = constants.k0 * exp(-constants.Ea / (constants.R_gas_kcal * T_K));
    
    C_S = (y(1) * P_MPa) / (constants.R_gas_L_MPa * T_K); 
    
    P_H2 = y(2) * P_MPa; 
    
    if C_S > 0 && P_H2 > 0
        rate_mol_gcat_h = K_rate_const * (C_S^1.0) * (P_H2^0.4);
    else
        rate_mol_gcat_h = 0;
    end
    
    rate_kmol_gcat_h = rate_mol_gcat_h / 1000;
    
    r_kmol_gcat_h = constants.stoich_vec * rate_kmol_gcat_h;
    
    dF_dW = r_kmol_gcat_h(:); 
end
function [stream_out, W_cat, Q_reactor, chi_S] = ...
         unit_pfr(stream_in, specs, constants)
    
    V_liq_feed_L_h = get_liquid_vol_flow(stream_in, constants);
    W_cat = V_liq_feed_L_h / specs.LHSV; 
    
    F_in_kmol_h = stream_in.F_comp;
    T_K = stream_in.T + 273.15;
    P_MPa = stream_in.P;
    
    ode_options = odeset('RelTol', 1e-6, 'NonNegative', 1:constants.num_comps);
    [~, F_out_mat] = ode45(@pfr_odes, [0 W_cat], F_in_kmol_h, ode_options, ...
                           T_K, P_MPa, constants);
    
    F_out_kmol_h = F_out_mat(end, :);
    
    stream_out = stream_in;
    stream_out.F_comp = F_out_kmol_h;
    
    stream_out.H = get_stream_enthalpy(stream_out, 'V', constants, specs);
    
    F_S_in = F_in_kmol_h(1);
    F_S_out = F_out_kmol_h(1);
    chi_S = (F_S_in - F_S_out) / F_S_in;
    
    moles_reacted_kmol_h = F_S_in - F_S_out;
    Q_kJ_h = moles_reacted_kmol_h * constants.H_rxn;
    Q_reactor = Q_kJ_h / 3600; 
end
function Q_cooler = unit_cooler(stream_in, stream_vap_out, stream_liq_out, constants, specs)
    
    H_in_kJ_h = get_stream_enthalpy(stream_in, 'V', constants, specs) * sum(stream_in.F_comp);
    
    H_vap_out_kJ_h = get_stream_enthalpy(stream_vap_out, 'V', constants, specs) * sum(stream_vap_out.F_comp);
    H_liq_out_kJ_h = get_stream_enthalpy(stream_liq_out, 'L', constants, specs) * sum(stream_liq_out.F_comp);
    
    H_out_total_kJ_h = H_vap_out_kJ_h + H_liq_out_kJ_h;
    
    Q_kJ_h = H_out_total_kJ_h - H_in_kJ_h;
    Q_cooler = Q_kJ_h / 3600; 
end
function [vapor_out, liquid_out] = unit_flash_sep(stream_in, constants, specs)
    
    T_C = stream_in.T;
    P_MPa = specs.P_flash;
    z_feed = stream_in.F_comp / sum(stream_in.F_comp);
    
    max_iter = 50;
    tol = 1e-5;
    
    K = get_K_values_Raoult(T_C, P_MPa, constants);
    
    for i = 1:max_iter
        [V_F, x, y] = solve_rachford_rice(z_feed, K);
        
        K_new = get_K_values_ChaoSeader(T_C, P_MPa, x, y, constants);
        
        err = sum(abs(K_new - K) ./ K);
        if err < tol
            break;
        end
        
        K = 0.5 * K + 0.5 * K_new; 
        
        if i == max_iter
           warning('Flash separator did not converge on K-values');
        end
    end
    
    F_total_in = sum(stream_in.F_comp);
    V_total = V_F * F_total_in;
    L_total = (1 - V_F) * F_total_in;
    
    vapor_out.F_comp = y * V_total;
    vapor_out.T = T_C;
    vapor_out.P = P_MPa;
    vapor_out.H = get_stream_enthalpy(vapor_out, 'V', constants, specs);
    
    liquid_out.F_comp = x * L_total;
    liquid_out.T = T_C;
    liquid_out.P = P_MPa;
    liquid_out.H = get_stream_enthalpy(liquid_out, 'L', constants, specs);
end
function [purge_out, recycle_out] = unit_splitter(stream_in, S, constants)
    recycle_out.F_comp = stream_in.F_comp * S;
    purge_out.F_comp = stream_in.F_comp * (1 - S);
    
    recycle_out.T = stream_in.T;
    recycle_out.P = stream_in.P;
    recycle_out.H = stream_in.H;
    
    purge_out.T = stream_in.T;
    purge_out.P = stream_in.P;
    purge_out.H = stream_in.H;
end
function [stream_out, W_comp] = unit_compressor(stream_in, P_out, eta, constants, specs)
    stream_out = stream_in;
    stream_out.P = P_out;
    
    F_in_kmol_h = sum(stream_in.F_comp);
    if F_in_kmol_h == 0
        W_comp = 0;
        return;
    end
    
    Cp_avg_vap = 29.0; 
    Cv_avg_vap = Cp_avg_vap - 8.314;
    k = Cp_avg_vap / Cv_avg_vap;
    
    T_in_K = stream_in.T + 273.15;
    P_in_MPa = stream_in.P;
    P_out_MPa = P_out;
    if P_out_MPa <= P_in_MPa
        W_comp = 0;
        stream_out.T = stream_in.T;
        stream_out.H = stream_in.H;
        return;
    end
    
    T_out_ideal_K = T_in_K * (P_out_MPa / P_in_MPa)^((k - 1) / k);
    
    H_in_kJ_kmol = get_stream_enthalpy(stream_in, 'V', constants, specs);
    
    stream_ideal_out = stream_in;
    stream_ideal_out.T = T_out_ideal_K - 273.15;
    H_out_ideal_kJ_kmol = get_stream_enthalpy(stream_ideal_out, 'V', constants, specs);
    
    W_ideal_kJ_h = (H_out_ideal_kJ_kmol - H_in_kJ_kmol) * F_in_kmol_h;
    W_real_kJ_h = W_ideal_kJ_h / eta;
    W_comp = W_real_kJ_h / 3600; 
    
    H_out_real_kJ_kmol = H_in_kJ_kmol + (W_real_kJ_h / F_in_kmol_h);
    stream_out.H = H_out_real_kJ_kmol;
    stream_out.T = get_temp_from_enthalpy(stream_out, H_out_real_kJ_kmol, 'V', constants, specs);
end
function K = get_K_values_ChaoSeader(T_C, P_MPa, x, y, constants)
    
    T_K = T_C + 273.15;
    P_bar = P_MPa * 10;
    
    gamma = calculate_gamma_SH(T_K, x, constants);
    
    phi_hat = calculate_phi_hat_RK(T_K, P_bar, y, constants);
    
    K = zeros(1, constants.num_comps);
    for i = 1:constants.num_comps
        Psat_bar = calculate_Psat_Antoine(T_C, i, constants);
        
        phi_sat = calculate_phi_sat_RK(T_K, Psat_bar, i, constants);
        
        if phi_hat(i) > 1e-9 
            K(i) = (gamma(i) * phi_sat * Psat_bar) / (phi_hat(i) * P_bar);
        else
            K(i) = 1e-9;
        end
        
        if i == 2 || i == 3
            K(i) = 1000 / P_MPa; 
        end
    end
end
function gamma = calculate_gamma_SH(T_K, x, constants)
    
    n_comps = constants.num_comps;
    gamma = ones(1, n_comps);
    R_J_mol_K = constants.R_gas_J_mol_K;
    
    V_L = [constants.comps.V_L]; 
    delta = [constants.comps.delta]; 
    
    R_L_MPa = R_J_mol_K / 1e6 * 1e3; 
    
    V_m = sum(x .* V_L);
    if V_m == 0
        return; 
    end
    
    vol_frac = (x .* V_L) / V_m;
    
    delta_avg = sum(vol_frac .* delta);
    
    for i = 1:n_comps
        ln_gamma_i = (V_L(i) / (R_L_MPa * T_K)) * (delta_avg - delta(i))^2;
        gamma(i) = exp(ln_gamma_i);
    end
end
function phi_hat = calculate_phi_hat_RK(T_K, P_bar, y, constants)
    
    R = constants.R_gas_L_bar; 
    n_comps = constants.num_comps;
    
    a_i = zeros(1, n_comps);
    b_i = zeros(1, n_comps);
    for i = 1:n_comps
        Tc = constants.comps(i).Tc;
        Pc = constants.comps(i).Pc;
        Tr = T_K / Tc;
        
        a_i(i) = (0.42748 * R^2 * Tc^2.5) / Pc;
        b_i(i) = (0.08664 * R * Tc) / Pc;
    end
    
    a_m = (sum(y .* sqrt(a_i)))^2;
    b_m = sum(y .* b_i);
    
    A_m = (a_m * P_bar) / (R^2 * T_K^2.5);
    B_m = (b_m * P_bar) / (R * T_K);
    
    p = [1, -1, (A_m - B_m - B_m^2), (-A_m * B_m)];
    r = roots(p);
    
    Z = max(real(r(imag(r) == 0)));
    
    phi_hat = zeros(1, n_comps);
    for i = 1:n_comps
        term1 = (Z - 1) * b_i(i) / b_m;
        term2 = -log(Z - B_m);
        term3 = -(A_m / B_m) * (2 * sqrt(a_i(i) / a_m) - (b_i(i) / b_m)) * log(1 + (B_m / Z));
        
        ln_phi_hat_i = term1 + term2 + term3;
        phi_hat(i) = exp(ln_phi_hat_i);
    end
end
function phi_sat = calculate_phi_sat_RK(T_K, Psat_bar, i, constants)
    
    R = constants.R_gas_L_bar; 
    
    Tc = constants.comps(i).Tc;
    Pc = constants.comps(i).Pc;
    
    a_i = (0.42748 * R^2 * Tc^2.5) / Pc;
    b_i = (0.08664 * R * Tc) / Pc;
    
    A_i = (a_i * Psat_bar) / (R^2 * T_K^2.5);
    B_i = (b_i * Psat_bar) / (R * T_K);
    
    p = [1, -1, (A_i - B_i - B_i^2), (-A_i * B_i)];
    r = roots(p);
    
    Z = max(real(r(imag(r) == 0)));
    
    ln_phi_sat = (Z - 1) - log(Z - B_i) - (A_i / B_i) * log(1 + (B_i / Z));
    phi_sat = exp(ln_phi_sat);
end
function Psat_bar = calculate_Psat_Antoine(T_C, i, constants)
    
    if i == 2 || i == 3 
        Psat_bar = 1e9; 
        return;
    end
    
    ant = constants.comps(i).Antoine;
    if isempty(ant)
        Psat_bar = 1e9;
        return;
    end
    
    Psat_bar = 10^(ant(1) - ant(2) / (T_C + ant(3)));
end
function K = get_K_values_Raoult(T_C, P_MPa, constants)
    
    K = zeros(1, constants.num_comps);
    P_bar = P_MPa * 10;
    
    for i = 1:constants.num_comps
        comp = constants.comps(i);
        if i == 2 || i == 3 
            K(i) = 1e3; 
        else
            Psat_bar = calculate_Psat_Antoine(T_C, i, constants);
            K(i) = Psat_bar / P_bar;
        end
    end
end