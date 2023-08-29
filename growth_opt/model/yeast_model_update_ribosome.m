function [dydt, sig, met_reac, prot_syn_rate, beta, alpha, rib, tRNA, eIF_a_s, eIF_a_tau, other_met_reac, g_rate, ribo_rate] = yeast_model_update_ribosome(t, cell_state, par, mutant_type, snf1_vals, jgy_vals, mgl)

%% ------------------------------------------------------------------------
%                               Variables 
%--------------------------------------------------------------------------
%disp(t) 

% metabolites 
met   = cell_state(1:6); 
aa_ex = met(1); 
gl    = met(2); 
pc    = met(3); 
eh    = met(4); 
aa_in = met(5); 
ae    = met(6); 

% proteins
prot = cell_state(7:14);
p_r  = prot(1); %7
%{
p_z  = prot(2); %8
p_gy = prot(3); %9
p_fe = prot(4); %10
p_gn = prot(5); %11
p_mt = prot(6); %12
p_as = prot(7); %13
p_at = prot(8); %14
%}

% free ribosomes 
r0 = cell_state(15);

% number of cells 
n_cells = cell_state(16);

%% ------------------------------------------------------------------------
%                            Signaling network
%--------------------------------------------------------------------------

[sig] = module_signaling(met, par, mutant_type, snf1_vals);

%% ------------------------------------------------------------------------
%                            Metabolic network
%--------------------------------------------------------------------------

[met_reac] = module_met(met, prot, sig.snf1, par, mutant_type, jgy_vals, mgl); 

% fluxes
J_gy = met_reac.flux(1);
J_fe = met_reac.flux(2);
J_gn = met_reac.flux(3);
J_mt = met_reac.flux(4);
J_as = met_reac.flux(5);
J_at = met_reac.flux(6);

%% ------------------------------------------------------------------------
%                           Gene expression network
%--------------------------------------------------------------------------

[beta, alpha, rib, tRNA, eIF_a_s, eIF_a_tau] = module_gene_expression_update_ribosome(r0, met, prot, sig.snf1, sig.tor, par, mutant_type);

%% ------------------------------------------------------------------------
%                               growth rate
%--------------------------------------------------------------------------
 
prot_syn_rate = sum(par.l .* beta); % aa units 
g_rate = prot_syn_rate./par.pro_den;

%% ------------------------------------------------------------------------
%                              other fluxes
%--------------------------------------------------------------------------

J_po = 0;
J_eo = (par.q_eo * g_rate) + (par.d_em * ae) + (par.d_em * J_gn); 

other_met_reac.flux.po = J_po; 
other_met_reac.flux.eo = J_eo;  

%% ------------------------------------------------------------------------
%                              mass balances
%--------------------------------------------------------------------------
J_ATP = (par.q_gy * J_gy) ...
      + (par.q_mt * J_mt) ...
      - (par.q_gn * J_gn) ...
      - (par.q_as * J_as) ...
      - (par.q_p * prot_syn_rate) ...
      - (par.q_at * J_at) ...
      - J_eo ...
      - (g_rate * ae);

if strcmp(mutant_type,'none') ...
   || strcmp(mutant_type,'const_snf1') ...
   || strcmp(mutant_type,'const_tor') ...   
   || strcmp(mutant_type,'const_jgy') ...

    met_rate = [ (par.D * (par.aaex_in - aa_ex)) - ((par.v_cell * n_cells/par.v_ex) * J_at);                             % Aa_e:  extracellular amino acid                                                                                                % Ma
                 (par.D * (par.gl_in - gl)) - ((par.v_cell * n_cells/par.v_ex) * J_gy);                                  % Gl_ex: glucose                      
                 (par.n_gy * J_gy) + (par.n_gn * J_gn) - J_fe - J_mt - (par.n_as * J_as) - J_po - (g_rate * pc);         % Pc:    precursor 
                 (par.D * (par.eh_in - eh)) + (par.v_cell * n_cells/par.v_ex * (J_fe - J_gn));                           % Eh:    extracelluar EtOH                                       
                  par.q_aa * J_at + par.q_aa * J_as - par.q_aa * prot_syn_rate  - (g_rate * aa_in);                      % Aa_i:  amino acid; Jao = 0
                  J_ATP;                                                                                                 % Ae:    ATP
               ];       
           
elseif strcmp(mutant_type,'const_snf1_gl') ...
       || strcmp(mutant_type,'const_tor_gl') ...
       || strcmp(mutant_type,'const_gl') ...
       || strcmp(mutant_type,'const_jgy_gl')
   
    met_rate = [ (par.D * (par.aaex_in - aa_ex)) - ((par.v_cell * n_cells/par.v_ex) * J_at);                             % Aa_e:  extracellular amino acid                                                                                                % Ma
                  0;                                                                                                     % Gl_ex: glucose                      
                 (par.n_gy * J_gy) + (par.n_gn * J_gn) - J_fe - J_mt - (par.n_as * J_as)  - J_po - (g_rate * pc);        % Pc:    precursor
                 (par.D * (par.eh_in - eh)) + (par.v_cell * n_cells/par.v_ex * (J_fe - J_gn));                           % Eh:    extracelluar EtOH                                       
                  par.q_aa * J_at + par.q_aa * J_as - par.q_aa * prot_syn_rate  - (g_rate * aa_in);                       % Aa_i: amino acid
                  J_ATP;                                                                                                 % Ae:    ATP
               ];      

elseif strcmp(mutant_type,'const_snf1_eh') ...
       || strcmp(mutant_type,'const_tor_eh') ...
       || strcmp(mutant_type,'const_eh') ...
       || strcmp(mutant_type,'const_jgy_eh')

    met_rate = [ (par.D * (par.aaex_in - aa_ex)) - ((par.v_cell * n_cells/par.v_ex) * J_at);                             % Aa_e:  extracellular amino acid                                                                                                % Ma
                 (par.D * (par.gl_in - gl)) - ((par.v_cell * n_cells/par.v_ex) * J_gy);                                  % Gl_ex: glucose                      
                 (par.n_gy * J_gy) + (par.n_gn * J_gn) - J_fe - J_mt - (par.n_as * J_as) - J_po - (g_rate * pc);         % Pc:    precursor 
                  0;                                                                                                     % Eh:    extracelluar EtOH                                       
                  par.q_aa * J_at + par.q_aa * J_as - par.q_aa * prot_syn_rate  - (g_rate * aa_in);                      % Aa_i:  amino acid; Jao = 0
                  J_ATP;                                                                                                 % Ae:    ATP
               ];         

elseif strcmp(mutant_type,'const_gl_eh_aaex_cell') ...
       || strcmp(mutant_type,'const_snf1_gl_eh_aaex_cell')...
       || strcmp(mutant_type,'const_tor_gl_eh_aaex_cell')...
   
    met_rate = [ 0;                                                                                                      % Aa_e:  extracellular amino acid                                                                                                % Ma
                 0;                                                                                                      % Gl_ex: glucose                      
                (par.n_gy * J_gy) + (par.n_gn * J_gn) - J_fe - J_mt - (par.n_as * J_as) - J_po - (g_rate * pc);          % Pc:    precursor
                 0;                                                                                                      % Eh:   extracelluar EtOH                                       
                 par.q_aa * J_at + par.q_aa * J_as - par.q_aa * prot_syn_rate - (g_rate * aa_in);                        % Aa_i: amino acid
                 J_ATP;                                                                                                  % Ae:   ATP
                ];                    
end 

prot_rate = beta - (g_rate * prot); % protein units 
ribo_rate = par.k_ro * (p_r - (r0 * par.n_ro)) - (par.d_ro * r0) - (g_rate * r0); % par.n_ro = 1

if strcmp(mutant_type,'none') ...
        || strcmp(mutant_type,'const_snf1') ...
        || strcmp(mutant_type,'const_snf1_gl') ...
        || strcmp(mutant_type,'const_tor') ...
        || strcmp(mutant_type,'const_tor_gl') ...        
        || strcmp(mutant_type,'const_gl') ...
        || strcmp(mutant_type,'const_jgy') ...
        || strcmp(mutant_type,'const_snf1_eh')...
        || strcmp(mutant_type,'const_tor_eh') ...
        || strcmp(mutant_type,'const_eh') ...
        || strcmp(mutant_type,'const_jgy_eh')     

    cell_rate = - (par.D * n_cells) + (g_rate * n_cells);
    
elseif strcmp(mutant_type,'const_gl_eh_aaex_cell') ...
        || strcmp(mutant_type,'const_snf1_gl_eh_aaex_cell')...
        || strcmp(mutant_type,'const_tor_gl_eh_aaex_cell')...
    cell_rate = 0; 
    
end 

%% ------------------------------------------------------------------------
%                                  outputs
%--------------------------------------------------------------------------

dydt = [met_rate;
        prot_rate;
        ribo_rate;
        cell_rate];

end 