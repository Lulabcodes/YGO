
clc; close all; 

addpath('../../')
addpath(genpath('../../'))

% load simulation set, data
shared_setup;

%% Simulation and plotting for Figure 4 related SI 
%{\
disp('Figure 4')
% Simulation for closed- and open-loop experiments with different methylglyoxyal 
tic 

% initialize arrays for wt constant glu simulation
mgl_selec_5 = [0 0.1 3 6 30]* 10^3;
mgl_selec_5_wt_y_const_gl                  = ones(numel(mgl_selec_5), numel(fieldnames(num_y)));
mgl_selec_5_wt_met_reac_const_gl.prot      = ones(numel(mgl_selec_5), 6); 
mgl_selec_5_wt_met_reac_const_gl.substrate = ones(numel(mgl_selec_5), 6);
mgl_selec_5_wt_met_reac_const_gl.atp       = ones(numel(mgl_selec_5), 6); 
mgl_selec_5_wt_met_reac_const_gl.sig       = ones(numel(mgl_selec_5), 6);
mgl_selec_5_wt_met_reac_const_gl.flux      = ones(numel(mgl_selec_5), 6); 
mgl_selec_5_wt_prot_syn_const_gl.alpha     = ones(numel(mgl_selec_5), 8);
mgl_selec_5_wt_prot_syn_const_gl.tc        = ones(numel(mgl_selec_5), 1);  
mgl_selec_5_wt_prot_syn_const_gl.eIF_a     = ones(numel(mgl_selec_5), 1);  
mgl_selec_5_wt_sig_const_gl.snf1           = ones(numel(mgl_selec_5), 1);
mgl_selec_5_wt_sig_const_gl.tor            = ones(numel(mgl_selec_5), 1); 
mgl_selec_5_wt_g_rate_const_gl             = ones(numel(mgl_selec_5), 1);
    
% initialize arrays for snf1 mutant constant glu simulation 
mgl_selec_5_y_snf1_s_val                  = ones(numel(mgl_selec_5), numel(snf1_vals_25), numel(fieldnames(num_y))); 
mgl_selec_5_met_reac_snf1_s_val.prot      = ones(numel(mgl_selec_5), numel(snf1_vals_25), 6); 
mgl_selec_5_met_reac_snf1_s_val.substrate = ones(numel(mgl_selec_5), numel(snf1_vals_25), 6); 
mgl_selec_5_met_reac_snf1_s_val.atp       = ones(numel(mgl_selec_5), numel(snf1_vals_25), 6);
mgl_selec_5_met_reac_snf1_s_val.sig       = ones(numel(mgl_selec_5), numel(snf1_vals_25), 6); 
mgl_selec_5_met_reac_snf1_s_val.flux      = ones(numel(mgl_selec_5), numel(snf1_vals_25), 6); 
mgl_selec_5_prot_syn_snf1_s_val.alpha     = ones(numel(mgl_selec_5), numel(snf1_vals_25), 8); 
mgl_selec_5_prot_syn_snf1_s_val.tc        = ones(numel(mgl_selec_5), numel(snf1_vals_25), 1);  
mgl_selec_5_sig_snf1_s_val.snf1           = ones(numel(mgl_selec_5), numel(snf1_vals_25), 1);
mgl_selec_5_sig_snf1_s_val.tor            = ones(numel(mgl_selec_5), numel(snf1_vals_25), 1); 
mgl_selec_5_g_rate_snf1_s_val             = ones(numel(mgl_selec_5), numel(snf1_vals_25), 1);

% run simulation for different methylglyoxal concentrations ---------------
for i = 1: numel(mgl_selec_5)

    % WT simulation with constant extracellular glucose concentration -----
    mutant = 'const_gl'; % glucose concentration is held constant
    
    % initial conditions 
    met.ex_amino_acids = 2.3*10^5;                  % 1 % rich media 
    met.glucose        = 20*5550.7;                 % 2  
    met.precursor      = 1.00E+03;                  % 3  
    met.ethanol        = 0; 	                    % 4  
    met.in_amino_acids = 1.50E+05;                  % 5        
    met.atp            = 1.00E+03;                  % 6 
    num.met = numel(fieldnames(met));
    
    cell_state = [table2array(struct2table(met))';  % 1 - 7
                  table2array(struct2table(prot))'; % 8 - 17 
                  R0;                               % 18
                  cells;                            % 19
                  ];
              
    % ODE options           
    Rel_tol  = 1.0E-03; 
    Abs_tol  = 1.0E-06; 
    options  = odeset('RelTol',Rel_tol, ...
                      'AbsTol',Abs_tol, ...
                      'NonNegative',(1:length(cell_state)));
    
    t_batch_start = 0; 
    t_batch_final = 1;   
     
    % find reasonable initial cell state 
    cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_wt, jgy_vals, mgl_selec_5(i));
    disp('Calculate initial conditions for const_gl batch simulation: done') 
    
    % run simulation 
    disp('Running const_gl batch simulation...') 
    [t, y]  = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_wt, jgy_vals, mgl_selec_5(i)), ...
                                [t_batch_start t_batch_final], ...
                                 cell_state, ...
                                 options); 
    
    mgl_selec_5_wt_y          = real(y);           
    mgl_selec_5_wt_y_const_gl(i,:) = y(end,:);
                             
    % get intermediate values   
    [~, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, rib_t, tRNA_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t, g_rate_t, ribo_rate_t] = yeast_model_update_ribosome(t(end),  y(end,:)', par, mutant, snf1_vals_wt, jgy_vals, mgl_selec_5(i));
        
    mgl_selec_5_wt_met_reac_const_gl.prot(i,:)      = real(met_reac_t.prot)';
    mgl_selec_5_wt_met_reac_const_gl.substrate(i,:) = real(met_reac_t.substrate)';
    mgl_selec_5_wt_met_reac_const_gl.atp(i,:)       = real(met_reac_t.atp)';
    mgl_selec_5_wt_met_reac_const_gl.sig(i,:)       = real(met_reac_t.sig)';
    mgl_selec_5_wt_met_reac_const_gl.flux(i,:)      = real(met_reac_t.flux)';
        
    mgl_selec_5_wt_prot_syn_const_gl.alpha(i,:) = real(table2array(struct2table(alpha_t))); 
    mgl_selec_5_wt_prot_syn_const_gl.tc(i,:)    = real(tRNA_t.tc)'; 
        
    mgl_selec_5_wt_sig_const_gl.snf1(i,:) = real(sig_t.snf1)';
    mgl_selec_5_wt_sig_const_gl.tor(i,:)  = real(sig_t.tor)'; 
        
    mgl_selec_5_wt_g_rate_const_gl(i,:) = real(g_rate_t)';
    
    disp('WT simulation with constant glucose: done')
    
    
    

    % PKA/Snf1 mutant simulation with constant extracellular glucose
    % concentration -------------------------------------------------------
    
    mutant = 'const_snf1_gl'; % s = const and glucose concentration is held constant 

    % initial conditions 
    met.ex_amino_acids = 2.3*10^5;                  % 1 %minimal media = 0, YPD = some value 
    met.glucose        = 20*5550.7;                 % 2  
    met.precursor      = 1.00E+03;                  % 3  
    met.ethanol        = 0;                         % 4  
    met.in_amino_acids = 1.50E+05;                  % 5                
    met.atp            = 1.00E+03;                  % 6
    num.met = numel(fieldnames(met));
    
    cell_state = [table2array(struct2table(met))';  % 1 - 7
                  table2array(struct2table(prot))'; % 8 - 17 
                  R0;                               % 18
                  cells;                            % 19
                  ];    
    
    
    t_batch_start = 0; 
    t_batch_final = 1;   
    
    % read parameters
    par = read_parametersv2;
    
    for m = 1:numel(snf1_vals_25)
        disp(m)
      
        % get reasonable initial cell state  
        cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_25(m), jgy_vals, mgl_selec_5(i));
        disp('Calculate initial conditions for const_gl batch simulation: done') 
    
        % run simulation 
        disp('Running const_gl batch simulation...') 
        [t, y]  = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_25(m), jgy_vals, mgl_selec_5(i)), ...
                                    [t_batch_start t_batch_final], ...
                                     cell_state, ...
                                     options); 
                                            
        y = real(y);  
        mgl_selec_5_y_snf1_s_val(i,m,:) = y(end,:);
        
        % get intermediate values 
        [~, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, rib_t, tRNA_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t, g_rate_t, ribo_rate_t] = yeast_model_update_ribosome(t(end), y(end,:)', par, mutant, snf1_vals_25(m), jgy_vals, mgl_selec_5(i));
        mgl_selec_5_met_reac_snf1_s_val.prot(i,m,:)      = real(met_reac_t.prot)';
        mgl_selec_5_met_reac_snf1_s_val.substrate(i,m,:) = real(met_reac_t.substrate)';
        mgl_selec_5_met_reac_snf1_s_val.atp(i,m,:)       = real(met_reac_t.atp)';
        mgl_selec_5_met_reac_snf1_s_val.sig(i,m,:)       = real(met_reac_t.sig)';
        mgl_selec_5_met_reac_snf1_s_val.flux(i,m,:)      = real(met_reac_t.flux)';
    
        mgl_selec_5_prot_syn_snf1_s_val.alpha(i,m,:) = real(table2array(struct2table(alpha_t))); 
        mgl_selec_5_prot_syn_snf1_s_val.tc(i,m,:)    = real(tRNA_t.tc)'; 
    
        mgl_selec_5_sig_snf1_s_val.snf1(i,m,:) = real(sig_t.snf1)';
        mgl_selec_5_sig_snf1_s_val.tor(i,m,:)  = real(sig_t.tor)'; 
            
        mgl_selec_5_g_rate_snf1_s_val(i,m,:) = g_rate_t; 
        
    end 
    
    toc  
    disp('Constant Snf1 and with mgl constant glucose simulation: done')
    
end 

%% plotting ---------------------------------------------------------------

% Figure SI 7: simulated growth rate vs snf1 (camp equiv.) ----------------
figure; 
for i = 1: numel(mgl_selec_5)
    subplot(ceil(numel(mgl_selec_5)/5), 5, i)
    hold on;
    
    plot(snf1_to_camp(snf1_vals_25), mgl_selec_5_g_rate_snf1_s_val(i,:,:), 'Color', plt_clrs.lightgray)
    
    [~,indxpeak(i)] = max(mgl_selec_5_g_rate_snf1_s_val(i,:,:)); % find peak
    plot(snf1_to_camp(snf1_vals_25(indxpeak(i))), mgl_selec_5_g_rate_snf1_s_val(i,indxpeak(i),:), 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 10) % use darker color for largest growth rate
    
    yline(mgl_selec_5_wt_g_rate_const_gl(i),'--') % wt 
    
    hold off;
    xlabel('model cAMP eq.')
    ylabel('growth rate (h^{-1})')
    ylim([0 0.6])
    set(gca, 'XScale', 'log')
    xlim([0.5 1])
    xticks([0.5 0.7 1])
    axis square;
    box on;
end 

%}
