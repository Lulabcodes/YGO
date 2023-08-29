
clc; close all; clear;

addpath('../../')
addpath(genpath('../../'))

% load simulation setup, data
shared_setup;

%% Figure 2 - constant glucose snf1 simulation (wt)
%
tic
disp('Figure 2') 

% read parameters
par = read_parametersv2;

mutant = 'const_gl'; % glucose concentration is held constant

% start and end times 
t_batch_start = 0; 
t_batch_final = 1; % can choose small number because initial state will already be the steady-state

% find reasonable initial cell state (what the steady-state for the given initial conditions is)
cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_wt, jgy_vals, mgl_0);
disp('Calculate initial conditions for const_gl batch simulation: done') 

% run simulation 
disp('Running const_gl batch simulation...') 
[t, y]  = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_wt, jgy_vals, mgl_0), ...
                            [t_batch_start t_batch_final], ...
                             cell_state, ...
                             options); 

% save simulation values 
y = real(y);           
y_const_gl = y(end,:);

% initialize arrays to hold intermediate values 
met_reac_const_gl.prot      = ones(1, 6); 
met_reac_const_gl.substrate = ones(1, 6);
met_reac_const_gl.atp       = ones(1, 6); 
met_reac_const_gl.sig       = ones(1, 6);
met_reac_const_gl.flux      = ones(1, 6); 

prot_syn_const_gl.alpha     = ones(1, 8);
prot_syn_const_gl.tc        = 1;  
prot_syn_const_gl.eIF_a     = 1;  

sig_const_gl.snf1           = 1;
sig_const_gl.tor            = 1; 

%g_rate_const_gl            = 1;

% get intermediate values   
[~, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, rib_t, tRNA_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t, g_rate_t, ribo_rate_t] = yeast_model_update_ribosome(t(end), y_const_gl', par, mutant, snf1_vals_wt, jgy_vals, mgl_0);
    
met_reac_const_gl.prot(:)      = real(met_reac_t.prot)';
met_reac_const_gl.substrate(:) = real(met_reac_t.substrate)';
met_reac_const_gl.atp(:)       = real(met_reac_t.atp)';
met_reac_const_gl.sig(:)       = real(met_reac_t.sig)';
met_reac_const_gl.flux(:)      = real(met_reac_t.flux)';
    
prot_syn_const_gl.alpha(:) = real(table2array(struct2table(alpha_t))); 
prot_syn_const_gl.tc(:)    = real(tRNA_t.tc)'; 
    
sig_const_gl.snf1(:) = real(sig_t.snf1)';
sig_const_gl.tor(:)  = real(sig_t.tor)'; 
    
g_rate_const_gl = real(g_rate_t)';
 
toc 
disp('WT simulation with constant glucose: done')

%% Fig 2: snf1 mutant simulation with constant extracellular glucose concentration
tic 

% read parameters
par = read_parametersv2;

mutant = 'const_snf1_gl';

% % start and end times 
t_batch_start = 0; 
t_batch_final = 1; % can choose small number because initial state will already be the steady-state

y_snf1 = [];

for m = 1:numel(snf1_vals_25)
    disp(m)
    % get specific snf1 structure name 
    s_val = sprintf('snf1_%g',snf1_vals_25(m));
    s_val(s_val  == '.') = '_'; 
    
    % find reasonable initial cell state (what the steady-state for the given initial conditions is)
    cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_25(m), jgy_vals, mgl_0);
    disp('Calculate initial conditions for const_gl batch simulation: done') 
    
    % run simulation 
    disp('Running const_gl batch simulation...') 
    [t, y]  = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_25(m), jgy_vals, mgl_0), ...
                                [t_batch_start t_batch_final], ...
                                 cell_state, ...
                                 options);            
                             
    y = real(y);  
    y_snf1.(s_val) = y(end,:);
    
    % initialize arrays to hold intermediate values 
    met_reac_snf1.(s_val).prot      = ones(1, 6); 
    met_reac_snf1.(s_val).substrate = ones(1, 6); 
    met_reac_snf1.(s_val).atp       = ones(1, 6);
    met_reac_snf1.(s_val).sig       = ones(1, 6); 
    met_reac_snf1.(s_val).flux      = ones(1, 6); 

    prot_syn_snf1.(s_val).alpha     = ones(1, 8); 
    prot_syn_snf1.(s_val).tc        = 1;  
    
    sig_snf1.(s_val).snf1 = 1;
    sig_snf1.(s_val).tor  = 1; 
    
    g_rate_snf1.(s_val) = 1;

    % get intermediate values 
    [~, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, rib_t, tRNA_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t, g_rate_t, ribo_rate_t] = yeast_model_update_ribosome(t(end), y_snf1.(s_val)', par, mutant, snf1_vals_25(m), jgy_vals, mgl_0);
    met_reac_snf1.(s_val).prot(:)      = real(met_reac_t.prot)';
    met_reac_snf1.(s_val).substrate(:) = real(met_reac_t.substrate)';
    met_reac_snf1.(s_val).atp(:)       = real(met_reac_t.atp)';
    met_reac_snf1.(s_val).sig(:)       = real(met_reac_t.sig)';
    met_reac_snf1.(s_val).flux(:)      = real(met_reac_t.flux)';

    prot_syn_snf1.(s_val).alpha(:) = real(table2array(struct2table(alpha_t))); 
    prot_syn_snf1.(s_val).tc(:)    = real(tRNA_t.tc)'; 
    prot_syn_snf1.(s_val).rib.others(:) = real(rib_t.others)';
    prot_syn_snf1.(s_val).rib.others(:) = real(rib_t.others)';
        
    sig_snf1.(s_val).snf1(:) = real(sig_t.snf1)';
    sig_snf1.(s_val).tor(:)  = real(sig_t.tor)'; 
        
    g_rate_snf1.(s_val) = g_rate_t; 
        
end 

toc  
disp('Constant Snf1 and constant glucose simulation: done')
 

%% Figure 2 and related SI plotting ---------------------------------------
close all 
disp('Figure 2')

%simu_dots     = 1;
%sumu_dots_num = 9;

% plot format 
x_lim_camp   = [0.5 1]; 
x_ticks_camp = [0.5 0.7 1];

% Figure 2 ---------------------------------------------------------------
figure; 
a = 1; b = 2; 

% Simulation: gr vs cAMP 
% select snf1 mutants to plot 
x_camp         = snf1_to_camp(snf1_vals_25);
[~,start_indx] = min(abs(x_camp-x_lim_camp(1)));
[~,end_indx]   = min(abs(x_camp-x_lim_camp(2)));
int_index      = [1, 3, 7, 11, 14, 17, 20, 23, 25];

gr = table2array(struct2table(g_rate_snf1));
[~, max_idx] = max(gr(int_index));

subplot(a, b, 1)
hold on;
plot(snf1_to_camp(snf1_vals_25),table2array(struct2table(g_rate_snf1)), '-', 'Color', camp_gr) % plot snf1 mutant growth rates vs snf1 (camp equiv.), entire simulation 
plot(x_camp(int_index),         gr(int_index),         'o','MarkerFaceColor', camp_gr) % plot snf1 mutant growth rates vs snf1 (camp equiv.), selected points
plot(x_camp(int_index(max_idx)),gr(int_index(max_idx)),'o','MarkerFaceColor', 'k') % make snf1 mutant max growth rate black

yline(g_rate_const_gl, '--') % wt growth rate 

hold off;
ylabel('growth rate (h^{-1})')
xlabel('cAMP equiv.')
ylim([0 1.15*max(ylim)])
xlim(x_lim_camp)
xticks(x_ticks_camp)
set(gca, 'XScale', 'log')
set(gca,'XMinorTick','off')
box on; 
axis square


%  legend 
sb = subplot(a, b, 2);
set(sb, 'Visible', 'off')
yyaxis left
hold on; 
plot(0, 0, '-o', 'MarkerFaceColor', camp_gr, 'Color', camp_gr) 
plot(0, 0, 'o',  'MarkerFaceColor', 'k',     'Color', 'k') 
plot(0, 0,  'k--')
hold off; 
legend({'cEAC', 'cEAC max.', 'cIAC'})








