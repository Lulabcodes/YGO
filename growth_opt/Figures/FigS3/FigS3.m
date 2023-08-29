
clc; close all; 

addpath('../../')
addpath(genpath('../../'))

% load simulation setup, data
shared_setup;
%% Figure 2 - constant glucose snf1 simulation 25 points (wt)
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
 

%% Plotting ---------------------------------------------------------------

% Supplementary Figure 2 --------------------------------------------------

y_ticks = [10^-5 10^-3 10^-1]; % y_ticks for E fractions.

% convert y_snf1 to y_snf1_mat 
% y_snf1: a struct of model variable length(snf1_vals_25) number of fields
% y_snf1_mat: a matrix: [num_y by length(snf1_vals_25)]
y_snf1_size  = [length(fieldnames(num_y)),length(snf1_vals_25)];
y_snf1_mat   = reshape(table2array(struct2table(y_snf1)),y_snf1_size)';
[~, indx_WT] = min(abs(snf1_vals_25-sig_const_gl.snf1)); 

% calculate R, E, Z fraction:
total_protein_con = sum(y_snf1_mat(:,num_y.r:num_y.gn).* par.l(num_prot.r:num_prot.gn)',2);
R_snf1_25   = y_snf1_mat(:,num_y.r).*par.l(num_prot.r)'./total_protein_con(1);
Z_snf1_25   = y_snf1_mat(:,num_y.z).*par.l(num_prot.z)'./total_protein_con(1);
E_snf1_25   = sum(y_snf1_mat(:,num_y.gy:num_y.at).*par.l(num_prot.gy:num_prot.at)',2)./par.pro_den;
Egy_snf1_25 = sum(y_snf1_mat(:,num_y.gy).*par.l(num_prot.gy)',2)./par.pro_den;
Efe_snf1_25 = sum(y_snf1_mat(:,num_y.fe).*par.l(num_prot.fe)',2)./par.pro_den;
Egn_snf1_25 = sum(y_snf1_mat(:,num_y.gn).*par.l(num_prot.gn)',2)./par.pro_den;
Emt_snf1_25 = sum(y_snf1_mat(:,num_y.mt).*par.l(num_prot.mt)',2)./par.pro_den;
Eas_snf1_25 = sum(y_snf1_mat(:,num_y.as).*par.l(num_prot.as)',2)./par.pro_den;
Eat_snf1_25 = sum(y_snf1_mat(:,num_y.at).*par.l(num_prot.at)',2)./par.pro_den;
E_all = [Egy_snf1_25 Efe_snf1_25 Egn_snf1_25 Emt_snf1_25 Eas_snf1_25 Eat_snf1_25]; 

figure;
a = 1; b = 2; 

% sim: REZ vs cAMP 
subplot(a, b, 1)
hold on;
plot(snf1_to_camp(snf1_vals_25), R_snf1_25, '-', 'Color', plt_clrs.yellow) % R
plot(snf1_to_camp(snf1_vals_25), E_snf1_25, '-', 'Color', plt_clrs.green)  % E
plot(snf1_to_camp(snf1_vals_25), Z_snf1_25, '-', 'Color', plt_clrs.gray)   % Z
xline(snf1_to_camp(sig_const_gl.snf1),'--') % wt snf1 (camp equiv.) value 
hold off;
xlabel('model cAMP equiv.')
ylabel({'proteins fraction'})
xlim([0.4 1])
xticks([0.4 0.6 0.9])
ylim([0 1.05*max(ylim)])
set(gca, 'XScale', 'log')
set(gca,'XMinorTick','off')
legend('R ','E ','Z','WT cAMP')
legend box off;
axis square; 
box on; 

% sim: plot Emt vs cAMP 
subplot(a, b, 2)
hold on 
plot(snf1_to_camp(snf1_vals_25), Emt_snf1_25, '-', 'Color', plt_clrs.green)
xline(snf1_to_camp(sig_const_gl.snf1), '--') % wt snf1 (cAMP equiv.) value
hold off
xlabel('model cAMP equiv.')
ylabel({'E_{mt} fraction'})
xlim([0.4 1])
xticks([0.4 0.6 0.9])
ylim([0 1.05*max(ylim)])
set(gca, 'XScale', 'log')
set(gca,'XMinorTick','off')
axis square; 
box on; 

disp("SI2: done")



%}
