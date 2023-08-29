
clc; close all; clear

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

 
%% Fig 2: batch fermentation with wt 
tic

% read parameters
par = read_parametersv2;

% wt ----------------------------------------------------------------------
mutant = 'none';

% start and end times 
t_batch_start = 0; 
t_batch_final = 20; 

% get reasonable initial cell state
cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_wt, jgy_vals, mgl_0);
disp('Calculate initial conditions for  batch simulation: done') 

% run batch fermentation 
disp('Running  batch simulation...')                 
[t_wt, y_wt]  = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_wt, jgy_vals, mgl_0), ...
                        [t_batch_start t_batch_final], ...
                         cell_state, ...
                         options);                                                           

% snf1 mutant -------------------------------------------------------------

mutant = 'const_snf1';

for m = 1:numel(snf1_vals)
    disp(m)
    % get specific snf1 structure name 
    s_val = sprintf('snf1_%g',snf1_vals(m));
    s_val(s_val  == '.') = '_'; 
    
    % get cell state
    cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals(m), jgy_vals, mgl_0);
    
    % run batch fermentation 
    [t, y]  = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals(m), jgy_vals, mgl_0), ...
                            [t_batch_start t_batch_final], ...
                             cell_state, ...
                             options);  
                         
    [t10_s, y10_s]  = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals(m), jgy_vals, mgl_0), ...
                            [t_batch_start t_batch_log], ...
                             cell_state, ...
                             options);  
                         
    y10(m,:) = y10_s(end,:) ;                    
    t_batch_10.(s_val) = t;                   
                             
    y = real(y);   
    y_batch_10.(s_val) = y; 
    
    % initialize arrays to hold intermediate values 
    met_reac_10.(s_val).prot      = ones(numel(t_batch_10.(s_val)), 6); 
    met_reac_10.(s_val).substrate = ones(numel(t_batch_10.(s_val)), 6); 
    met_reac_10.(s_val).atp       = ones(numel(t_batch_10.(s_val)), 6); 
    met_reac_10.(s_val).sig       = ones(numel(t_batch_10.(s_val)), 6); 
    met_reac_10.(s_val).flux      = ones(numel(t_batch_10.(s_val)), 6); 

    prot_syn_10.(s_val).alpha     = ones(numel(t_batch_10.(s_val)), 8); 
    prot_syn_10.(s_val).tc        = ones(numel(t_batch_10.(s_val)), 1);  
    prot_syn_10.(s_val).eIF_a     = ones(numel(t_batch_10.(s_val)), 1);  

    sig_10.(s_val).snf1 = ones(numel(t_batch_10.(s_val)), 1);
    sig_10.(s_val).tor  = ones(numel(t_batch_10.(s_val)), 1); 

    % get intermediate values 
    for k = 1:length(t_batch_10.(s_val))
       [~, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, rib_t, tRNA_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t, g_rate_t, ribo_rate_t] = yeast_model_update_ribosome(t_batch_10.(s_val)(k), y_batch_10.(s_val)(k,:)', par, mutant, snf1_vals(m), jgy_vals, mgl_0);

        met_reac_10.(s_val).prot(k,:)      = real(met_reac_t.prot)';
        met_reac_10.(s_val).substrate(k,:) = real(met_reac_t.substrate)';
        met_reac_10.(s_val).atp(k,:)       = real(met_reac_t.atp)';
        met_reac_10.(s_val).sig(k,:)       = real(met_reac_t.sig)';
        met_reac_10.(s_val).flux(k,:)      = real(met_reac_t.flux)';

        prot_syn_10.(s_val).alpha(k,:) = real(table2array(struct2table(alpha_t))); 
        prot_syn_10.(s_val).tc(k,:)    = real(tRNA_t.tc)'; 

        sig_10.(s_val).snf1(k,:) = real(sig_t.snf1)';
        sig_10.(s_val).tor(k,:)  = real(sig_t.tor)'; 
    end 
end 
                         
toc

%% batch fermentation with constant snf1 (figure 1, S1) - 25 points 
%
tic

% read parameters
par = read_parametersv2;

% wt ----------------------------------------------------------------------

mutant        = 'none';

% start and end times 
t_batch_start = 0; 
t_batch_final = 20; 

% get reasonable initial cell state
cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_wt, jgy_vals, mgl_0);
disp('Calculate initial conditions for  batch simulation: done') 

% run batch fermentation 
disp('Running  batch simulation...') 
[t_wt_25, y_wt_25]  = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_wt, jgy_vals, mgl_0), ...
                        [t_batch_start t_batch_final], ...
                         cell_state, ...
                         options);                                                           

% snf1 mutant -------------------------------------------------------------

mutant = 'const_snf1';

for m = 1:numel(snf1_vals_25)
    disp(m)
    % get specific snf1 structure name 
    s_val = sprintf('snf1_%g',snf1_vals_25(m));
    s_val(s_val  == '.') = '_'; 
    
    % get reasonable initial cell state
    cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_25(m), jgy_vals, mgl_0);
    
    % run batch fermentation 
    [t25_s, y25_s]  = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_25(m), jgy_vals, mgl_0), ...
                        [t_batch_start t_batch_log], ...
                         cell_state, ...
                         options);  

    y25(m,:) = y25_s(end,:) ;                    
    t_batch_25.(s_val) = t;                   
                             
end 
                         
toc
disp('Batch fermentation: done')


%% Figure 2 and related SI plotting ---------------------------------------
close all 
disp('Figure 2')

%simu_dots     = 1;
%sumu_dots_num = 9;

% plot format 
x_lim_time   = [0 15];




% Sim: cells vs time
snf1_val = flip(snf1_vals); % flip snf1 from increasing order to decreasing order
                            % to match 

hold on; 
for k = 1:numel(snf1_val)
    s_val = sprintf('snf1_%g',snf1_val(k));
    s_val(s_val  == '.') = '_'; 
    plot(t_batch_10.(s_val), y_batch_10.(s_val)(:,end),'Color',camp_colors(k,:)) 
end 

if length(y_wt) ~= 1
    plot(t_wt, y_wt(:,end),'--k')
end

hold off;
ylabel('cells')
xlabel(x_label_batch)
xlim(x_lim_time)
xticks([0 7 14])
xlim([0 14])
ylim([0, 1.1*max(ylim)])
axis square; 
box on; 

