
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
%x_lim_time   = [0 15];
x_lim_camp   = [0.5 1]; 
x_ticks_camp = [0.5 0.7 1];

% Figure 2 ---------------------------------------------------------------
figure; 
a = 2; b = 1; 

% % Simulation: gr vs cAMP 
% % select snf1 mutants to plot 
x_camp         = snf1_to_camp(snf1_vals_25);
int_index      = [1, 3, 7, 11, 14, 17, 20, 23, 25];


% Sim: gl, eh, cells vs cAMP
[minValue,closestIndex] = min(abs(snf1_to_camp(snf1_vals_25)-snf1_to_camp(sig_const_gl.snf1))); % find the snf1 mutant that most resembles the wt 

subplot(a, b, 1)
hold on;
plot(snf1_to_camp(snf1_vals_25), y25(:,num_y.gl), '-', 'Color', glu_color) % glucose, entire simulation 
plot(snf1_to_camp(snf1_vals_25), y25(:,num_y.eh), '-', 'Color', eth_color) % ethanol, entire simulation 

plot(x_camp(int_index), y25(int_index,num_y.gl), 'o', 'MarkerFaceColor', glu_color) % glucose, selected points
plot(x_camp(int_index), y25(int_index,num_y.eh), 'o', 'MarkerFaceColor', eth_color) % ethanol, selected points 

plot(x_camp(int_index(3)), y25(int_index(3),num_y.gl), 'o', 'MarkerFaceColor', wt_camp_gl) % use darker color to indicate lowest mutant residual glucose concentration
plot(x_camp(int_index(3)), y25(int_index(3),num_y.eh), 'o', 'MarkerFaceColor', wt_camp_eh) % use darker color to indicate highest mutant ethanol concentration 

yline(y25(closestIndex,num_y.gl), '--', 'Color', wt_camp_gl) % wt residual glucose concentration 
yline(y25(closestIndex,num_y.eh), '--', 'Color', wt_camp_eh) % wt ethanol concentration 

hold off;
ylabel({'glucose (uM)','ethanol (uM)'})
xlabel('cAMP equiv.')  
ylim([0 1.2*max(ylim)])
xlim(x_lim_camp)
xticks(x_ticks_camp)
set(gca, 'XScale', 'log')
set(gca,'XMinorTick','off')
axis square; 
box on; 



%  legend 
sb = subplot(a, b, 1+b);
set(sb, 'Visible', 'off')
yyaxis left
hold on 
plot(0,  0, '-o', 'Color', glu_color, 'MarkerFaceColor', glu_color)
plot(0,  0, '-o', 'Color', eth_color, 'MarkerFaceColor', eth_color)
hold off 
legend({'glucose', 'ethanol'})


disp("Fig.2: done")


