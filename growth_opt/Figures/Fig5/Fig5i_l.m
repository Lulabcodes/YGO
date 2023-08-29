clc; close all; clear; 

addpath('../../')
addpath(genpath('../../'))

% load simulation set, data
shared_setup;

%% Simulation and plotting for Figure 5: YPE --> YPD and related SI
%{\
disp('Figure 5: YPE --> YPD')
tic

%cells  = cell_mason_unit_conv; 

% batch wt fermentation experiment start and end times     
t_batch_eh_start = 0; 
t_batch_eh_final = 24; 

t_batch_gl_start = 0; 
t_batch_gl_final = 32 - 24; 4; 

% simulations 
% batch fermentation with wt, high and low constant snf1 on ethanol -------

% wt 
disp('WT batch fermentation - ethanol')

mutant = 'none';

% initial conditions 
met.ex_amino_acids = 2.3*10^6;               % 1 %minimal media = 0, YPD = some value 
met.glucose        = 0;                      % 2  
%met.precursor      = 1.00E+03;              % 3  
met.ethanol        = 20 * 21706.0994139353;  % 4  
%met.in_amino_acids = 1.50E+05;              % 5  
%met.atp            = 1.00E+03;              % 6 
num.met = numel(fieldnames(met));

cell_state = [table2array(struct2table(met))';  
              table2array(struct2table(prot))'; 
              R0;                               
              cells;                         
              ];

% get a reasonable initial cell state 
cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_wt, jgy_vals, mgl_0);
disp('Calculated initial conditions: done')

% run simulation of batch fermentation with ethanol 
[t_wt_eh, y_wt_eh] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_wt, jgy_vals, mgl_0), ...
                            [t_batch_eh_start, t_batch_eh_final], ...
                             cell_state, ...
                             options); 
                         
y_wt_gl_start = y_wt_eh(end,:);

% switch to glucose 
y_wt_gl_start(2) = 20*5550.7; % reset glucose concentration to 2% glucose
y_wt_gl_start(4) = 0; % reset ethanol concentration to 0
                         
% get intermediate values 
sig_wt_eh.snf1 = ones(numel(t_wt_eh, 1)); 
gr_wt_eh       = ones(numel(t_wt_eh, 1)); 

for k = 1:numel(t_wt_eh)
    [~, sig_t, ~, ~, ~, ~, ~, ~, ~, ~, ~, g_rate_t, ~] = yeast_model_update_ribosome(t_wt_eh, y_wt_eh(k,:)', par, mutant, snf1_vals_wt, jgy_vals, mgl_0);

    sig_wt_eh.snf1(k) = real(sig_t.snf1); 
    gr_wt_eh(k)       = real(g_rate_t);   
end
            
% run simulation for batch fermentation - glucose phase
% note that the final cell state for the ethanol simulation is the initial
% cell state for the glucose simulation except the glucose and ethanol
% concentrations are reset 
[t_wt_gl, y_wt_gl] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_wt, jgy_vals, mgl_0), ...
                            [t_batch_gl_start, t_batch_gl_final], ...
                             y_wt_gl_start, ...
                             options); 
                                                  
% get intermediate values 
sig_wt_gl.snf1 = ones(numel(t_wt_gl, 1)); 
gr_wt_gl       = ones(numel(t_wt_gl, 1)); 

for k = 1:numel(t_wt_gl)
    [~, sig_t, ~, ~, ~, ~, ~, ~, ~, ~, ~, g_rate_t, ~] = yeast_model_update_ribosome(t_wt_gl, y_wt_gl(k,:)', par, mutant, snf1_vals_wt, jgy_vals, mgl_0);

    sig_wt_gl.snf1(k) = real(sig_t.snf1); 
    gr_wt_gl(k)       = real(g_rate_t);   
end

disp('WT YPD --> YPE: done')                         


% batch high low snf1 fermentation experiment conditions ------------------    
disp('Constant high and low Snf1 batch fermentations - ethanol')

mutant = 'const_snf1';

snf1_vals_eh = [snf1_low snf1_high];
snf1_vals_eh = snf1_vals_eh * par.s_tot; 

% ethanol phase 
for m = 1:numel(snf1_vals_eh)
    
    % get specific snf1 structure name 
    s_val = sprintf('snf1_%g',snf1_vals_eh(m));
    s_val(s_val  == '.') = '_'; 
    
    % get reasonable initial cell state   
    cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_eh(m), jgy_vals, mgl_0);

    % run batch fermentation 
    [t, y] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_eh(m), jgy_vals, mgl_0), ...
                                [t_batch_eh_start, t_batch_eh_final], ...
                                 cell_state, ...
                                 options); 
                        
    t_batch_eh.(s_val) = t;                   
                             
    y = real(y);   
    y_batch_eh.(s_val) = y; 
    
    y_batch_gl_start.(s_val) = y_batch_eh.(s_val)(end,:);
    y_batch_gl_start.(s_val)(2) = 20*5550.7;
    y_batch_gl_start.(s_val)(4) = 0;
    
    % initialize arrays to hold intermediate values 
    met_reac_eh.(s_val).prot      = ones(numel(t_batch_eh.(s_val)), 6);
    met_reac_eh.(s_val).substrate = ones(numel(t_batch_eh.(s_val)), 6); 
    met_reac_eh.(s_val).atp       = ones(numel(t_batch_eh.(s_val)), 6); 
    met_reac_eh.(s_val).sig       = ones(numel(t_batch_eh.(s_val)), 6); 
    met_reac_eh.(s_val).flux      = ones(numel(t_batch_eh.(s_val)), 6); 
    prot_syn_eh.(s_val).alpha     = ones(numel(t_batch_eh.(s_val)), 8); 
    prot_syn_eh.(s_val).tc        = ones(numel(t_batch_eh.(s_val)), 1);  
    prot_syn_eh.(s_val).eIF_a     = ones(numel(t_batch_eh.(s_val)), 1);  
    sig_eh.(s_val).snf1           = ones(numel(t_batch_eh.(s_val)), 1);
    sig_eh.(s_val).tor            = ones(numel(t_batch_eh.(s_val)), 1); 
    gr_eh.(s_val)                 = ones(numel(t_batch_eh.(s_val)), 1); 

    % get intermediate values 
    for k = 1:length(t_batch_eh.(s_val))

        [dydt_t, sig_t, met_reac_t, ~, ~, alpha_t, ~, tRNA_t, ~, ~, ~, g_rate_t, ~] = yeast_model_update_ribosome(t_batch_eh.(s_val)(k), y_batch_eh.(s_val)(k,:)', par, mutant, snf1_vals_eh(m), jgy_vals, mgl_0);
        
        met_reac_eh.(s_val).prot(k,:)      = real(met_reac_t.prot)';
        met_reac_eh.(s_val).substrate(k,:) = real(met_reac_t.substrate)';
        met_reac_eh.(s_val).atp(k,:)       = real(met_reac_t.atp)';
        met_reac_eh.(s_val).sig(k,:)       = real(met_reac_t.sig)';
        met_reac_eh.(s_val).flux(k,:)      = real(met_reac_t.flux)';

        prot_syn_eh.(s_val).alpha(k,:) = real(table2array(struct2table(alpha_t))); 
        prot_syn_eh.(s_val).tc(k,:)    = real(tRNA_t.tc)'; 

        sig_eh.(s_val).snf1(k,:) = real(sig_t.snf1)';
        sig_eh.(s_val).tor(k,:)  = real(sig_t.tor)'; 
        
        gr_eh.(s_val)(k) = real(g_rate_t); 
    end 
end 

% glucose phase
for m = 1:numel(snf1_vals_eh)
    
    % get specific snf1 structure name 
    s_val = sprintf('snf1_%g',snf1_vals_eh(m));
    s_val(s_val  == '.') = '_'; 
    
    % run batch fermentation 
    [t, y] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_eh(m), jgy_vals, mgl_0), ...
                                [t_batch_gl_start, t_batch_gl_final], ...
                                 y_batch_gl_start.(s_val), ...
                                 options); 
                        
    t_batch_gl.(s_val) = t;                   
                             
    y = real(y);   
    y_batch_gl.(s_val) = y; 
    
    % initialize arrays to hold intermediate values 
    met_reac_gl.(s_val).prot      = ones(numel(t_batch_gl.(s_val)), 6); 
    met_reac_gl.(s_val).substrate = ones(numel(t_batch_gl.(s_val)), 6); 
    met_reac_gl.(s_val).atp       = ones(numel(t_batch_gl.(s_val)), 6); 
    met_reac_gl.(s_val).sig       = ones(numel(t_batch_gl.(s_val)), 6); 
    met_reac_gl.(s_val).flux      = ones(numel(t_batch_gl.(s_val)), 6); 
    prot_syn_gl.(s_val).alpha     = ones(numel(t_batch_gl.(s_val)), 8); 
    prot_syn_gl.(s_val).tc        = ones(numel(t_batch_gl.(s_val)), 1);  
    prot_syn_gl.(s_val).eIF_a     = ones(numel(t_batch_gl.(s_val)), 1);  
    sig_gl.(s_val).snf1           = ones(numel(t_batch_gl.(s_val)), 1);
    sig_gl.(s_val).tor            = ones(numel(t_batch_gl.(s_val)), 1); 
    gr_gl.(s_val)                 = ones(numel(t_batch_gl.(s_val)), 1); 

    % get intermediate
    for k = 1:length(t_batch_gl.(s_val))

        [dydt_t, sig_t, met_reac_t, ~, ~, alpha_t, ~, tRNA_t, ~, ~, ~, g_rate_t, ~] = yeast_model_update_ribosome(t_batch_gl.(s_val)(k), y_batch_gl.(s_val)(k,:)', par, mutant, snf1_vals_eh(m), jgy_vals, mgl_0);

        met_reac_gl.(s_val).prot(k,:)      = real(met_reac_t.prot)';
        met_reac_gl.(s_val).substrate(k,:) = real(met_reac_t.substrate)';
        met_reac_gl.(s_val).atp(k,:)       = real(met_reac_t.atp)';
        met_reac_gl.(s_val).sig(k,:)       = real(met_reac_t.sig)';
        met_reac_gl.(s_val).flux(k,:)      = real(met_reac_t.flux)';

        prot_syn_gl.(s_val).alpha(k,:) = real(table2array(struct2table(alpha_t))); 
        prot_syn_gl.(s_val).tc(k,:)    = real(tRNA_t.tc)'; 

        sig_gl.(s_val).snf1(k,:) = real(sig_t.snf1)';
        sig_gl.(s_val).tor(k,:)  = real(sig_t.tor)'; 
        
        gr_gl.(s_val)(k) = real(g_rate_t); 
    end 
end 
                         
disp('Batch fermentation: done')


% batch fermentation with high snf1 under glucose phase, low snf1 under
% ethanol phase ******************************************************
disp('Batch fermentation - low SNF1 during glucose phase, high SNF1 during ethanol phase')

disp('Starting high SNF1 ethanol phase...')

mutant = 'const_snf1';

snf1_vals_switch = [snf1_high snf1_low];
snf1_vals_switch = snf1_vals_switch * par.s_tot; 

% get specific snf1 structure name 
s_val = sprintf('snf1_%g_%g', snf1_vals_switch(1), snf1_vals_switch(2));
s_val(s_val  == '.') = '_'; 

% get reasonable initial cell state for high snf1 mutant on ethanol
cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_switch(1), jgy_vals, mgl_0);

% run batch fermentation for high snf1 mutant on ethanol 
[t, y] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_switch(1), jgy_vals, mgl_0), ...
                            [t_batch_eh_start, t_batch_eh_final], ...
                             cell_state, ...
                             options); 

t_batch_switch_eh.(s_val) = t;                   

y = real(y);   
y_batch_switch_eh.(s_val)    = y; 
y_switch_gl_start.(s_val)    = y(end,:);

% reset glucose and ethanol concentrations for glucose phase 
y_switch_gl_start.(s_val)(2) = 20*5550.7; % 2% glucose
y_switch_gl_start.(s_val)(4) = 0; % no ethanol

% initialize arrays to hold intermediate values 
met_reac_switch_eh.(s_val).prot      = ones(numel(t_batch_switch_eh.(s_val)), 6); 
met_reac_switch_eh.(s_val).substrate = ones(numel(t_batch_switch_eh.(s_val)), 6); 
met_reac_switch_eh.(s_val).atp       = ones(numel(t_batch_switch_eh.(s_val)), 6); 
met_reac_switch_eh.(s_val).sig       = ones(numel(t_batch_switch_eh.(s_val)), 6); 
met_reac_switch_eh.(s_val).flux      = ones(numel(t_batch_switch_eh.(s_val)), 6); 

prot_syn_switch_eh.(s_val).alpha     = ones(numel(t_batch_switch_eh.(s_val)), 8); 
prot_syn_switch_eh.(s_val).tc        = ones(numel(t_batch_switch_eh.(s_val)), 1);  
prot_syn_switch_eh.(s_val).eIF_a     = ones(numel(t_batch_switch_eh.(s_val)), 1);  

sig_switch_eh.(s_val).snf1 = ones(numel(t_batch_switch_eh.(s_val)), 1);
sig_switch_eh.(s_val).tor  = ones(numel(t_batch_switch_eh.(s_val)), 1); 

gr_switch_eh.(s_val) = ones(numel(t_batch_switch_eh.(s_val)), 1);

% get intermediate values 
for k = 1:length(t_batch_switch_eh.(s_val))

    [dydt_t, sig_t, met_reac_t, ~, ~, alpha_t, ~, tRNA_t, ~, ~, ~, g_rate_t, ~] = yeast_model_update_ribosome(t_batch_switch_eh.(s_val)(k), y_batch_switch_eh.(s_val)(k,:)', par, mutant, snf1_vals_switch(1), jgy_vals, mgl_0);

    met_reac_switch_eh.(s_val).prot(k,:)      = real(met_reac_t.prot)';
    met_reac_switch_eh.(s_val).substrate(k,:) = real(met_reac_t.substrate)';
    met_reac_switch_eh.(s_val).atp(k,:)       = real(met_reac_t.atp)';
    met_reac_switch_eh.(s_val).sig(k,:)       = real(met_reac_t.sig)';
    met_reac_switch_eh.(s_val).flux(k,:)      = real(met_reac_t.flux)';

    prot_syn_switch_eh.(s_val).alpha(k,:) = real(table2array(struct2table(alpha_t))); 
    prot_syn_switch_eh.(s_val).tc(k,:)    = real(tRNA_t.tc)'; 

    sig_switch_eh.(s_val).snf1(k,:) = real(sig_t.snf1)';
    sig_switch_eh.(s_val).tor(k,:)  = real(sig_t.tor)'; 
    
    gr_switch_eh.(s_val)(k) = real(g_rate_t);
end 

disp('Starting low SNF1 glucose phase...')

% batch fermentation experiment conditions

% run batch fermentation 
% note the initial cell state for the glucose phase is the same as the
% final cell state of the ethanol phase except the glucose and ethanol
% concentrations were reset
[t, y] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_switch(2), jgy_vals, mgl_0), ...
                            [t_batch_gl_start, t_batch_gl_final], ...
                             y_switch_gl_start.(s_val), ...
                             options); 

t_batch_switch_gl.(s_val) = t;                   

y = real(y);   
y_batch_switch_gl.(s_val) = y; 

% initialize arrays to hold intermediate values 
met_reac_switch_gl.(s_val).prot      = ones(numel(t_batch_switch_gl.(s_val)), 6); 
met_reac_switch_gl.(s_val).substrate = ones(numel(t_batch_switch_gl.(s_val)), 6); 
met_reac_switch_gl.(s_val).atp       = ones(numel(t_batch_switch_gl.(s_val)), 6); 
met_reac_switch_gl.(s_val).sig       = ones(numel(t_batch_switch_gl.(s_val)), 6); 
met_reac_switch_gl.(s_val).flux      = ones(numel(t_batch_switch_gl.(s_val)), 6); 

prot_syn_switch_gl.(s_val).alpha     = ones(numel(t_batch_switch_gl.(s_val)), 8); 
prot_syn_switch_gl.(s_val).tc        = ones(numel(t_batch_switch_gl.(s_val)), 1);  
prot_syn_switch_gl.(s_val).eIF_a     = ones(numel(t_batch_switch_gl.(s_val)), 1);  

sig_switch_gl.(s_val).snf1 = ones(numel(t_batch_switch_gl.(s_val)), 1);
sig_switch_gl.(s_val).tor  = ones(numel(t_batch_switch_gl.(s_val)), 1); 

gr_switch_gl.(s_val) = ones(numel(t_batch_switch_gl.(s_val)), 1);

% get intermediate values 
for k = 1:length(t_batch_switch_gl.(s_val))

    [dydt_t, sig_t, met_reac_t, ~, ~, alpha_t, ~, tRNA_t, ~, ~, ~, g_rate_t, ~] = yeast_model_update_ribosome(t_batch_switch_gl.(s_val)(k), y_batch_switch_gl.(s_val)(k,:)', par, mutant, snf1_vals_switch(2), jgy_vals, mgl_0);

    met_reac_switch_gl.(s_val).prot(k,:)      = real(met_reac_t.prot)';
    met_reac_switch_gl.(s_val).substrate(k,:) = real(met_reac_t.substrate)';
    met_reac_switch_gl.(s_val).atp(k,:)       = real(met_reac_t.atp)';
    met_reac_switch_gl.(s_val).sig(k,:)       = real(met_reac_t.sig)';
    met_reac_switch_gl.(s_val).flux(k,:)      = real(met_reac_t.flux)';

    prot_syn_switch_gl.(s_val).alpha(k,:) = real(table2array(struct2table(alpha_t))); 
    prot_syn_switch_gl.(s_val).tc(k,:)    = real(tRNA_t.tc)'; 

    sig_switch_gl.(s_val).snf1(k,:) = real(sig_t.snf1)';
    sig_switch_gl.(s_val).tor(k,:)  = real(sig_t.tor)'; 
    
    gr_switch_gl.(s_val)(k) = real(g_rate_t);
end 
                          
toc
disp('Batch fermentation: done')


%% plotting ----------------------------------------------------------------

% simulation 

% plot format
alpha1 =  0;       
alpha2 =  0;  

%y_lim_gr      = [0 0.55];
x_lim_time    = [0 28]; 
x_ticks_time  = [0 28/2 28]; 
y_lim_cell    = [0 7]*10^11;
y_ticks_cell  = [0 3 6]*10^11; 

y_lim_gl_eh   = [0 5.5]*10^5;
y_ticks_gl_eh = [0 2 4]*10^5; 

% Figure 5i-p -------------------------------------------------------------
figure;
a = 1; b = 4; 

% simulation: wt in ethanol
subplot(a, b, 1) 
hold on; 

yyaxis left 
plot([t_wt_eh; t_wt_eh(end) + t_wt_gl], [y_wt_eh(:,2); y_wt_gl(:,2)], '-', 'Color', glu_color) % glucose
plot([t_wt_eh; t_wt_eh(end) + t_wt_gl], [y_wt_eh(:,4); y_wt_gl(:,4)], '-', 'Color', eth_color) % ethanol
ylabel({'glucose,'; 'ethanol (\muM)'}, 'Color','k')
ylim(y_lim_gl_eh)
yticks(y_ticks_gl_eh)

yyaxis right 
plot([t_wt_eh; t_wt_eh(end) + t_wt_gl], [y_wt_eh(:,end); y_wt_gl(:,end)],'Color',plt_clrs.yellow) % cells
ylim(y_lim_cell)
yticks(y_ticks_cell)
ylabel('cells')

hold off;
xlim(x_lim_time)
xlabel(x_label_batch) 
xticks(x_ticks_time)

ax = gca; % set axis color as black
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.YLabel.Color = 'k'; 
ax.XLabel.Color = 'k';  

patch_background(0, t_wt_eh(end),  max(ylim), max(xlim), min(ylim), color_high_camp, color_low_camp, alpha1, alpha2)
axis square; 
box on; 


% simulation: high cAMP mutant in ethanol
snf1_val = flip(snf1_vals_eh); % flip snf1 from increasing order to decreasing order
s_val = sprintf('snf1_%g',snf1_val(2));
s_val(s_val  == '.') = '_'; 

subplot(a, b, 2) 
hold on; 
yyaxis left 

plot([t_batch_eh.(s_val); t_batch_eh.(s_val)(end) + t_batch_gl.(s_val)], [y_batch_eh.(s_val)(:,2); y_batch_gl.(s_val)(:,2)], '-', 'Color', glu_color) % glucose
plot([t_batch_eh.(s_val); t_batch_eh.(s_val)(end) + t_batch_gl.(s_val)], [y_batch_eh.(s_val)(:,4); y_batch_gl.(s_val)(:,4)], '-', 'Color', eth_color) % ethanol
ylabel({'glucose,'; 'ethanol (\muM)'})
ylim(y_lim_gl_eh)
yticks(y_ticks_gl_eh)

yyaxis right 
plot([t_batch_eh.(s_val); t_batch_eh.(s_val)(end) + t_batch_gl.(s_val)], [y_batch_eh.(s_val)(:,end); y_batch_gl.(s_val)(:,end)],'Color',plt_clrs.yellow) % cells 
ylim(y_lim_cell)
yticks(y_ticks_cell)
ylabel('cells')

hold off;
xlabel(x_label_batch) 
xlim(x_lim_time)
xticks(x_ticks_time)

ax = gca; % set axis color as black
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.YLabel.Color = 'k'; 
ax.XLabel.Color = 'k'; 

patch_background(0, t_wt_eh(end),  max(ylim), max(xlim), min(ylim), color_high_camp, color_low_camp, alpha1, alpha2)
axis square; 
box on; 



% simulation: low cAMP mutant in ethanol
subplot(a, b, 3) 
s_val = sprintf('snf1_%g',snf1_val(1));
s_val(s_val  == '.') = '_'; 
yyaxis left 
hold on; 
plot([t_batch_eh.(s_val); t_batch_eh.(s_val)(end) + t_batch_gl.(s_val)], [y_batch_eh.(s_val)(:,2); y_batch_gl.(s_val)(:,2)], '-', 'Color', glu_color) % glucose 
plot([t_batch_eh.(s_val); t_batch_eh.(s_val)(end) + t_batch_gl.(s_val)], [y_batch_eh.(s_val)(:,4); y_batch_gl.(s_val)(:,4)], '-', 'Color', eth_color) % ethanol
ylabel({'glucose,'; 'ethanol (\muM)'})
ylim(y_lim_gl_eh)
yticks(y_ticks_gl_eh)
yyaxis right 
plot([t_batch_eh.(s_val); t_batch_eh.(s_val)(end) + t_batch_gl.(s_val)], [y_batch_eh.(s_val)(:,end);  y_batch_gl.(s_val)(:,end)], '-', 'Color', plt_clrs.yellow) % cells 
ylim(y_lim_cell)
yticks(y_ticks_cell)
hold off
xlabel(x_label_batch) 
xlim(x_lim_time)
xticks(x_ticks_time)
ylabel('cells')
%xline(t_batch_eh_final)
ax = gca; % set axis color as black
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.YLabel.Color = 'k'; 
ax.XLabel.Color = 'k';  
patch_background(0, t_wt_eh(end),  max(ylim), max(xlim), min(ylim), color_high_camp, color_low_camp, alpha1, alpha2)
axis square; 
box on; 


% simulation: shift cAMP mutant in ethanol
s_val = sprintf('snf1_%g_%g', snf1_vals_switch(1), snf1_vals_switch(2));
s_val(s_val  == '.') = '_'; 

subplot(a, b, 4) 
hold on; 

yyaxis left 
plot([t_batch_switch_eh.(s_val); t_batch_switch_eh.(s_val)(end) + t_batch_switch_gl.(s_val)], [y_batch_switch_eh.(s_val)(:,2); y_batch_switch_gl.(s_val)(:,2)], '-', 'Color', glu_color) % glucose 
plot([t_batch_switch_eh.(s_val); t_batch_switch_eh.(s_val)(end) + t_batch_switch_gl.(s_val)], [y_batch_switch_eh.(s_val)(:,4); y_batch_switch_gl.(s_val)(:,4)], '-', 'Color', eth_color) % ethanol
ylabel({'glucose,'; 'ethanol (\muM)'})
ylim(y_lim_gl_eh)
yticks(y_ticks_gl_eh)
yyaxis right 

plot([t_batch_switch_eh.(s_val); t_batch_switch_eh.(s_val)(end) + t_batch_switch_gl.(s_val)], [y_batch_switch_eh.(s_val)(:,end); y_batch_switch_gl.(s_val)(:,end)], 'Color', plt_clrs.yellow) % cells 
ylim(y_lim_cell)
yticks(y_ticks_cell)
ylabel('cells')

hold off;
xlabel(x_label_batch) 
xlim(x_lim_time)
xticks(x_ticks_time)

ax = gca; % set axis color as black
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.YLabel.Color = 'k'; 
ax.XLabel.Color = 'k'; 

patch_background(0, t_wt_eh(end),  max(ylim), max(xlim), min(ylim), color_high_camp, color_low_camp, alpha1, alpha2)
axis square; 
box on; 



