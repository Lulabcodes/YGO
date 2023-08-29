clc; close all; clear; 

addpath('../../')
addpath(genpath('../../'))

% load simulation set, data
shared_setup;

%% Simulation and plotting for Figure 5: YPD --> YPE and related SI
%{\
disp('Figure 5: YPD --> YPE')

% not change parameter normalized time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
YPD_time = 12+2 ; %12
YPE_time = 36-YPD_time; %22.5;
cells  = cell_mason_unit_conv; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  batch fermentation with wt, high and low constant snf1 - GLUCOSE
%  ------------------------------------------------------------------------
disp('WT batch fermentation - glucose')
tic
met.ex_amino_acids = 2.0*10^5;          % 1  % rich media 
met.glucose        = 20*gl_gPerL_to_uM; % 2  
%met.precursor      = 1.00E+03;         % 3  
met.ethanol        = 0;                 % 4  
%met.in_amino_acids = 1.50E+05;         % 5  
%met.atp            = 1.00E+03;         % 6 
num.met = numel(fieldnames(met));

cell_state = [table2array(struct2table(met))';  
              table2array(struct2table(prot))'; 
              R0;                               
              cells;                            
              ];

% batch wt fermentation experiment conditions          
t_batch_gl_start = 0; 
t_batch_gl_final = YPD_time; 


mutant = 'none';
cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_wt, jgy_vals, mgl_0);
disp('Calculated initial conditions: done')

% batch fermentation - glucose phase
[t_wt_gl, y_wt_gl] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_wt, jgy_vals, mgl_0), ...
                            [t_batch_gl_start, t_batch_gl_final], ...
                             cell_state, ...
                             options); 
                         
y_wt_eh_start = y_wt_gl(end,:);
y_wt_eh_start(2) = 0;
y_wt_eh_start(4) = 20 * 21706.0994139353;
                         
% get intermediate values 
sig_wt_gl.snf1 = ones(numel(t_wt_gl, 1)); 
gr_wt_gl       = ones(numel(t_wt_gl, 1)); 

for k = 1:numel(t_wt_gl)
    [~, sig_t, ~, ~, ~, ~, ~, ~, ~, ~, ~, g_rate_t, ~] = yeast_model_update_ribosome(t_wt_gl, y_wt_gl(k,:)', par, mutant, snf1_vals_wt, jgy_vals, mgl_0);

    sig_wt_gl.snf1(k) = real(sig_t.snf1); 
    gr_wt_gl(k)       = real(g_rate_t);   
end


t_batch_eh_final = YPE_time;              
% batch fermentation - ethanol phase
[t_wt_eh, y_wt_eh] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_wt, jgy_vals, mgl_0), ...
                            [t_batch_gl_start, t_batch_eh_final], ...
                             y_wt_eh_start, ...
                             options); 
                                                  
% get intermediate values 
sig_wt_eh.snf1 = ones(numel(t_wt_eh, 1)); 
gr_wt_eh       = ones(numel(t_wt_eh, 1)); 

for k = 1:numel(t_wt_eh)
    [~, sig_t, ~, ~, ~, ~, ~, ~, ~, ~, ~, g_rate_t, ~] = yeast_model_update_ribosome(t_wt_eh, y_wt_eh(k,:)', par, mutant, snf1_vals_wt, jgy_vals, mgl_0);

    sig_wt_eh.snf1(k) = real(sig_t.snf1); 
    gr_wt_eh(k)       = real(g_rate_t);   
end

disp('WT YPD --> YPE: done')                         
toc 

% batch high low snf1 fermentation experiment conditions ------------------    
disp('High and low Snf1 batch fermentations - glucose')
tic

mutant = 'const_snf1';

snf1_vals_gl = [snf1_low snf1_high];
snf1_vals_gl = snf1_vals_gl * par.s_tot; 

% read parameters
%par = read_parametersv2;

% glucose phase 
for m = 1:numel(snf1_vals_gl)
    
    % get specific snf1 structure name 
    s_val = sprintf('snf1_%g',snf1_vals_gl(m));
    s_val(s_val  == '.') = '_'; 
    
    % get cell state   
    cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_gl(m), jgy_vals, mgl_0);

    % batch fermentation 
    [t, y] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_gl(m), jgy_vals, mgl_0), ...
                                [t_batch_gl_start, t_batch_gl_final], ...
                                 cell_state, ...
                                 options); 
                        
    t_batch_gl.(s_val) = t;                   
                             
    y = real(y);   
    y_batch_gl.(s_val) = y; 
    
    y_batch_eh_start.(s_val) = y_batch_gl.(s_val)(end,:);
    y_batch_eh_start.(s_val)(2) = 0;
    y_batch_eh_start.(s_val)(4) = 20 * 21706.0994139353;
    
    % get intermediate values 
    met_reac_gl.(s_val).prot      = ones(numel(t_batch_gl.(s_val)), 6); 
    met_reac_gl.(s_val).substrate = ones(numel(t_batch_gl.(s_val)), 6); 
    met_reac_gl.(s_val).atp       = ones(numel(t_batch_gl.(s_val)), 6); 
    met_reac_gl.(s_val).sig       = ones(numel(t_batch_gl.(s_val)), 6); 
    met_reac_gl.(s_val).flux      = ones(numel(t_batch_gl.(s_val)), 6); 

    prot_syn_gl.(s_val).alpha     =  ones(numel(t_batch_gl.(s_val)), 8); 
    prot_syn_gl.(s_val).tc        = ones(numel(t_batch_gl.(s_val)), 1);  
    prot_syn_gl.(s_val).eIF_a     = ones(numel(t_batch_gl.(s_val)), 1);  

    sig_gl.(s_val).snf1           = ones(numel(t_batch_gl.(s_val)), 1);
    sig_gl.(s_val).tor            = ones(numel(t_batch_gl.(s_val)), 1); 
    
    gr_gl.(s_val)                 = ones(numel(t_batch_gl.(s_val)), 1); 

    for k = 1:length(t_batch_gl.(s_val))

        [dydt_t, sig_t, met_reac_t, ~, ~, alpha_t, ~, tRNA_t, ~, ~, ~, g_rate_t, ~] = yeast_model_update_ribosome(t_batch_gl.(s_val)(k), y_batch_gl.(s_val)(k,:)', par, mutant, snf1_vals_gl(m), jgy_vals, mgl_0);

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

t_batch_eh_final = YPE_time;

% ethanol phase
for m = 1:numel(snf1_vals_gl)
    
    % get specific snf1 structure name 
    s_val = sprintf('snf1_%g',snf1_vals_gl(m));
    s_val(s_val  == '.') = '_'; 
    
    % batch fermentation 
    [t, y] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_gl(m), jgy_vals, mgl_0), ...
                                [t_batch_gl_start, t_batch_eh_final], ...
                                 y_batch_eh_start.(s_val), ...
                                 options); 
                        
    t_batch_eh.(s_val) = t;                   
                             
    y = real(y);   
    y_batch_eh.(s_val) = y; 
    
    % get intermediate values 
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

    for k = 1:length(t_batch_eh.(s_val))

        [dydt_t, sig_t, met_reac_t, ~, ~, alpha_t, ~, tRNA_t, ~, ~, ~, g_rate_t, ~] = yeast_model_update_ribosome(t_batch_eh.(s_val)(k), y_batch_eh.(s_val)(k,:)', par, mutant, snf1_vals_gl(m), jgy_vals, mgl_0);

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
                         
disp('Batch fermentation: done')
toc 




% batch fermentation with low snf1 under glucose phase, high snf1 under ethanol phase

disp('Starting low SNF1 glucose phase...')

snf1_vals_switch = [snf1_low snf1_high];
snf1_vals_switch = snf1_vals_switch * par.s_tot; 

% batch fermentation experiment conditions

t_batch_switch_start = 0; 
t_batch_switch_final = YPD_time;   

% get specific snf1 structure name 
s_val = sprintf('snf1_%g_%g', snf1_vals_switch(1), snf1_vals_switch(2));
s_val(s_val  == '.') = '_'; 

mutant = 'const_snf1';

% get cell state   
cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_switch(1), jgy_vals, mgl_0);

% batch fermentation 
[t, y] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_switch(1), jgy_vals, mgl_0), ...
                            [t_batch_switch_start, t_batch_switch_final], ...
                             cell_state, ...
                             options); 

t_batch_switch_gl.(s_val) = t;                   

y = real(y);   
y_batch_switch_gl.(s_val)   = y; 
y_switch_eh_start.(s_val) = y(end,:);
y_switch_eh_start.(s_val)(2) = 0;
y_switch_eh_start.(s_val)(4) = 20 * 21706.0994139353;

% get intermediate values 
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

for k = 1:length(t_batch_switch_gl.(s_val))

    [dydt_t, sig_t, met_reac_t, ~, ~, alpha_t, ~, tRNA_t, ~, ~, ~, g_rate_t, ~] = yeast_model_update_ribosome(t_batch_switch_gl.(s_val)(k), y_batch_switch_gl.(s_val)(k,:)', par, mutant, snf1_vals_switch(1), jgy_vals, mgl_0);

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

disp('Starting high SNF1 ethanol phase...')

% batch fermentation start and end times

t_batch_switch_start = 0; 
t_batch_switch_final = YPE_time;   

% batch fermentation 
[t, y] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_switch(2), jgy_vals, mgl_0), ...
                            [t_batch_switch_start, t_batch_switch_final], ...
                             y_switch_eh_start.(s_val), ...
                             options); 

t_batch_switch_eh.(s_val) = t;                   

y = real(y);   
y_batch_switch_eh.(s_val) = y; 

% get intermediate values 
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

for k = 1:length(t_batch_switch_eh.(s_val))

    [dydt_t, sig_t, met_reac_t, ~, ~, alpha_t, ~, tRNA_t, ~, ~, ~, g_rate_t, ~] = yeast_model_update_ribosome(t_batch_switch_eh.(s_val)(k), y_batch_switch_eh.(s_val)(k,:)', par, mutant, snf1_vals_switch(2), jgy_vals, mgl_0);

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
                          
toc
disp('Batch fermentation: done')

%% plotting ----------------------------------------------------------------

% plot format
alpha1    =   0.9;              
alpha2    =   0.9; 

x_lim_sim      = [0 28];
x_ticks_sim    = [0 28/2 28]; 
y_lim_cell     = [0 9.5*10^11]; 
y_ticks_cell   = [0 4 8]*10^11; 

y_lim_glu_eh   = [-100 50*10^4];
y_ticks_glu_eh = [0 2 4]*10^5;


% Fig 5a-h  ---------------------------------------------------------------
figure;
a = 1; b = 4; 

% simulation --------------------------------------------------------------
                            
% simulation: wt in gl
subplot(a, b, 1) 
hold on; 

yyaxis left 
plot([t_wt_gl; t_wt_gl(end) + t_wt_eh], [y_wt_gl(:,2); y_wt_eh(:,2)], '-', 'Color', glu_color) % glucose
plot([t_wt_gl; t_wt_gl(end) + t_wt_eh], [y_wt_gl(:,4); y_wt_eh(:,4)], '-', 'Color', eth_color) % ethanol 
ylabel({'glucose,'; 'ethanol (\muM)'})
ylim(y_lim_glu_eh)
yticks(y_ticks_glu_eh)

yyaxis right 
plot([t_wt_gl; t_wt_gl(end) + t_wt_eh], [y_wt_gl(:,end); y_wt_eh(:,end)], 'Color', plt_clrs.yellow) % cells 
ylim(y_lim_cell)
yticks(y_ticks_cell)

hold off;
xlabel(x_label_batch) 
ylabel('cells')

xlim(x_lim_sim); 
xticks(x_ticks_sim);  

ax = gca; % set axis color as black
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.YLabel.Color   = 'k';
ax.XLabel.Color   = 'k';

%patch_background(0, t_wt_gl(end),  max(ylim), max(xlim), min(ylim), color_low_camp, color_high_camp, alpha1, alpha2)
axis square; 
box on;  


%  simulation: high cAMP mutant in gl
snf1_val = flip(snf1_vals_gl); % flip snf1 from increasing order to decreasing order
s_val = sprintf('snf1_%g',snf1_val(2)); % select high camp mutant
s_val(s_val  == '.') = '_'; 

subplot(a, b, 2) 
hold on; 
yyaxis left 

plot([t_batch_gl.(s_val); t_batch_gl.(s_val)(end) + t_batch_eh.(s_val)], [y_batch_gl.(s_val)(:,2); y_batch_eh.(s_val)(:,2)], '-', 'Color', glu_color) % glucose
plot([t_batch_gl.(s_val); t_batch_gl.(s_val)(end) + t_batch_eh.(s_val)], [y_batch_gl.(s_val)(:,4); y_batch_eh.(s_val)(:,4)], '-', 'Color', eth_color) % ethanol 
ylabel({'glucose,'; 'ethanol (\muM)'})
ylim(y_lim_glu_eh)
yticks(y_ticks_glu_eh)

yyaxis right 
plot([t_batch_gl.(s_val); t_batch_gl.(s_val)(end) + t_batch_eh.(s_val)], [y_batch_gl.(s_val)(:,end); y_batch_eh.(s_val)(:,end)], 'Color', plt_clrs.yellow) % cells
ylim(y_lim_cell)
yticks(y_ticks_cell)

hold off;
xlabel(x_label_batch) 
ylabel('cells')

xlim(x_lim_sim); 
xticks(x_ticks_sim);

ax = gca; % set axis color as black
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.YLabel.Color = 'k'; 
ax.XLabel.Color = 'k'; 

%patch_background(0, t_wt_gl(end),  max(ylim), max(xlim), min(ylim), color_low_camp, color_high_camp, alpha1, alpha2)
axis square; 
box on; 


%  simulation: low cAMP mutant in gl
s_val = sprintf('snf1_%g',snf1_val(1));
s_val(s_val  == '.') = '_'; % select low camp mutant 

subplot(a, b, 3) 
hold on; 

yyaxis left 
plot([t_batch_gl.(s_val); t_batch_gl.(s_val)(end) + t_batch_eh.(s_val)], [y_batch_gl.(s_val)(:,2); y_batch_eh.(s_val)(:,2)], '-', 'Color', glu_color) % glucose
plot([t_batch_gl.(s_val); t_batch_gl.(s_val)(end) + t_batch_eh.(s_val)], [y_batch_gl.(s_val)(:,4); y_batch_eh.(s_val)(:,4)], '-', 'Color', eth_color) % ethanol 
ylabel({'glucose,'; 'ethanol (\muM)'})
ylim(y_lim_glu_eh)
yticks(y_ticks_glu_eh)

yyaxis right 
plot([t_batch_gl.(s_val); t_batch_gl.(s_val)(end) + t_batch_eh.(s_val)], [y_batch_gl.(s_val)(:,end);  y_batch_eh.(s_val)(:,end)], '-', 'Color', plt_clrs.yellow) % cells
ylim(y_lim_cell)
yticks(y_ticks_cell)

hold off;
xlabel(x_label_batch) 
ylabel('cells')

xlim(x_lim_sim); 
xticks(x_ticks_sim);

ax = gca; % set axis color as black
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.YLabel.Color = 'k'; 
ax.XLabel.Color = 'k'; 

%patch_background(0, t_wt_gl(end),  max(ylim), max(xlim), min(ylim), color_low_camp, color_high_camp, alpha1, alpha2)
axis square; 
box on; 


% simulation: switch cAMP 
s_val = sprintf('snf1_%g_%g', snf1_vals_switch(1), snf1_vals_switch(2)); % select switch camp mutant 
s_val(s_val  == '.') = '_'; 

subplot(a, b, 4) 
hold on; 

yyaxis left 
plot([t_batch_switch_gl.(s_val); t_batch_switch_gl.(s_val)(end) + t_batch_switch_eh.(s_val)], [y_batch_switch_gl.(s_val)(:,2); y_batch_switch_eh.(s_val)(:,2)],'-', 'Color', glu_color) % glucose
plot([t_batch_switch_gl.(s_val); t_batch_switch_gl.(s_val)(end) + t_batch_switch_eh.(s_val)], [y_batch_switch_gl.(s_val)(:,4); y_batch_switch_eh.(s_val)(:,4)],'-', 'Color', eth_color) % ethanol
ylabel({'glucose,'; 'ethanol (\muM)'})
ylim(y_lim_glu_eh)
yticks(y_ticks_glu_eh)

yyaxis right 
plot([t_batch_switch_gl.(s_val); t_batch_switch_gl.(s_val)(end) + t_batch_switch_eh.(s_val)], [y_batch_switch_gl.(s_val)(:,end); y_batch_switch_eh.(s_val)(:,end)],'Color',plt_clrs.yellow) % cells 
ylim(y_lim_cell)
yticks(y_ticks_cell)
ylabel('cells')

hold off;
xlabel(x_label_batch) 
xlim(x_lim_sim);
xticks(x_ticks_sim);

ax = gca; % set axis color as black
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.YLabel.Color = 'k'; 
ax.XLabel.Color = 'k'; 

%patch_background(0, t_wt_gl(end),  max(ylim), max(xlim), min(ylim), color_low_camp, color_high_camp, alpha1, alpha2)
axis square; 
box on; 








%}
