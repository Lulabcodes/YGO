clc; close all; clear

addpath('../../')
addpath(genpath('../../'))

% load simulation set, data
shared_setup;

%% Simulation and plotting for Figure 6 and related SI
%
tic
disp('Figure 6: YPD diauxic shift')

par = read_parametersv2;
cells  = cell_mason_unit_conv; 

%  batch fermentation with wt, high and low constant snf1 - GLUCOSE
%  ------------------------------------------------------------------------
disp('WT batch fermentation - glucose')

mutant = 'none';

met.ex_amino_acids = 2*10^5;             % 1 % rich media
met.glucose        = 20*gl_gPerL_to_uM;  % 2  
met.ethanol        = 0;                  % 4  
%met.in_amino_acids = 1.50E+05;          % 5  
%met.atp            = 1.00E+03;          % 6
num.met = numel(fieldnames(met));

% initial cell state
cell_state = [table2array(struct2table(met))';  
              table2array(struct2table(prot))'; 
              R0;                               
              cells;                            
              ];

% batch wt fermentation experiment conditions          
t_batch_gl_start = 0; 
t_batch_gl_final = 24;

% get reasonable initial cell state 
cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_wt, jgy_vals, mgl_0);
disp('Calculated initial conditions: done')

% run batch fermentation 
[t_wt_gl, y_wt_gl] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_wt, jgy_vals, mgl_0), ...
                            [t_batch_gl_start, t_batch_gl_final], ...
                             cell_state, ...
                             options); 
                         
% get intermediate values 
sig_wt_gl.snf1 = ones(numel(t_wt_gl, 1)); 
gr_wt_gl       = ones(numel(t_wt_gl, 1)); 

for k = 1:numel(t_wt_gl)
    [~, sig_t, ~, ~, ~, ~, ~, ~, ~, ~, ~, g_rate_t, ~] = yeast_model_update_ribosome(t_wt_gl, y_wt_gl(k,:)', par, mutant, snf1_vals_wt, jgy_vals, mgl_0);

    sig_wt_gl.snf1(k) = real(sig_t.snf1); 
    gr_wt_gl(k)       = real(g_rate_t);   
end
                         
disp('Batch simulation: done')                         


% batch high low snf1 fermentation experiment conditions ------------------    
disp('High and low Snf1 batch fermentations - glucose')

mutant = 'const_snf1';
snf1_vals_gl = [snf1_low snf1_high];
snf1_vals_gl = snf1_vals_gl * par.s_tot; 


for m = 1:numel(snf1_vals_gl)
    
    % get specific snf1 structure name 
    s_val = sprintf('snf1_%g',snf1_vals_gl(m));
    s_val(s_val  == '.') = '_'; 
    
    % get reasonable cell state   
    cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_gl(m), jgy_vals, mgl_0);

    % run batch fermentation 
    [t, y] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_gl(m), jgy_vals, mgl_0), ...
                                [t_batch_gl_start, t_batch_gl_final], ...
                                 cell_state, ...
                                 options); 
                        
    t_batch_gl.(s_val) = t;                   
                             
    y = real(y);   
    y_batch_gl.(s_val) = y; 
    
    % get intermediate values 
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
                         
disp('Batch fermentation: done')
toc 


% batch fermentation with wt, high and low constant snf1 - ethanol --------
% batch fermentation with low snf1 under glucose phase, high snf1 under ethanol phase
%--------------------------------------------------------------------------

% wt 
disp('Batch fermentation - low SNF1 during glucose phase, high SNF1 during ethanol phase')

tic
disp('WT batch fermentation...')

mutant = 'none';

% start and end times 
t_batch_switch_start = 0; 
t_batch_switch_final = t_batch_gl_final; 

% initial conditions 
met.ex_amino_acids = 2.0*10^5;                   % 1 % rich media
met.glucose        = 20*gl_gPerL_to_uM;          % 2  
%met.precursor      = 1.00E+03;                  % 3  
met.ethanol        = 0; %21 * 21706.0994139353;  % 4  
%met.in_amino_acids = 1.50E+05;                  % 5  
%met.atp            = 1.00E+03;                  % 6 
num.met = numel(fieldnames(met));

% initial cell state 
cell_state = [table2array(struct2table(met))';  
              table2array(struct2table(prot))';  
              R0;                               
              cells;                            
              ];
          

% get reasonable initial cell state
cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_wt, jgy_vals, mgl_0);
disp('Calculated initial conditions for WT batch fermentation: done')

% run batch fermentation 
[t_wt_switch, y_wt_switch] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_wt, jgy_vals, mgl_0), ...
                            [t_batch_switch_start, t_batch_switch_final], ...
                             cell_state, ...
                             options);  
                        

                         
% low snf1 mutant 
disp('Starting low SNF1 glucose phase...')

snf1_vals_switch = [snf1_low snf1_high];
snf1_vals_switch = snf1_vals_switch * par.s_tot; 

mutant = 'const_snf1';

% get specific snf1 structure name 
s_val = sprintf('snf1_%g_%g', snf1_vals_switch(1), snf1_vals_switch(2));
s_val(s_val  == '.') = '_'; 


% start and end times 
t_batch_switch_start = 0; 
t_batch_switch_final = t_batch_gl_final;   

% get reasonable cell state   
cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_switch, jgy_vals, mgl_0);

% run batch fermentation 
[t, y] = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_switch, jgy_vals, mgl_0), ...
                            [t_batch_switch_start, t_batch_switch_final], ...
                             cell_state, ...
                             options); 

t_batch_switch.(s_val) = t;                   

y = real(y);   
y_batch_switch.(s_val) = y; 

% get intermediate values 
met_reac_switch.(s_val).prot      = ones(numel(t_batch_switch.(s_val)), 6); 
met_reac_switch.(s_val).substrate = ones(numel(t_batch_switch.(s_val)), 6); 
met_reac_switch.(s_val).atp       = ones(numel(t_batch_switch.(s_val)), 6); 
met_reac_switch.(s_val).sig       = ones(numel(t_batch_switch.(s_val)), 6); 
met_reac_switch.(s_val).flux      = ones(numel(t_batch_switch.(s_val)), 6); 

prot_syn_switch.(s_val).alpha     = ones(numel(t_batch_switch.(s_val)), 8); 
prot_syn_switch.(s_val).tc        = ones(numel(t_batch_switch.(s_val)), 1);  
prot_syn_switch.(s_val).eIF_a     = ones(numel(t_batch_switch.(s_val)), 1);  

sig_switch.(s_val).snf1 = ones(numel(t_batch_switch.(s_val)), 1);
sig_switch.(s_val).tor  = ones(numel(t_batch_switch.(s_val)), 1); 

gr_switch.(s_val) = ones(numel(t_batch_switch.(s_val)), 1);

for k = 1:length(t_batch_switch.(s_val))

    [dydt_t, sig_t, met_reac_t, ~, ~, alpha_t, ~, tRNA_t, ~, ~, ~, g_rate_t, ~] = yeast_model_update_ribosome(t_batch_switch.(s_val)(k), y_batch_switch.(s_val)(k,:)', par, mutant, snf1_vals_switch, jgy_vals, mgl_0);

    met_reac_switch.(s_val).prot(k,:)      = real(met_reac_t.prot)';
    met_reac_switch.(s_val).substrate(k,:) = real(met_reac_t.substrate)';
    met_reac_switch.(s_val).atp(k,:)       = real(met_reac_t.atp)';
    met_reac_switch.(s_val).sig(k,:)       = real(met_reac_t.sig)';
    met_reac_switch.(s_val).flux(k,:)      = real(met_reac_t.flux)';

    prot_syn_switch.(s_val).alpha(k,:) = real(table2array(struct2table(alpha_t))); 
    prot_syn_switch.(s_val).tc(k,:)    = real(tRNA_t.tc)'; 

    sig_switch.(s_val).snf1(k,:) = real(sig_t.snf1)';
    sig_switch.(s_val).tor(k,:)  = real(sig_t.tor)'; 
    
    gr_switch.(s_val)(k) = real(g_rate_t);
end 
                          
toc
disp('Batch fermentation: done')





% plotting ----------------------------------------------------------------

% plot format 
alpha1 = 0.0;     
alpha2 = 0.0;

x_lim_time   = [0 24]; 
x_ticks_time = [0 12 24];
y_lim_gr     = [0 0.59]; 
y_lim_cell   = [0 10*10^11]; 
y_lim_glu_eh = [0 2*10^5];

% Figure 6:  diauxic shift ------------------------------------------------
figure;
a = 1; b = 4; 

% simulation
                            
% simulation: wt in gl
subplot(a, b, 1) 
hold on; 

yyaxis left 
plot(t_wt_gl, y_wt_gl(:,2), '-', 'Color', glu_color) 
plot(t_wt_gl, y_wt_gl(:,4), '-', 'Color', eth_color) 
ylabel({'glucose,'; 'ethanol (\muM)'})
ylim(y_lim_glu_eh)

yyaxis right 
plot(t_wt_gl, y_wt_gl(:,end), 'Color', plt_clrs.yellow) 
ylim(y_lim_cell)
ylabel('cells')

hold off;
xlabel(x_label_batch) 
xlim(x_lim_time)
xticks(x_ticks_time)

ax = gca; % set axis color as black
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

axis square; 
box on; 


% simulation:  'high cAMP mutant' in gl
snf1_val = flip(snf1_vals_gl); % flip snf1 from increasing order to decreasing order
s_val = sprintf('snf1_%g',snf1_val(2));
s_val(s_val  == '.') = '_'; 

subplot(a, b, 2) 
hold on; 

yyaxis left 
plot(t_batch_gl.(s_val), y_batch_gl.(s_val)(:,2), '-', 'Color', glu_color) 
plot(t_batch_gl.(s_val), y_batch_gl.(s_val)(:,4), '-', 'Color', eth_color) 
ylabel({'glucose,'; 'ethanol (\muM)'})
ylim(y_lim_glu_eh)

yyaxis right 
plot(t_batch_gl.(s_val), y_batch_gl.(s_val)(:,end), 'Color', plt_clrs.yellow) 
ylim(y_lim_cell)
ylabel('cells')
ylim(y_lim_cell)

hold off;
xlabel(x_label_batch) 
xlim(x_lim_time)
xticks(x_ticks_time)

ax = gca; % set axis color as black
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

axis square; 
box on; 


% simulation:  low cAMP mutant in gl
s_val = sprintf('snf1_%g',snf1_val(1));
s_val(s_val  == '.') = '_'; 

subplot(a, b, 3) 
hold on; 

yyaxis left 
plot(t_batch_gl.(s_val), y_batch_gl.(s_val)(:,2), '-', 'Color', glu_color) 
plot(t_batch_gl.(s_val), y_batch_gl.(s_val)(:,4), '-', 'Color', eth_color) 
ylabel({'glucose,'; 'ethanol (\muM)'})
ylim(y_lim_glu_eh)

yyaxis right 
plot(t_batch_gl.(s_val), y_batch_gl.(s_val)(:,end), 'Color', plt_clrs.yellow) 
ylim(y_lim_cell)
ylabel('cells')
hold off;

xlabel(x_label_batch) 
xticks(x_ticks_time)
xlim(x_lim_time)

ax = gca; % set axis color as black
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

axis square; 
box on; 


% simulation: shift cAMP mutant
s_val = sprintf('snf1_%g_%g', snf1_vals_switch(1), snf1_vals_switch(2));
s_val(s_val  == '.') = '_'; 

subplot(a, b, 4) 
yyaxis left 
hold on; 
plot(t_batch_switch.(s_val), y_batch_switch.(s_val)(:,2), '-', 'Color', glu_color) 
plot(t_batch_switch.(s_val), y_batch_switch.(s_val)(:,4), '-', 'Color', eth_color) 
ylabel({'glucose,'; 'ethanol (\muM)'})
ylim(y_lim_glu_eh)

yyaxis right 
plot(t_batch_switch.(s_val), y_batch_switch.(s_val)(:,end),'Color',plt_clrs.yellow) 
ylim(y_lim_cell)
ylabel('cells')

hold off;
xlim(x_lim_time)
xticks(x_ticks_time)
xlabel(x_label_batch) 
xline(t_batch_gl_final)

ax = gca; % set axis color as black
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

axis square; 
box on; 





%}
