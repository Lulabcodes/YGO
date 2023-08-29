clc; close all; clear

addpath('../../')
addpath(genpath('../../'))

% load simulation set, data
shared_setup;

%% constant extracellular glucose concentration

%
disp('Figure 7')
tic 

% read parameters
par = read_parametersv2; 

% initial conditions 
met.ex_amino_acids  = 0;          % 1
met.glucose         = 20*5550.7;  % 2
%met.precursor       = 1.00E+03;  % 3
met.ethanol         = 0;          % 4  
%met.in_amino_acids  = 5*10^4;    % 5  
%met.atp             = 3.00E+03;  % 6 
num.met = numel(fieldnames(met));

% proteins 
prot.r  = 0.203     * par.pro_den / par.l(1);    % 11  1 % converted amino acid fraction to protein fraction
prot.z  = 0.538     * par.pro_den / par.l(2);    % 12  2
prot.gy = 0.0736    * par.pro_den / par.l(3);    % 13  3
prot.fe = 0.0186    * par.pro_den / par.l(4);    % 14  4
prot.gn = 0.00550   * par.pro_den / par.l(5);    % 15  5
prot.mt = 0.0487    * par.pro_den / par.l(6);    % 16  6
prot.as = 0.0981    * par.pro_den / par.l(7);    % 17  7
prot.at = 0.00372   * par.pro_den / par.l(8);    % 18  10
num.prot = numel(fieldnames(prot));

% total ribosomes 
R0 = 14; % 19

cells = cells_ours_unit_conv;

% temp cell state
cell_state = [table2array(struct2table(met))';  
              table2array(struct2table(prot))'; 
              R0;                               
              cells;                            
              ];

% set ode options
Rel_tol  = 1.0E-03; 
Abs_tol  = 1.0E-06; 
options  = odeset('RelTol',Rel_tol, ...
                  'AbsTol',Abs_tol, ...
                  'NonNegative',(1:length(cell_state)));
              
              
              
% WT simulation with constant extracellular glucose concentration ---------
disp('WT simulation with constant extracellular glucose concentration: starting...') 

mutant = {'const_gl'}; % glucose concentration is held constant

% start and end times 
t_batch_start = 0; 
t_batch_final = 10;   

% get reasonable cell state 
cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, tor_vals_wt, jgy_vals, mgl_0);
disp('Calculate initial conditions: done') 

% run simulation 
disp('Running simulation...') 
[t, y]  = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, tor_vals_wt, jgy_vals, mgl_0), ...
                            [t_batch_start t_batch_final], ...
                             cell_state, ...
                             options); 

y = real(y);           
y_const_gl = y(end,:);
                         
% initialize arrays to hold intermediate variables 
met_reac_const_gl.prot      = ones(1, 6); 
met_reac_const_gl.substrate = ones(1, 6);
met_reac_const_gl.atp       = ones(1, 6); 
met_reac_const_gl.sig       = ones(1, 6);
met_reac_const_gl.flux      = ones(1, 6); 
prot_syn_const_gl.alpha     = ones(1, 8);
prot_syn_const_gl.tc        = 1;  
prot_syn_const_gl.eIF_a     = 1;  
sig_const_gl.tor            = 1;
sig_const_gl.tor            = 1; 

% get intermediate values   

[~, sig_t, met_reac_t, ~, ~, alpha_t, ~, tRNA_t, ~, ~, ~, g_rate_t, ~] = yeast_model_update_ribosome(t(end), y_const_gl', par, mutant, tor_vals_wt, jgy_vals, mgl_0);
    
met_reac_const_gl.prot(:)      = real(met_reac_t.prot)';
met_reac_const_gl.substrate(:) = real(met_reac_t.substrate)';
met_reac_const_gl.atp(:)       = real(met_reac_t.atp)';
met_reac_const_gl.sig(:)       = real(met_reac_t.sig)';
met_reac_const_gl.flux(:)      = real(met_reac_t.flux)';
    
prot_syn_const_gl.alpha(:) = real(table2array(struct2table(alpha_t))); 
prot_syn_const_gl.tc(:)    = real(tRNA_t.tc)'; 
    
sig_const_gl.tor(:) = real(sig_t.tor)';
sig_const_gl.tor(:) = real(sig_t.tor)'; 
    
g_rate_const_gl = real(g_rate_t)';


disp('WT simulation with constant glucose: done')

%% Closed- and open-loop experiments with LOW 3AT -------------------------- 
%
% initial conditions 
met.ex_amino_acids = 0;           % 1 %minimal media 
met.glucose        = 20*5550.7;   % 2  
%met.precursor      = 1.00E+03;   % 3  
met.ethanol        = 0;           % 4  
%met.in_amino_acids = 1.50E+05;   % 5  
%met.atp            = 1.00E+03;   % 6 
num.met = numel(fieldnames(met));

% initial cell state 
cell_state = [table2array(struct2table(met))';  
              table2array(struct2table(prot))'; 
              R0;                               
              cells;                            
              ];

% set ODE options 
Rel_tol  = 1.0E-03; 
Abs_tol  = 1.0E-06; 
options  = odeset('RelTol',Rel_tol, ...
                  'AbsTol',Abs_tol, ...
                  'NonNegative',(1:length(cell_state)));
              
              

% WT simulation with constant extracellular glucose concentration 
par       = read_parametersv2;
mutant    = 'const_gl';  % glucose concentration is held constant
par.c_3AT = con_3AT_low; % 3ATP concentration

% start and end times 
t_batch_start = 0; 
t_batch_final = 1; % can use small value because the initial state will already be the steady-state 
 
% get reasonable initial cell state (which will be the steady-state)
cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_wt, jgy_vals, mgl_0);
disp('Calculate initial conditions for const_gl batch simulation: done') 

% run simulation 
disp('Running const_gl batch simulation...') 
[t, y]  = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_wt, jgy_vals, mgl_0), ...
                            [t_batch_start t_batch_final], ...
                             cell_state, ...
                             options); 

c_3AT_low_wt_y          = real(y);           
c_3AT_low_wt_y_const_gl = y(end,:);

% initialize arrays to hold intermediate values
c_3AT_low_wt_met_reac_const_gl.prot      = ones(1, 6); 
c_3AT_low_wt_met_reac_const_gl.substrate = ones(1, 6);
c_3AT_low_wt_met_reac_const_gl.atp       = ones(1, 6); 
c_3AT_low_wt_met_reac_const_gl.sig       = ones(1, 6);
c_3AT_low_wt_met_reac_const_gl.flux      = ones(1, 6); 
c_3AT_low_wt_prot_syn_const_gl.alpha     = ones(1, 8);
c_3AT_low_wt_prot_syn_const_gl.tc        = 1;  
c_3AT_low_wt_prot_syn_const_gl.eIF_a     = 1;  
c_3AT_low_wt_sig_const_gl.snf1           = 1;
c_3AT_low_wt_sig_const_gl.tor            = 1; 


% get intermediate values   
[~, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, rib_t, tRNA_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t, g_rate_t, ribo_rate_t] = yeast_model_update_ribosome(t(end), c_3AT_low_wt_y_const_gl', par, mutant, snf1_vals_wt, jgy_vals, mgl_0);
    
c_3AT_low_wt_met_reac_const_gl.prot(:)      = real(met_reac_t.prot)';
c_3AT_low_wt_met_reac_const_gl.substrate(:) = real(met_reac_t.substrate)';
c_3AT_low_wt_met_reac_const_gl.atp(:)       = real(met_reac_t.atp)';
c_3AT_low_wt_met_reac_const_gl.sig(:)       = real(met_reac_t.sig)';
c_3AT_low_wt_met_reac_const_gl.flux(:)      = real(met_reac_t.flux)';
    
c_3AT_low_wt_prot_syn_const_gl.alpha(:) = real(table2array(struct2table(alpha_t))); 
c_3AT_low_wt_prot_syn_const_gl.tc(:)    = real(tRNA_t.tc)'; 
    
c_3AT_low_wt_sig_const_gl.snf1(:) = real(sig_t.snf1)';
c_3AT_low_wt_sig_const_gl.tor(:)  = real(sig_t.tor)'; 
    
c_3AT_low_wt_g_rate_const_gl      = real(g_rate_t)';
 
disp('WT simulation with constant glucose: done')



% TORC1 mutant simulation with constant extracellular glucose
% concentration -----------------------------------------------------------
mutant = 'const_tor_gl'; % s = const and glucose concentration is held constant 

% start and end times 
t_batch_start = 0; 
t_batch_final = 1; % can choose small value because the initial state will already be the steady state 

% run simulation for different TORC1 mutants 
for m = 1:numel(tor_vals)
    disp(m)
    
    % get specific snf1 structure name 
    t_val = sprintf('tval_%g',tor_vals(m));
    t_val(t_val  == '.') = '_'; 
  
    % get reasonable initial cell state (which will be the steady-state for
    % the given conditons)
    cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, tor_vals(m), jgy_vals, mgl_0);
    disp('Calculate initial conditions for const_gl batch simulation: done') 

    % run simulation 
    disp('Running const_gl batch simulation...') 
    [t, y]  = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, tor_vals(m), jgy_vals, mgl_0), ...
                                [t_batch_start t_batch_final], ...
                                 cell_state, ...
                                 options); 
                                        
    y = real(y);  
    c_3AT_low_y_snf1.(t_val) = y(end,:);
    
    % initialize arrays to hold intermediate values 
    c_3AT_low_met_reac_snf1.(t_val).prot      = ones(1, 6); 
    c_3AT_low_met_reac_snf1.(t_val).substrate = ones(1, 6); 
    c_3AT_low_met_reac_snf1.(t_val).atp       = ones(1, 6);
    c_3AT_low_met_reac_snf1.(t_val).sig       = ones(1, 6); 
    c_3AT_low_met_reac_snf1.(t_val).flux      = ones(1, 6); 
    c_3AT_low_prot_syn_snf1.(t_val).alpha     = ones(1, 8); 
    c_3AT_low_prot_syn_snf1.(t_val).tc        = 1;  
    c_3AT_low_sig_snf1.(t_val).snf1           = 1;
    c_3AT_low_sig_snf1.(t_val).tor            = 1; 
    c_3AT_low_g_rate_snf1.(t_val)             = 1;

    % get intermediate values 
    [~, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, rib_t, tRNA_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t, g_rate_t, ribo_rate_t] = yeast_model_update_ribosome(t(end), c_3AT_low_y_snf1.(t_val)', par, mutant, tor_vals(m), jgy_vals, mgl_0);
    c_3AT_low_met_reac_snf1.(t_val).prot(:)      = real(met_reac_t.prot)';
    c_3AT_low_met_reac_snf1.(t_val).substrate(:) = real(met_reac_t.substrate)';
    c_3AT_low_met_reac_snf1.(t_val).atp(:)       = real(met_reac_t.atp)';
    c_3AT_low_met_reac_snf1.(t_val).sig(:)       = real(met_reac_t.sig)';
    c_3AT_low_met_reac_snf1.(t_val).flux(:)      = real(met_reac_t.flux)';

    c_3AT_low_prot_syn_snf1.(t_val).alpha(:) = real(table2array(struct2table(alpha_t))); 
    c_3AT_low_prot_syn_snf1.(t_val).tc(:)    = real(tRNA_t.tc)'; 

    c_3AT_low_sig_snf1.(t_val).snf1(:) = real(sig_t.snf1)';
    c_3AT_low_sig_snf1.(t_val).tor(:)  = real(sig_t.tor)'; 
        
    c_3AT_low_g_rate_snf1.(t_val) = g_rate_t; 
end 
  

% Closed- and open-loop experiments with HIGH 3AT -------------------------

% initial conditions 
met.ex_amino_acids = 0;            % 1 %minimal media = 0
met.glucose        = 20*5550.7;    % 2  
%met.precursor      = 1.00E+03;    % 3  
met.ethanol        = 0;            % 4  
%met.in_amino_acids = 1.50E+05;    % 5  
%met.atp            = 1.00E+03;    % 7 
num.met = numel(fieldnames(met));

cell_state = [table2array(struct2table(met))';  % 1 - 7
              table2array(struct2table(prot))'; % 8 - 17 
              R0;                               % 18
              cells;                            % 19
              ];
          
% set ODE options           
Rel_tol  = 1.0E-03; 
Abs_tol  = 1.0E-06; 
options  = odeset('RelTol',Rel_tol, ...
                  'AbsTol',Abs_tol, ...
                  'NonNegative',(1:length(cell_state)));

% WT simulation with constant extracellular glucose concentration 
par = read_parametersv2;
mutant    = 'const_gl'; % glucose concentration is held constant
par.c_3AT = con_3AT_high;  % 3ATP concentration

% start and end times 
t_batch_start = 0; 
t_batch_final = 1; % can choose small value because the initial cell state will already be the steady state
 
% get reasonable initial cell state
cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, snf1_vals_wt, jgy_vals, mgl_0);
disp('Calculate initial conditions for const_gl batch simulation: done') 

% run simulation
disp('Running const_gl batch simulation...') 
[t, y]  = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, snf1_vals_wt, jgy_vals, mgl_0), ...
                            [t_batch_start t_batch_final], ...
                             cell_state, ...
                             options); 

c_3AT_high20_wt_y          = real(y);           
c_3AT_high20_wt_y_const_gl = y(end,:);
                        
% initialize arrays to hold intermediate values 
c_3AT_high20_wt_met_reac_const_gl.prot      = ones(1, 6); 
c_3AT_high20_wt_met_reac_const_gl.substrate = ones(1, 6);
c_3AT_high20_wt_met_reac_const_gl.atp       = ones(1, 6); 
c_3AT_high20_wt_met_reac_const_gl.sig       = ones(1, 6);
c_3AT_high20_wt_met_reac_const_gl.flux      = ones(1, 6); 
c_3AT_high20_wt_prot_syn_const_gl.alpha     = ones(1, 8);
c_3AT_high20_wt_prot_syn_const_gl.tc        = 1;  
c_3AT_high20_wt_prot_syn_const_gl.eIF_a     = 1;  
c_3AT_high20_wt_sig_const_gl.snf1           = 1;
c_3AT_high20_wt_sig_const_gl.tor            = 1; 


% get intermediate values   
[~, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, rib_t, tRNA_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t, g_rate_t, ribo_rate_t] = yeast_model_update_ribosome(t(end), c_3AT_high20_wt_y_const_gl', par, mutant, snf1_vals_wt, jgy_vals, mgl_0);
    
c_3AT_high20_wt_met_reac_const_gl.prot(:)      = real(met_reac_t.prot)';
c_3AT_high20_wt_met_reac_const_gl.substrate(:) = real(met_reac_t.substrate)';
c_3AT_high20_wt_met_reac_const_gl.atp(:)       = real(met_reac_t.atp)';
c_3AT_high20_wt_met_reac_const_gl.sig(:)       = real(met_reac_t.sig)';
c_3AT_high20_wt_met_reac_const_gl.flux(:)      = real(met_reac_t.flux)';
    
c_3AT_high20_wt_prot_syn_const_gl.alpha(:) = real(table2array(struct2table(alpha_t))); 
c_3AT_high20_wt_prot_syn_const_gl.tc(:)    = real(tRNA_t.tc)'; 
    
c_3AT_high20_wt_sig_const_gl.snf1(:) = real(sig_t.snf1)';
c_3AT_high20_wt_sig_const_gl.tor(:)  = real(sig_t.tor)'; 
    
c_3AT_high20_wt_g_rate_const_gl = real(g_rate_t)';
 
disp('WT simulation with constant glucose: done')


% TORC1 mutant simulation with constant extracellular glucose
% concentration -----------------------------------------------------------
mutant = 'const_tor_gl'; % s = const and glucose concentration is held constant 

% start and end times 
t_batch_start = 0; 
t_batch_final = 1; % can choose small value since the initial state will already be the steady state

% run simulation for different TORC1 mutants 
for m = 1:numel(tor_vals)
    disp(m)
    
    % get specific snf1 structure name 
    t_val = sprintf('tval_%g',tor_vals(m));
    t_val(t_val  == '.') = '_'; 
  
    % get reasonable cell state (will be steady state for the given conditions )
    cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, tor_vals(m), jgy_vals, mgl_0);
    disp('Calculate initial conditions for const_gl batch simulation: done') 

    % run simulation 
    disp('Running const_gl batch simulation...') 
    [t, y]  = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, tor_vals(m), jgy_vals, mgl_0), ...
                                [t_batch_start t_batch_final], ...
                                 cell_state, ...
                                 options); 
                                        
    y = real(y);  
    c_3AT_high20_y_snf1.(t_val) = y(end,:);
    
    % initialize arrays to hold intermediate values 
    c_3AT_high20_met_reac_snf1.(t_val).prot      = ones(1, 6); 
    c_3AT_high20_met_reac_snf1.(t_val).substrate = ones(1, 6); 
    c_3AT_high20_met_reac_snf1.(t_val).atp       = ones(1, 6);
    c_3AT_high20_met_reac_snf1.(t_val).sig       = ones(1, 6); 
    c_3AT_high20_met_reac_snf1.(t_val).flux      = ones(1, 6); 

    c_3AT_high20_prot_syn_snf1.(t_val).alpha     = ones(1, 8); 
    c_3AT_high20_prot_syn_snf1.(t_val).tc        = 1;  
    
    c_3AT_high20_sig_snf1.(t_val).snf1 = 1;
    c_3AT_high20_sig_snf1.(t_val).tor  = 1; 
    
    c_3AT_high20_g_rate_snf1.(t_val) = 1;

    % get intermediate values 
    [~, sig_t, met_reac_t, prot_syn_rate_t, beta_t, alpha_t, rib_t, tRNA_t, eIF_a_s_t, eIF_a_tau_t, other_met_reac_t, g_rate_t, ribo_rate_t] = yeast_model_update_ribosome(t(end), c_3AT_high20_y_snf1.(t_val)', par, mutant, tor_vals(m), jgy_vals, mgl_0);
    c_3AT_high20_met_reac_snf1.(t_val).prot(:)      = real(met_reac_t.prot)';
    c_3AT_high20_met_reac_snf1.(t_val).substrate(:) = real(met_reac_t.substrate)';
    c_3AT_high20_met_reac_snf1.(t_val).atp(:)       = real(met_reac_t.atp)';
    c_3AT_high20_met_reac_snf1.(t_val).sig(:)       = real(met_reac_t.sig)';
    c_3AT_high20_met_reac_snf1.(t_val).flux(:)      = real(met_reac_t.flux)';

    c_3AT_high20_prot_syn_snf1.(t_val).alpha(:) = real(table2array(struct2table(alpha_t))); 
    c_3AT_high20_prot_syn_snf1.(t_val).tc(:)    = real(tRNA_t.tc)'; 

    c_3AT_high20_sig_snf1.(t_val).snf1(:) = real(sig_t.snf1)';
    c_3AT_high20_sig_snf1.(t_val).tor(:)  = real(sig_t.tor)'; 
        
    c_3AT_high20_g_rate_snf1.(t_val) = g_rate_t; 
end 
  
disp('Constant Snf1 and with high mgl constant glucose simulation: done')  
%}

%% plotting ---------------------------------------------------------------

% plot format
x_gcn2_lim       = [0.1 1]; 
x_gcn2_ticks     = [0.1 0.3 0.9]; 
y_lim_gnc2_3at   = [0.00 0.42];
y_ticks_gcn2_3at = [0 0.2 0.4]; 


% SIM: 3AT = 20 mM
gr = table2array(struct2table(c_3AT_high20_g_rate_snf1));

% select points to plot with dots 
int_index = [1     5     9    13    17   23    27    31   35   39  41];
% select points to plot with dots 
x_gcn2 = tor_to_gcn2(tor_vals); % convert torc1 values to gcn2 values 
sim_dots_x = x_gcn2(int_index);
sim_dots_y = gr(int_index);
[max_gr_val, max_gr_index] = max(sim_dots_y); 

figure; 

hold on;
plot(tor_to_gcn2(tor_vals),table2array(struct2table(c_3AT_high20_g_rate_snf1)), '-', 'Color', eth_color) % plot entire simulation
plot(x_gcn2(int_index),gr(int_index),     'o', 'MarkerFaceColor', eth_color) % plot selected points with dots 
plot(sim_dots_x(max_gr_index),max_gr_val, 'o', 'MarkerFaceColor', wt_camp_eh) % use a darker color for the largest growth rate
yline(c_3AT_high20_wt_g_rate_const_gl, '--', 'Color', wt_camp_eh) % wt 
yline(g_rate_const_gl, '--') % wt growth rate 
hold off;

ylabel('growth rate (h^{-1})')
xlabel('Gcn2 equiv.')
ylim(y_lim_gnc2_3at)
yticks(y_ticks_gcn2_3at)
xlim(x_gcn2_lim)  
xticks(x_gcn2_ticks)
set(gca, 'Xscale','log')
set(gca,'XMinorTick','off')
box on;
axis square; 
