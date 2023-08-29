clc; close all; clear

addpath('../../')
addpath(genpath('../../'))

% load simulation set, data
shared_setup;


%% Simulation and plotting for Figure 7 and related SI 
% constant extracellular glucose concentration

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
              


% tor mutant simulation with constant extracellular glucose concentration -
disp('TOR mutant simulation with constant extracellular glucose concentration: starting...') 

mutant = 'const_tor_gl'; % s = const and glucose concentration is held constant 

% start and end times 
t_batch_start = 0; 
t_batch_final = 1; % can choose small time because the initial state will already be the steady-state

%  run simulation for different torc1 mutants 
for m = 1:numel(tor_vals)
    disp(m)
    
    % get specific tor structure name 
    t_val = sprintf('tor_%g',tor_vals(m));
    t_val(t_val  == '.') = '_'; 
    
    % get reasonable cell state (will be the steady-state for the given conditions)
    cell_state = get_cell_state_update_ribosome(cell_state, par, mutant, tor_vals(m), jgy_vals, mgl_0);
    disp('Calculate initial conditions: done') 
    
    % run simulation 
    disp('Running simulation...') 
    [t, y]  = ode15s(@(t,z) yeast_model_update_ribosome(t, z, par, mutant, tor_vals(m), jgy_vals, mgl_0), ...
                                [t_batch_start t_batch_final], ...
                                 cell_state, ...
                                 options); 
                             
    y = real(y);  
    y_tor_gl.(t_val) = y(end,:);
    
    % initialize arrays to hold intermediate values 
    met_reac_tor_gl.(t_val).prot      = ones(1, 6); 
    met_reac_tor_gl.(t_val).substrate = ones(1, 6); 
    met_reac_tor_gl.(t_val).atp       = ones(1, 6);
    met_reac_tor_gl.(t_val).sig       = ones(1, 6); 
    met_reac_tor_gl.(t_val).flux      = ones(1, 6); 
    prot_syn_tor_gl.(t_val).alpha     = ones(1, 8); 
    prot_syn_tor_gl.(t_val).rib.others= ones(1, length(plt_labels.others)); 
    prot_syn_tor_gl.(t_val).tc        = 1;  
    sig_tor_gl.(t_val).tor = 1;
    sig_tor_gl.(t_val).tor = 1; 
    g_rate_tor_gl.(t_val) = 1;

    % get intermediate values 
    [~, sig_t, met_reac_t, ~, ~, alpha_t, rib_t, tRNA_t, ~, ~, ~, g_rate_t, ~] = yeast_model_update_ribosome(t(end), y_tor_gl.(t_val)', par, mutant, tor_vals(m), jgy_vals, mgl_0);
    met_reac_tor_gl.(t_val).prot(:)      = real(met_reac_t.prot)';
    met_reac_tor_gl.(t_val).substrate(:) = real(met_reac_t.substrate)';
    met_reac_tor_gl.(t_val).atp(:)       = real(met_reac_t.atp)';
    met_reac_tor_gl.(t_val).sig(:)       = real(met_reac_t.sig)';
    met_reac_tor_gl.(t_val).flux(:)      = real(met_reac_t.flux)';

    prot_syn_tor_gl.(t_val).alpha(:)     = real(table2array(struct2table(alpha_t))); 
    prot_syn_tor_gl.(t_val).tc(:)        = real(tRNA_t.tc)'; 
    prot_syn_tor_gl.(t_val).rib.others   = real(rib_t.others)'; 
    sig_tor_gl.(t_val).tor(:)            = real(sig_t.tor)';
    sig_tor_gl.(t_val).tor(:)            = real(sig_t.tor)'; 
        
    g_rate_tor_gl.(t_val) = g_rate_t; 
end 

disp('Constant TOR and constant glucose simulation: done')
%}


%% plotting ---------------------------------------------------------------

% Fig7 --------------------------------------------------------------------
% plot format
x_gcn2_lim       = [0.1 1]; 
x_gcn2_ticks     = [0.1 0.3 0.9]; 

y_lim_gnc2_gl    = [0 0.42];
y_ticks_gcn2_gl  = [0 0.2 0.4];


% Sim: gr vs Gcn2 in glucose
x_gcn2 = tor_to_gcn2(tor_vals); % convert torc1 values to gcn2 values 
gr     = table2array(struct2table(g_rate_tor_gl));

% select points to plot with dots 
int_index = [1     5     9    13    17   23    27    31   35   39  41];
sim_dots_x = x_gcn2(int_index);
sim_dots_y = gr(int_index);
[max_gr_val, max_gr_index] = max(sim_dots_y);

figure;
hold on;
plot(x_gcn2, table2array(struct2table(g_rate_tor_gl)), '-', 'Color', glu_gnc2_color) % plot entire simulation 
plot(x_gcn2(int_index), gr(int_index),   'o', 'MarkerFaceColor', glu_gnc2_color) % plot selected points with dots 
plot(sim_dots_x(max_gr_index),max_gr_val,'o','MarkerFaceColor',color_wt) % use a darker color for the largest growth rate 
yline( g_rate_const_gl, '--')
hold off;

ylabel('growth rate (h^{-1})')
xlabel('Gcn2 equiv.')
ylim(y_lim_gnc2_gl)
yticks(y_ticks_gcn2_gl)
xlim(x_gcn2_lim)
xlim([min(tor_to_gcn2(tor_vals)) max(tor_to_gcn2(tor_vals))])
xlim(x_gcn2_lim)
xticks(x_gcn2_ticks)
set(gca, 'Xscale','log')
set(gca,'XMinorTick','off')
box on; 
axis square; 

