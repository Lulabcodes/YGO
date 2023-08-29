clc; close all; clear;

addpath('../../')
addpath(genpath('../../'))

% load simulation set, data
shared_setup;

%% Figure 7 related SI
% constant extracellular glucose concentration

%

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

% Supplementary Figure 13 -------------------------------------------------
% REZ during growth on glucose 

% convert y_tor_gl to y_tor_gl_mat 
% y_tor_gl: a struct of model variable length(snf1_vals_25) number of fields
% y_tor_gl_mat: a matrix: [num_y by length(snf1_vals_25)]
y_tor_gl_size = [length(fieldnames(num_y)), length(tor_vals)];
y_tor_gl_mat  = reshape(table2array(struct2table(y_tor_gl)), y_tor_gl_size)';

% calculate R, E, Z fraction:
total_protein_con = sum(y_tor_gl_mat(1,num_y.r:num_y.gn).* par.l(num_prot.r:num_prot.gn)',2); % total protein concentration 

R_tor_gl_25 = y_tor_gl_mat(:,num_y.r).*par.l(num_prot.r)'./total_protein_con(1); % R
Z_tor_gl_25 = y_tor_gl_mat(:,num_y.z).*par.l(num_prot.z)'./total_protein_con(1); % E
E_tor_gl_25   = sum(y_tor_gl_mat(:,num_y.gy:num_y.at).*par.l(num_prot.gy:num_prot.at)',2)./total_protein_con(1); % Z
Eas_tor_gl_25 = sum(y_tor_gl_mat(:,num_y.as).*par.l(num_prot.as)',2)./total_protein_con(1); % Eas

% plot format 
x_lim_gnc2   = [0.1 0.9];
x_gcn2_ticks = [0.1 0.3 0.9];
y_ticks      = [10^-5 10^-3 10^-1]; % y_ticks for E fractions.


% plot REZ vs torc1 (gcn2 equiv.) 
figure;
a = 1; b = 2; 

% R, E, Z vs GCN2
subplot(a, b, 1)
hold on 
plot(tor_to_gcn2(tor_vals), R_tor_gl_25, '-', 'Color', plt_clrs.yellow)
plot(tor_to_gcn2(tor_vals), E_tor_gl_25, '-', 'Color', plt_clrs.green)
plot(tor_to_gcn2(tor_vals), Z_tor_gl_25, '-', 'Color', plt_clrs.gray)
xline(tor_to_gcn2(sig_const_gl.tor), '--')
hold off;

ylabel({'protein fraction'})
ylim([0 1.05*max(ylim)])
xlabel('Gcn2 equiv.')
xlim(x_lim_gnc2)
xticks(x_gcn2_ticks)
legend('R ','E ','Z','WT GCN2')
legend box off;
set(gca,'Xscale','log')
set(gca,'XMinorTick','off')
axis square; 
box on; 

% Eas vs GCN2
subplot(a, b, 2)
hold on;
plot(tor_to_gcn2(tor_vals),Eas_tor_gl_25,'-','Color',plt_clrs.green)
xline(tor_to_gcn2(sig_const_gl.tor), '--')
hold off;
ylabel({'E_{as} fraction'})
xlabel('Gcn2 equiv.')
xlim(x_lim_gnc2)
xticks(x_gcn2_ticks)
ylim([0.08 0.35])
set(gca,'Xscale','log')
set(gca,'XMinorTick','off')
axis square; 
box on; 
 
toc


%}