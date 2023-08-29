clc; close all; clear; 

%% plot format 
addpath('../')
addpath(genpath('./'))
set(0,'DefaultTextColor', 'k')
set(0, 'DefaultLineLineWidth',  1.2);
set(0, 'DefaultLineMarkerSize', 7.5);
set(0, 'DefaultLineMarkerEdgeColor', 'k');
set(0, 'DefaultAxesFontSize', 8)
set(0,'DefaultTextColor', 'k')
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

plt_labels = plot_labels;
paper_figure_color; 
plot_num; 
x_label_batch = 'time (h)';

%% import data
% data_gr_opt; 
% data_camp_gl_growth_log_growth_gl_camp_timecourse;
% data_gcn2_gl_growth_log_growth_gl_timecourse;
% data_gcn2_eh_growth_log_growth_timecourse;
% data_camp_methylgloxal;
% data_camp_other_sugars; 
% data_camp_YPD_E;
% data_camp_YPE_D;
% data_camp_YPD_diauxic_shift; 

%load_simulation_setting; 
data_unit_conv; 
par = read_parametersv2;

%% camp and gcn2 conversion to snf1 and torc1
snf1_to_camp = @(x) 1-x/par.s_tot; % snf1 to camp
tor_to_gcn2  = @(x) 1-x./par.tau_tot; 

%% mutants 
snf1_vals_wt = 0;               
snf1_vals    = linspace(0, 1, 8); %[0 0.001 0.139 0.4 0.5 0.7 0.8 1];
snf1_vals    = snf1_vals * par.s_tot; 

snf1_vals_25    = linspace(0, 1, 50); %[0 0.001 0.139 0.4 0.5 0.7 0.8 1];
snf1_vals_25    = snf1_vals_25 * par.s_tot; 

snf1_low    = 0.15;
snf1_high   = 1;

jgy_vals    = 0;

mgl_0       = 0;

mgl_selec   = linspace(0,20,20) * 10^3;
mgl_selec_5 = linspace(0,1,5) * 10^3;
mgl_low   = 0.1*10^3;
mgl_high  = 1*10^3;

tor_vals_wt  = 0; 

tor_vals = 0.1:0.02:1;
tor_vals = tor_vals * par.tau_tot; 

% 3AT
con_3AT_selec = linspace(0,20,20); %mM
con_3AT_low   = 6; 
con_3AT_high  = 20; 

label_x_shift = -0.3;
label_y_shift = 1.2;
fontsize      = 13;

tor_init = par.tau_tot; 

% initial conditions 
met.ex_amino_acids = 2*10^5;            % 1
met.glucose        = 20*gl_gPerL_to_uM; % 2
met.precursor      = 1.00E+03;          % 3
met.ethanol        = 0;                 % 4  
met.in_amino_acids = 5*10^4;            % 5  
met.atp            = 3.00E+03;          % 7 
num.met            = numel(fieldnames(met));

% proteins 
prot.r   = 0.203     * par.pro_den / par.l(1);    % 11  1 % converted amino acid fraction to protein fraction
prot.z   = 0.538     * par.pro_den / par.l(2);    % 12  2
prot.gy  = 0.0736    * par.pro_den / par.l(3);    % 13  3
prot.fe  = 0.0186    * par.pro_den / par.l(4);    % 14  4
prot.gn  = 0.00550   * par.pro_den / par.l(5);    % 15  5
prot.mt  = 0.0487    * par.pro_den / par.l(6);    % 16  6
prot.as  = 0.0981    * par.pro_den / par.l(7);    % 17  7
prot.at  = 0.00372   * par.pro_den / par.l(8);    % 18  8

num.prot = numel(fieldnames(prot));

% total ribosomes 
R0 = 14; % 19

% cells 
cells_ours_unit_conv   = 0.1* OD_to_gdw .* gdw_to_cell;
cell_mason_unit_conv   = 2.600e+08;

cells = cell_mason_unit_conv; %cells_ours_unit_conv; % 24
t_batch_log = 13.5;


% temp cell state
cell_state = [
              table2array(struct2table(met))';  
              table2array(struct2table(prot))'; 
              R0;                               
              cells;                           
              ];
         
% set ode options
Rel_tol  = 1.0E-06; 
Abs_tol  = 1.0E-02* ones(1,length(cell_state));
options  = odeset('RelTol',Rel_tol, ...
                  'AbsTol',Abs_tol, ...
                  'NonNegative',(1:length(cell_state)));


