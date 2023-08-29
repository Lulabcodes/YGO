function [met_reac] = module_met(met, prot, snf1, par, mutant_type, jht, mgl)

%% metabolites and proteins

aa_ex = met(1); 
gl    = met(2); 
pc    = met(3); 
eh    = met(4); 
aa_in = met(5); 
ae    = met(6); 

%p_r  = prot(1); 
%p_z  = prot(2);
p_gy = prot(3);
p_fe = prot(4);
p_gn = prot(5);
p_mt = prot(6);
p_as = prot(7);
p_at = prot(8);

%% metabolic reactions 
% enzymes
met_reac.prot = [p_gy; % gy    1
                 p_fe; % fe    2
                 p_gn; % gn    3
                 p_mt; % mt    4
                 p_as; % as    5
                 p_at; % at    6
                 ];
             
met_reac.rate = [par.k_gy; % gy    1
                 par.k_fe; % fe    2
                 par.k_gn; % gn    3
                 par.k_mt; % mt    4
                 par.k_as; % as    5
                 par.k_at; % at    6
                 ];             

% substrate contribution - hill(S, K, n)
met_reac.substrate = [hill(gl,    par.K_gy + gl/par.K_gy_gl, 1)* hill(par.K_mgl, mgl, 1);  % gy    1
                      hill(pc,    par.K_fe,                  par.n_fe_p);                  % fe    2
                      hill(eh,    par.K_gn_eh,               1);                           % gn    3
                      hill(pc,    par.K_mt,                  par.n_mt_p);                  % mt    4
                      hill(pc,    par.K_as_p,                par.n_as_p);                  % as    5
                      hill(aa_ex, par.K_at_aa,               1);                           % at    6
                     ];
                  
% ATP contribution                  
met_reac.atp = [mmki(ae, par.I_gy_ae, par.n_gy_ae);                % gy    1
                1;                                                 % fe    2
                1;                                                 % gn    3
                mmki(ae, par.I_mt_ae, par.n_mt_ae);                % mt    4
                hill(ae, par.K_as_ae, 1);                          % as    5
                hill(ae, par.K_at_ae, par.n_at_ae);                % at    6
                ];   

% signaling contribution            
met_reac.sig = [mmki(snf1,  par.I_gy_s,  par.n_gy_s)  * mmki(pc, par.I_gy_pc, par.n_gy_pc);  % gy    1
                1;                                                           % fe    2
                hill(snf1,  par.D_gn_s,  par.n_gn_s);                        % gn    3
                mmki(pc,    par.I_mt_pc, par.n_mt_pc);                       % mt    4
                mmki(aa_in, par.I_as_ac, par.n_as_ac);                       % as    5
                mmki(aa_in, par.I_at_aa, par.n_at_aa);                       % at    6
                ];      

met_reac.snf1 = [mmki(snf1,  par.I_gy_s,  par.n_gy_s);  % gy    1
                1;                                      % fe    2
                hill(snf1,  par.D_gn_s,  par.n_gn_s);   % gn    3
                1;                                      % mt    4
                1;                                      % as    5
                1;                                      % at    6
                ];   

met_reac.pc   = [mmki(pc, par.I_gy_pc, par.n_gy_pc);         % gy    1
                1;                                           % fe    2
                1;                                           % gn    3
                mmki(pc,    par.I_mt_pc, par.n_mt_pc);       % mt    4
                mmki(aa_in, par.I_as_ac, par.n_as_ac);       % as    5
                mmki(aa_in, par.I_at_aa, par.n_at_aa);       % at    6
                ];  

% fluxes            
met_reac.flux = met_reac.prot .* met_reac.rate .* met_reac.substrate .* met_reac.atp .* met_reac.sig;
met_reac.flux(5)  =  hill(par.K_3AT, par.c_3AT, 1)*met_reac.flux(5); 


