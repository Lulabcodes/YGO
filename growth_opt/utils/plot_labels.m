function [label] = plot_labels

label.met = {{'amino'; 'acids_{ex}'};
             'glucose_{ex}';
             'precursor';
             'ethanol';
             {'amino'; 'acids_{in}'};
             'atp'};
         
label.prot = {'R';
              'Z';
              'gy';
              'fe';
              'gn';
              'mt';
              'as';
              'at'};   
          
label.rib_cells = {'R_{0}';
                  'number of cells'};
              
label.rib = {'R_{af}';
             'R_{i}';
             'R_{at_j}';
             'R_{at}';
             'R_{as_j}';
             'R_{as}';
             'R_{f}'};              
              
label.met_reac.flux = {'J_{gy}';
                       'J_{fe}';
                       'J_{gn}';
                       'J_{mt}';
                       'J_{as}';
                       'J_{at}'};   

label.met_reac.prot = {'E_{gy}';
                       'E_{fe}';
                       'E_{gn}';
                       'E_{mt}';
                       'E_{as}';
                       'E_{at}'};                    
                   
label.met_reac.substrate = {{'k*[G_{l}]/' , '(K+[G_{l}])'};
                            {'k*[P_{c}]/' , '(K+[P_{c}])'};
                            {'k*[E_{h}]/' , '(K+[E_{h}])'};
                            {'k*[P_{c}]/' , '(K+[P_{c}])'};
                            {'k*[P_{c}]/' , '(K+[P_{c}])'};
                            {'k*[A_{a}^{e}]/' , '(K+[A_{a}^{e}])'}};  

label.met_reac.atp = {{'(I^{n})/' , '(I^{n} + [A_{e}]^{n})'};
                      'none';
                      {'([A_{e}]^{n})/' , '(I^{n} + [A_{e}]^{n})'};
                      {'(I^{n})/' , '(I^{n} + [A_{e}]^{n})'};
                      {'([A_{e}]^{n})/' , '(I^{n} + [A_{e}]^{n})'};
                      {'([A_{e}]^{n})/' , '(I^{n} + [A_{e}]^{n})'};
                      {'([A_{e}]^{n})/' , '(I^{n} + [A_{e}]^{n})'}};                            
                      
label.met_reac.sig = {{'(I^{n})/' , '(I_{gy,s}^{n_{gy,s}} + [s^{*}]^{n_{gy,s}})'};
                      'none';
                      {'([s^{*}]^{n})/' , '(I^{n} + [s^{*}]^{n})'};
                      {'([s^{*}]^{n})/' , '(I^{n} + [s^{*}]^{n})'};
                      {'(I^{n})/' , '(I^{n} + [A_{a}^{i}]^{n})'};
                      {'(I^{n})/' , '(I^{n} + [A_{a}^{i}]^{n})'}}; 


label.met_reac.sig = {{'(I_{gy,s}^{n_{gy,s}})/'    , '(I_{gy,s}^{n_{gy,s}} + [s^{*}]^{n_{gy,s}})'};
                      'none';
                      {'([s^{*}]^{n_{gn,s}})/'     , '(I_{gn,s}^{n_{gn,s}} + [s^{*}]^{n_{gn,s}})'};
                      {'( [s^{*}]^{n_{mt,s}})/'     , '(I_{mt,s}^{n_{mt,s}} + [s^{*}]^{n_{mt,s}})'};
                      {'(I_{as,ai}^{n_{as,ai}})/' , '(I_{as,ai}^{n_{as,ai}} + [A_{a}^{i}]^{n_{as,s}})'};
                      {'(I_{at,ai}^{n_{at,ai}})/' , '(I_{at,ai}^{n_{at,ai}} + [A_{a}^{i}]^{n_{at,s}})'}}; 


label.met_reac.sig = {{'[Pc],[s^{*}] term'};
                      'none';
                      {'[s^{*}] term'};
                      {'[Pc],[s^{*}] term'};
                      {'[Pc],[Aa] term'};
                      {'(I_{at,ai}^{n_{at,ai}})/' , '(I_{at,ai}^{n_{at,ai}} + [A_{a}^{i}]^{n_{at,s}})'}}; 


label.met_reac.snf1 = {{'[s^{*}] term'};
                      'none';
                      {'[s^{*}] term'};
                      {'none'};
                      {'none'};
                      {'(I_{at,ai}^{n_{at,ai}})/' , '(I_{at,ai}^{n_{at,ai}} + [A_{a}^{i}]^{n_{at,s}})'}}; 

label.met_reac.pc = {{'[Pc] term'};
                      'none';
                      {'none'};
                      {'[Pc] term'};
                      {'[Aa] term'};
                      {'(I_{at,ai}^{n_{at,ai}})/' , '(I_{at,ai}^{n_{at,ai}} + [A_{a}^{i}]^{n_{at,s}})'}}; 

%}
                      
label.prot_syn.alpha = {'\alpha_{r}';
                        '\alpha_{z}';
                        '\alpha_{gy}';
                        '\alpha_{fe}';
                        '\alpha_{gn}';
                        '\alpha_{mt}';
                        '\alpha_{as}';
                        '\alpha_{at}'};   
                    
label.prot_syn.beta  = {'\beta_{r}';
                        '\beta_{z}';
                        '\beta_{gy}';
                        '\beta_{fe}';
                        '\beta_{gn}';
                        '\beta_{mt}';
                        '\beta_{as}';
                        '\beta_{at}'};                      

label.prot_syn.prot_syn_rate = {'l_{j} \beta_{j}'};   
label.prot_syn.eIF_a_s       = {'K^{n}/','(K^{n} + [s^{*}]^{n})'}; 
label.prot_syn.eIF_a_tau     = {'[\tau^{*}]^{n})/','(K^{n} + [\tau^{*}]^{n})'};
                    
label.prot_syn.tc = {'t_{c}'};

label.rat_j = {'R_{at,r}';
               'R_{at,z}';
               'R_{at,gy}';
               'R_{at,fe}';
               'R_{at,gn}';
               'R_{at,mt}';
               'R_{at,as}';
               'R_{at,at}'};

label.sig.snf1 = {'SNF1 activity'};
label.sig.tor  = {'TORC1 activity'};

label.others = {'trna_tc';
                'trna_tu';
                'trna_contri_tc';
                'trna_contri_tu';
                'atp_contri';
                'sum \alpha';
                'eIF_a_s'; 
                'eIF_a_tau';
                'Rf';
                'f_to_at';
                'Rat';
                'f_to_as';                    
                'Ras';
                'prot_sys_rate'; 
                'g_rate';
                'fat';
                't_0';
                'aa'; 
                'aa term'};
               
label.atp_flux = {'q_{gy} J_{gy}';
                  'q_{mt} J_{mt}';
                  'q_{gn} J_{gn}';
                  'q_{as} J_{as}';
                  'q_{p} l_{j}\beta_{j}'};  
              
label.atp_flux_chem = {'q_{gy} J_{gy}';
                       'q_{mt} J_{mt}';
                       'q_{gn} J_{gn}';
                       'q_{as} J_{as}';
                       'q_{p} l_{j}\beta_{j}';
                       'J_{eo}'};                

end 