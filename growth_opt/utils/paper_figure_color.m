%% plotting set up 
plt_clrs       = plot_colors;

% patch two region color
color1 = [17  17  17]/255; % lighter gray      
color2 = [1   1   1]/255;  % darker gray 

%% ----------------------------------------------------------------------------------------------------------------
%                               growth opt colors 
% ----------------------------------------------------------------------------------------------------------------- 

color_camp1  = [158,   1,  66] ./ 255;
color_camp2  = [213,  62,  79] ./ 255;
color_camp3  = [244, 109,  67] ./ 255;
color_camp4  = [253, 174,  97] ./ 255;
color_camp5  = [254, 224, 139] ./ 255;
color_camp6  = [230, 245, 152] ./ 255;
color_camp7  = [171, 221, 164] ./ 255;
color_camp8  = [102, 194, 165] ./ 255;
color_camp9  = [50,  136, 189] ./ 255;
color_camp10 = [30,   82, 113] ./ 255;
color_camp11 = [127,115,180]   ./ 255; 
color_camp12 = [94,   79, 162] ./ 255; 

camp_colors = [color_camp1;
               color_camp2;
               color_camp3; 
               color_camp4;
               color_camp5;
               color_camp6;
               color_camp7;
               color_camp8;
               color_camp9;
               color_camp10;];

color_wt   = 'k';

% background for figure 5,6 -----------------------------------------------
color_low_camp  = '#F8F5F1'; 
color_high_camp = '#C7EAFB'; 
color_high_wt   = '#EAE0D6';

% cAMP colors -------------------------------------------------------------

% open-loop max color: 
camp_gr_max   = plt_clrs.blue;               
camp_cell_max = plt_clrs.red;

% GCN2 colors -------------------------------------------------------
glu_gnc2_color = plt_clrs.lightgray; 
wt_gnc2_color  = 'k';                

eth_gnc2_color    = plt_clrs.lightgray; 
wt_eth_gnc2_color = 'k';                

glu_color_max = 'k';         
eth_color_max = plt_clrs.blue; 
fru_color_max = plt_clrs.gray;   
suc_color_max = plt_clrs.orange; 

%% Gr opt Figure 1 color 

plt_clrs.orange      = '#fbb03a'; 

% orange --> white 
orange_light_gradients = {'#fbb03a';
                         '#ffbd5f';
                         '#ffcb82';
                         '#ffd9a3';
                         '#ffe6c3';
                         '#fff2e1';
                         '#ffffff';};
% orange --> black 
orange_dark_gradients = {'#fbb03a';
                         '#cd9033';
                         '#a1722b';
                         '#775423';
                         '#4f391b';
                         '#2a1f12';
                         '#000000';};

% green --> white 
green_light_gradients = {'#4cbd94';
                         '#70c8a5';
                         '#8fd4b7';
                         '#acdfc8';
                         '#c8eada';
                         '#e3f4ed';
                         '#ffffff';};
% green --> black 
green_dark_gradients = {'#4cbd94';
                        '#419b7a';
                        '#377a60';
                        '#2c5a48';
                        '#203d31';
                        '#15211c';
                        '#000000';};

% blue --> white 
blue_light_gradients = {'#2078b4';
                        '#568dc1';
                        '#7ca3cd';
                        '#9eb9da';
                        '#bfd0e6';
                        '#dfe7f3';
                        '#ffffff';};

% blue --> black 
blue_dark_gradients = {'#2078b4'
                       '#226393'
                       '#214f74'
                       '#1e3c56'
                       '#18293a'
                       '#111820'
                       '#000000'};

% red --> white 
red_light_gradients = {'#e31f26';
                       '#ef5547';
                       '#f97a69';
                       '#ff9d8d';
                       '#F9B7B7';
                       '#ffe0d9';
                       '#ffffff';};

% red --> black 
red_dark_gradients = {'#e31f26';
                      '#ba2121';
                      '#93201c';
                      '#6e1d18';
                      '#4a1813';
                      '#29110b';
                      '#000000';};

% gray --> black 
gray_light_gradients = {'#b8b8b8';
                        '#c1c1c1';
                        '#c9c9c9';
                        '#d2d2d2';
                        '#dbdbdb';
                        '#e4e4e4';
                        '#ededed';
                        '#f6f6f6';
                        '#ffffff';};


gray_dark_gradients = {'#b8b8b8';
                       '#9f9f9f';
                       '#868686';
                       '#6f6f6f';
                       '#585858';
                       '#424242';
                       '#2d2d2d';
                       '#1a1a1a';
                       '#000000';};

plt_clrs.mgl0  = plt_clrs.lightgray; 
plt_clrs.mgl6  = plt_clrs.green; 
plt_clrs.mgl8  = plt_clrs.green;
plt_clrs.mgl10 = plt_clrs.green;

plt_clrs.lightorange = orange_light_gradients{3}; 
plt_clrs.darkorange  = orange_dark_gradients{2}; 
plt_clrs.lightgreen  = green_light_gradients{3};  
plt_clrs.darkgreen   = green_dark_gradients{3}; 
plt_clrs.lightblue   = blue_light_gradients{4}; 
plt_clrs.darkblue    = blue_dark_gradients{3}; 
plt_clrs.lightred    = red_light_gradients{5};
plt_clrs.darkred     = red_dark_gradients{3}; 

camp_gr      =  plt_clrs.lightgray; 
wt_camp_gr   =  'k';

glu_color    = plt_clrs.lightgray; 
wt_camp_gl   =  'k';

eth_color    = plt_clrs.lightgreen; 
wt_camp_eh   =  plt_clrs.darkgreen;

camp_cell    = plt_clrs.lightorange;
wt_camp_cell =  plt_clrs.darkorange;

fru_color    = plt_clrs.lightblue ; 
wt_camp_fr   = plt_clrs.darkblue; 

suc_color    = plt_clrs.lightred; 
wt_camp_su   = plt_clrs.darkred;

plt_clrs.gl_colors   = plt_clrs.gray;
plt_clrs.eh_colors   = plt_clrs.green;
plt_clrs.cell_colors = plt_clrs.yellow;

wt_maker_size  = 10;
wt_maker_shape = 'd'; 