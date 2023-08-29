function patch_background(x_left,x_middle, y_top, x_right, y_bottom, color1, color2, alpha1, alpha2)


% Function description: Add background color of selected area to a figure;    
% INPUT:
    % x_middle: Fill area1 from 0 to x_middle; y_bottom to y_top
    % y_bottom: Fill area1 from 0 to x_middle; y_bottom to y_top
    % y_top:    Fill area1 from 0 to x_middle; y_bottom to y_top
    % color1:   Fill area1 with color1; Set as white if nothing need to be filled
    % alpha1:   Fill area1 with transparency alpha1; 
    % x_right:  Fill area2 from x_middle to x_right; y_bottom to y_top
    % color2:   Fill area1 with color2; Set as white if nothing need to be filled
    % alpha2:   Fill area2 with transparency alpha1;   
    
    patch_x  = [x_left,   x_middle,   x_middle,    x_left];
    patch_y  = [y_bottom, y_bottom,   y_top,       y_top];
    p(1) = patch(patch_x, patch_y,'k');
    p(1).FaceColor = color1;
    p(1).FaceAlpha = alpha1;
    p(1).EdgeAlpha = 0;
    curr_ax = gca; 
    if strcmp(curr_ax.YAxisLocation,'left')
        uistack(p(1),'bottom');
    end
    patch_x  = [x_middle, x_right,    x_right,     x_middle];
    patch_y  = [y_bottom, y_bottom,   y_top,       y_top];
    p(2) = patch(patch_x, patch_y,'k');
    p(2).FaceColor = color2;
    p(2).FaceAlpha = alpha2;
    p(2).EdgeAlpha = 0;
    if strcmp(curr_ax.YAxisLocation,'left')
        uistack(p(2),'bottom');
    end
    
    
%{    
% Function description: Add background color of selected area to a figure;    
% INPUT:
    % x_middle: Fill area1 from 0 to x_middle; y_bottom to y_top
    % y_bottom: Fill area1 from 0 to x_middle; y_bottom to y_top
    % y_top:    Fill area1 from 0 to x_middle; y_bottom to y_top
    % color1:   Fill area1 with color1; Set as white if nothing need to be filled
    % alpha1:   Fill area1 with transparency alpha1; 
    % x_right:  Fill area2 from x_middle to x_right; y_bottom to y_top
    % color2:   Fill area1 with color2; Set as white if nothing need to be filled
    % alpha2:   Fill area2 with transparency alpha1;   
    
    patch_x  = [0,        x_middle,   x_middle,    0];
    patch_y  = [y_bottom, y_bottom,   y_top,       y_top];
    p(1) = patch(patch_x, patch_y,'k');
    p(1).FaceColor = color1;
    p(1).FaceAlpha = alpha1;
    p(1).EdgeAlpha = alpha1;
    uistack(p(1),'bottom');
    
    patch_x  = [x_middle, x_right,    x_right,     x_middle];
    patch_y  = [y_bottom, y_bottom,   y_top,       y_top];
    p(2) = patch(patch_x, patch_y,'k');
    p(2).FaceColor = color2;
    p(2).FaceAlpha = alpha2;
    p(2).EdgeAlpha = alpha2;
    uistack(p(2),'bottom');
%}
    
end 