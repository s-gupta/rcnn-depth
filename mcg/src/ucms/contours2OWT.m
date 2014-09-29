% ------------------------------------------------------------------------ 
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
%  University of California Berkeley (UCB) - USA
% 
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  Pablo Arbelaez <arbelaez@berkeley.edu>
%  June 2014
% ------------------------------------------------------------------------ 
% This file is part of the MCG package presented in:
%    Arbelaez P, Pont-Tuset J, Barron J, Marques F, Malik J,
%    "Multiscale Combinatorial Grouping,"
%    Computer Vision and Pattern Recognition (CVPR) 2014.
% Please consider citing the paper if you use this code.
% ------------------------------------------------------------------------

function [owt2, superpixels] = contours2OWT(pb_cont,pb_ori)
% oriented watershed transform
%
% Pablo Arbelaez <arbelaez@berkeley.edu>

pb_cont(pb_cont<0)=0; % requires positive signal.

% create finest partition and transfer contour strength
ws_wt = create_finest_partition(pb_cont, pb_ori);

% convert to ucm2 format
owt2 = double(super_contour_4c(ws_wt));
owt2 = clean_watersheds(owt2);
owt2(end+1, :) = owt2(end, :);
owt2(:, end+1) = owt2(:, end);

% pre-compute finest superpixels
 labels2 = bwlabel(owt2 == 0, 8);
 superpixels = labels2(2:2:end, 2:2:end) - 1; % labels begin at 0 in mex file.

 %%
function ws_wt = create_finest_partition(pb_cont,pb_ori)

ws = watershed(pb_cont);
ws_bw = (ws == 0);

contours = fit_contour(double(ws_bw));
angles = zeros(numel(contours.edge_x_coords), 1);

for e = 1 : numel(contours.edge_x_coords)
    if contours.is_completion(e), continue; end
    v1 = contours.vertices(contours.edges(e, 1), :);
    v2 = contours.vertices(contours.edges(e, 2), :);
    
    if v1(2) == v2(2),
        ang = pi/2;
    else
        ang = atan((v1(1)-v2(1)) / (v1(2)-v2(2)));
    end
    angles(e) = ang*180/pi;
end

ws_wt = zeros(size(ws_bw));
for e = 1 : numel(contours.edge_x_coords)
    if contours.is_completion(e), continue; end
    for p = 1 : numel(contours.edge_x_coords{e}),
        vl = pb_cont(contours.edge_x_coords{e}(p), contours.edge_y_coords{e}(p))*abs(cos( pb_ori(contours.edge_x_coords{e}(p), contours.edge_y_coords{e}(p)) - pi/2 - angles(e)*pi/180 ));
        if ws_wt(contours.edge_x_coords{e}(p), contours.edge_y_coords{e}(p))<vl,
            ws_wt(contours.edge_x_coords{e}(p), contours.edge_y_coords{e}(p)) = vl;
        end
    end
end

[tx, ty] = size(ws_wt);
% assign junctions to strongest neighbor
for e = 1 : numel(contours.edge_x_coords)
    if contours.is_completion(e), continue; end
    
    v1 = contours.vertices(contours.edges(e,1),:);
    if (v1(1)<tx) && (ws_wt(v1(1)+1, v1(2)) > ws_wt(v1(1), v1(2))),
        ws_wt(v1(1), v1(2)) = ws_wt(v1(1)+1, v1(2));
    end
    if (v1(1)>1) && (ws_wt(v1(1)-1, v1(2)) > ws_wt(v1(1), v1(2))),
        ws_wt(v1(1), v1(2)) = ws_wt(v1(1)-1, v1(2));
    end
    if (v1(2)<ty) && (ws_wt(v1(1), v1(2)+1) > ws_wt(v1(1), v1(2))),
        ws_wt(v1(1), v1(2)) = ws_wt(v1(1), v1(2)+1);
    end
    if (v1(2)>1) && (ws_wt(v1(1), v1(2)-1) > ws_wt(v1(1), v1(2))),
        ws_wt(v1(1), v1(2)) = ws_wt(v1(1), v1(2)-1);
    end
    
    v2 = contours.vertices(contours.edges(e,2),:);
    if (v2(1)<tx) && (ws_wt(v2(1)+1, v2(2)) > ws_wt(v2(1), v2(2))),
        ws_wt(v2(1), v2(2)) = ws_wt(v2(1)+1, v2(2));
    end
    if (v2(1)>1) && (ws_wt(v2(1)-1, v2(2)) > ws_wt(v2(1), v2(2))),
        ws_wt(v2(1), v2(2)) = ws_wt(v2(1)-1, v2(2));
    end
    if (v2(2)<ty) && (ws_wt(v2(1), v2(2)+1) > ws_wt(v2(1), v2(2))),
        ws_wt(v2(1), v2(2)) = ws_wt(v2(1), v2(2)+1);
    end
    if (v2(2)>1) && (ws_wt(v2(1), v2(2)-1) > ws_wt(v2(1), v2(2))),
        ws_wt(v2(1), v2(2)) = ws_wt(v2(1), v2(2)-1);
    end
end

%%
function [pb2, V, H] = super_contour_4c(pb)

V = min(pb(1:end-1,:), pb(2:end,:));
H = min(pb(:,1:end-1), pb(:,2:end));

[tx, ty] = size(pb);
pb2 = zeros(2*tx, 2*ty);
pb2(1:2:end, 1:2:end) = pb;
pb2(1:2:end, 2:2:end-2) = H;
pb2(2:2:end-2, 1:2:end) = V;
pb2(end,:) = pb2(end-1, :);
pb2(:,end) = max(pb2(:,end), pb2(:,end-1));

%%

function [ws_clean] = clean_watersheds(ws)
% remove artifacts created by non-thin watersheds (2x2 blocks) that produce
% isolated pixels in super_contour

ws_clean = ws;

c = bwmorph(ws_clean == 0, 'clean', inf);

artifacts = ( c==0 & ws_clean==0 );
%R = regionprops(bwlabel(artifacts), 'PixelList'); % back compatibility < R2009a
R = regionprops(artifacts, 'PixelList');

for r = 1 : numel(R),
    xc = R(r).PixelList(1,2);
    yc = R(r).PixelList(1,1);
    
    vec = [ max(ws_clean(xc-2, yc-1), ws_clean(xc-1, yc-2)) ...
        max(ws_clean(xc+2, yc-1), ws_clean(xc+1, yc-2)) ...
        max(ws_clean(xc+2, yc+1), ws_clean(xc+1, yc+2)) ...
        max(ws_clean(xc-2, yc+1), ws_clean(xc-1, yc+2)) ];
    
    [~, id] = min(vec);
    switch id,
        case 1,
            if ws_clean(xc-2, yc-1) < ws_clean(xc-1, yc-2),
                ws_clean(xc, yc-1) = 0;
                ws_clean(xc-1, yc) = vec(1);
            else
                ws_clean(xc, yc-1) = vec(1);
                ws_clean(xc-1, yc) = 0;
                
            end
            ws_clean(xc-1, yc-1) = vec(1);
        case 2,
            if ws_clean(xc+2, yc-1) < ws_clean(xc+1, yc-2),
                ws_clean(xc, yc-1) = 0;
                ws_clean(xc+1, yc) = vec(2);
            else
                ws_clean(xc, yc-1) = vec(2);
                ws_clean(xc+1, yc) = 0;
            end
            ws_clean(xc+1, yc-1) = vec(2);
            
        case 3,
            if ws_clean(xc+2, yc+1) < ws_clean(xc+1, yc+2),
                ws_clean(xc, yc+1) = 0;
                ws_clean(xc+1, yc) = vec(3);
            else
                ws_clean(xc, yc+1) = vec(3);
                ws_clean(xc+1, yc) = 0;
            end
            ws_clean(xc+1, yc+1) = vec(3);
        case 4,
            if ws_clean(xc-2, yc+1) < ws_clean(xc-1, yc+2),
                ws_clean(xc, yc+1) = 0;
                ws_clean(xc-1, yc) = vec(4);
            else
                ws_clean(xc, yc+1) = vec(4);
                ws_clean(xc-1, yc) = 0;
            end
            ws_clean(xc-1, yc+1) = vec(4);
    end
end

% extract contours and neighboring regions given non-max suppressed edge map
function contours = fit_contour(nmax)

% extract contours
[skel, labels, is_v, is_e, assign, vertices, edges, ...
    v_left, v_right, e_left, e_right, c_left, c_right, ...
    edge_equiv_ids, is_compl, e_x_coords, e_y_coords] = ...
    mex_contour_sides(nmax,true);


% store pixel assignment maps
contours.skel   = skel;
contours.labels = labels;
contours.is_v   = is_v;
contours.is_e   = is_e;
contours.assign   = assign + 1;     % adjust from 0 to 1 based indexing

% store vertices and edges
contours.vertices = vertices + 1;   % adjust from 0 to 1 based indexing
contours.edges    = edges + 1;      % adjust from 0 to 1 based indexing

% store coordinates of pixels on edges (excluding endpoints)
contours.edge_x_coords = e_x_coords;
contours.edge_y_coords = e_y_coords;
for n = 1:length(e_x_coords)
    contours.edge_x_coords{n} = contours.edge_x_coords{n} + 1; % 0 -> 1 indexing
    contours.edge_y_coords{n} = contours.edge_y_coords{n} + 1; % 0 -> 1 indexing
end

% store edge equiv ids
contours.edge_equiv_ids = edge_equiv_ids + 1;

% store completion flags
contours.is_completion = is_compl;

