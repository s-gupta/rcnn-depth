% ------------------------------------------------------------------------ 
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
% 
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  June 2013
% ------------------------------------------------------------------------ 

function [jaccards, inters, false_pos, false_neg, true_areas, obj_classes] = eval_cands(hier, cands, gt)
    % Get all objects
    obj_ids = unique(gt.object);
    obj_ids(obj_ids==0) = [];
    obj_ids(obj_ids==255) = [];
    n_objs = numel(obj_ids);

    % Store true_areas
    true_areas = zeros(n_objs,1);
    for kk=1:n_objs
        true_areas(kk) = sum(gt.object(:)==kk);
    end
            
    % Store the class of each object
    for jj=1:n_objs
        idx = find(gt.object==obj_ids(jj),1,'first');
        obj_classes(jj) = gt.class(idx); %#ok<AGROW>
    end
            
    % Get number of leaves
    leave_labels = unique(hier.leaves_part);
    n_leaves = length(leave_labels);

    % Overlap leaves with GT objects
    leave_fp   = zeros(n_objs,n_leaves);  % False positives
    leave_int  = zeros(n_objs,n_leaves);  % Intersection with GT
    hier.leaves_part(gt.object==255) = 0; % Erase 'difficult' areas
    for kk=1:n_objs
        tmp = hist(double(hier.leaves_part(:)).*(gt.object(:)==kk),(0:n_leaves));
        leave_int(kk,:) = tmp(2:end);
        tmp = hist(double(hier.leaves_part(:)).*(gt.object(:)~=kk),(0:n_leaves));
        leave_fp(kk,:) = tmp(2:end);
    end

    % Get number of regions
    n_mergings = length(hier.ms_struct);
    n_regs = n_leaves+n_mergings;

    % Propagate overlaps to all regions in the hierarchy
    reg_int = zeros(n_objs,n_regs);
    reg_fp  = zeros(n_objs,n_regs);

    % Add leaves
    reg_int(:,1:n_leaves) = leave_int;
    reg_fp(:,1:n_leaves)  = leave_fp;

    % Perform mergings (recursive computation)
    for jj=1:n_mergings
        parent   = hier.ms_struct(jj).parent;
        children = hier.ms_struct(jj).children;
        for kk=1:length(children)
            reg_int(:,parent) = reg_int(:,parent) + reg_int(:,children(kk));
            reg_fp(:,parent)  = reg_fp(:,parent) + reg_fp(:,children(kk));
        end
    end

    % Compute tp, fp for all candidates
    n_cands    = size(cands,1);
    inters     = zeros(n_objs,n_cands);
    false_pos  = zeros(n_objs,n_cands);
    false_neg  = zeros(n_objs,n_cands);
    jaccards   = zeros(n_objs,n_cands);
    if n_cands>0
        for kk=1:n_objs
            [inters(kk,:),false_pos(kk,:)] = mex_assess_one_sel(reg_int(kk,:),reg_fp(kk,:),cands);
            false_neg(kk,:) = true_areas(kk)-inters(kk,:);
            jaccards(kk,:) = inters(kk,:)./(false_pos(kk,:)+true_areas(kk));
        end
    end
end