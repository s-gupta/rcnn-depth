
% Simple
object = [1 1 0
          1 1 0
          1 1 1];
ground_truth = [1 1 1
                1 1 0
                0 0 0];
[J, inters, fp, fn] = jaccard( object, ground_truth );          
assert(J==0.5)
assert(isequal(inters, [1 1 0
                        1 1 0
                        0 0 0]))
assert(isequal(fp, [0 0 0
                    0 0 0
                    1 1 1]))
assert(isequal(fn, [0 0 1
                    0 0 0
                    0 0 0]))

             
% With ignore
object = [1 1 0
          1 1 0
          1 1 1];
ground_truth = [1 1 1
                1 0 0
                0 0 0];
valid_pixels = [1 1 1
                1 0 1
                1 0 1];
[J, inters, fp, fn] = jaccard( object, ground_truth, valid_pixels);          
assert(J==0.5)
assert(isequal(inters, [1 1 0
                        1 0 0
                        0 0 0]))
assert(isequal(fp, [0 0 0
                    0 0 0
                    1 0 1]))
assert(isequal(fn, [0 0 1
                    0 0 0
                    0 0 0]))
                
                
%% Eval masks
masks = cat(3, [1 0; 1 0], [1 1; 1 1], [1 1; 0 0]);
masks = masks>0;
gt = uint8([1 1; 2 0]);
objs = unique(gt);
objs(objs==0) = [];
n_obj = length(objs);

valid_pixels = ([1 1; 1 0]>0);

[areas,int,fn] = mex_eval_masks(masks,gt,valid_pixels);
fp = repmat(areas,n_obj,1)-int;
assert(size(fp,1)==n_obj)
assert(size(fp,2)==size(masks,3))
assert(isequal(int,[1 2 2; 1 1 0]))
assert(isequal(fp,[1 1 0; 1 2 2]))
assert(isequal(areas,[2 3 2]))

assert(isequal(fn,[1 0 0; 0 0 1]))
