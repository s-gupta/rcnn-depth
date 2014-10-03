function inst = makeContinuous(inst, begin)
  if(~exist('begin','var'))
    begin = 1;
  end
  u = sort(unique(inst(:)));
  if(isequal(u,begin:(length(u)+begin-1))),
    return;
  end
  t  = containers.Map(u,begin:(length(u)+begin-1));
  T = @(x)t(x);
  inst = arrayfun(T,inst);
end
