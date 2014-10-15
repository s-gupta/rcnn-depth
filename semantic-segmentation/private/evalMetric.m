function ap = evalMetric(gt, out, wt, typ)
	switch typ,
		case 'ap1',
			%gt is +1 or -1
			[P R ap acc maxAcc cM] = calcPR((gt+1)/2, out, wt);
		case 'ap1min',
			%gt is +1 or -1
			[P R ap1 acc maxAcc cM] = calcPR((gt+1)/2, out, wt);
			[P R ap2 acc maxAcc cM] = calcPR((-gt+1)/2, -out, wt);
			ap = min(ap1, ap2);
		case 'multiClassAccuracy',
			[gr, pred] = max(out,[],1);
			ap = sum((gt == pred) .* wt)./sum(wt);
		otherwise,
			ap = 0;
	end
end
