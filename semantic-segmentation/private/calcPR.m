function [P R ap acc maxAcc cM] = calcPR(gt, out, wt)
	if(~exist('wt','var'))
		wt = ones(size(gt));
	end
	assert(all(gt >= 0),'Ground truth should be between 0 and 1.\n');
	gt = wt(:).*gt(:);

	tog = [gt(:), wt(:), out(:)];
	tog = sortrows(tog,-3);
	sortgt = tog(:,1);
	cumsumsortgt = cumsum(sortgt);
	sortwt = tog(:,2);
	cumsumsortwt = cumsum(sortwt);
	P = cumsumsortgt./cumsumsortwt;
	R = cumsumsortgt./sum(sortgt);
	ap = VOCap(R,P);


	tp = cumsumsortgt;
	fp = cumsumsortwt - cumsumsortgt;
	tn = (sum(wt)-cumsumsortwt)-(sum(gt)-cumsumsortgt);
	fn = sum(gt)-cumsumsortgt;
	acc = (tp+tn)./(tp+fp+tn+fn);
	maxAcc = max(acc);
	[gr ind] = max(tp./(tp+fn)+tn./(tn+fp));
	cM(1,:) = [tp(ind) fn(ind)]./(tp(ind)+fn(ind));
	cM(2,:) = [fp(ind) tn(ind)]./(fp(ind)+tn(ind));
end
