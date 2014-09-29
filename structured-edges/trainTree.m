function dummyOut = trainTree( opts, stream, treeInd )
dummyOut = 0;
% Train a single tree in forest model.
if(opts.rgbd == 3), colorModel = load(opts.rgbd3opts.colorModelFile);end

%% EDIT BY SG
if(~isempty(opts.train_data))
  imgIds = opts.train_data.imlist;
  nImgs = length(imgIds);
  trnImgDir = opts.train_data.image_dir;
  trnDepDir = opts.train_data.dep_dir;
  trnGtDir = opts.train_data.ground_truth_dir;
  ext = opts.train_data.ext;
else
  trnImgDir = [opts.bsdsDir '/images/train/'];
  trnDepDir = [opts.bsdsDir '/depth/train/'];
  trnGtDir = [opts.bsdsDir '/groundTruth/train/'];
  imgIds=dir(trnImgDir); imgIds=imgIds([imgIds.bytes]>0);
  imgIds={imgIds.name}; ext=imgIds{1}(end-2:end);
  nImgs=length(imgIds); for i=1:nImgs, imgIds{i}=imgIds{i}(1:end-4); end
end
%% END EDIT BY SG

% extract commonly used options
imWidth=opts.imWidth; imRadius=imWidth/2;
gtWidth=opts.gtWidth; gtRadius=gtWidth/2;
nChns=opts.nChns; nTotFtrs=opts.nTotFtrs; rgbd=opts.rgbd;
nPos=opts.nPos; nNeg=opts.nNeg; shrink=opts.shrink;

% finalize setup
treeDir = [opts.modelDir '/tree/'];
treeFn = [treeDir opts.modelFnm '_tree'];
if(exist([treeFn int2str2(treeInd,3) '.mat'],'file'))
  fprintf('Reusing tree %d of %d\n',treeInd,opts.nTrees); return; end
fprintf('\n-------------------------------------------\n');
fprintf('Training tree %d of %d\n',treeInd,opts.nTrees); tStart=clock;

% set global stream to stream with given substream (will undo at end)
streamOrig = RandStream.getGlobalStream();
set(stream,'Substream',treeInd);
RandStream.setGlobalStream( stream );

% collect positive and negative patches and compute features
fids=sort(randperm(nTotFtrs,round(nTotFtrs*opts.fracFtrs)));
k = nPos+nNeg; nImgs=min(nImgs,opts.nImgs);
ftrs = zeros(k,length(fids),'single');
labels = zeros(gtWidth,gtWidth,k,'uint8'); k = 0;
tid = ticStatus('Collecting data',30,1);
for i = 1:nImgs
  % get image and compute channels
  gt=load([trnGtDir imgIds{i} '.mat']); gt=gt.groundTruth;
  Iin=imread([trnImgDir imgIds{i} '.' ext]); 
  I = Iin; siz=size(I);
  if(rgbd), 
    Din = imread([trnDepDir imgIds{i} '.png']);
    D=single(Din)/1e4; 
  end
  
  %% BEGIN EDIT BY SG
  ng = zeros(size(I,1), size(I,2), 0);

  if(rgbd==1), 
    I=D; 
  elseif(rgbd==2), 
    I=cat(3,single(I)/255,D); 
  elseif(rgbd==3), 
    I=cat(3,single(I)/255,D,1e-3./D); 
    cues = getAllCues(Iin, Din, opts.rgbd3opts.C, colorModel, opts.rgbd3opts.vars, ...
      fullfile_ext(opts.rgbd3opts.cueCacheDir, imgIds{i}, 'mat')); 
    ng = cat(3, cues{:});
  end

  p=zeros(1,4); p([2 4])=mod(4-mod(siz(1:2),4),4);
  if(any(p)), 
    I=imPad(I,p,'symmetric'); 
    ng=imPad(ng,p,'symmetric');
  end
  [chnsReg,chnsSim] = edgesChns(I, opts, ng);
  % END EDIT BY SG
  %%

  % sample positive and negative locations
  nGt=length(gt); xy=[]; k1=0; B=false(siz(1),siz(2));
  B(shrink:shrink:end,shrink:shrink:end)=1;
  B([1:imRadius end-imRadius:end],:)=0;
  B(:,[1:imRadius end-imRadius:end])=0;
  for j=1:nGt
    M=gt{j}.Boundaries; M(bwdist(M)<gtRadius)=1;
    [y,x]=find(M.*B); k2=min(length(y),ceil(nPos/nImgs/nGt));
    rp=randperm(length(y),k2); y=y(rp); x=x(rp);
    xy=[xy; x y ones(k2,1)*j]; k1=k1+k2; %#ok<AGROW>
    [y,x]=find(~M.*B); k2=min(length(y),ceil(nNeg/nImgs/nGt));
    rp=randperm(length(y),k2); y=y(rp); x=x(rp);
    xy=[xy; x y ones(k2,1)*j]; k1=k1+k2; %#ok<AGROW>
  end
  if(k1>size(ftrs,1)-k), k1=size(ftrs,1)-k; xy=xy(1:k1,:); end
  % crop patches and ground truth labels
  psReg=zeros(imWidth/shrink,imWidth/shrink,nChns,k1,'single');
  lbls=zeros(gtWidth,gtWidth,k1,'uint8');
  psSim=psReg; ri=imRadius/shrink; rg=gtRadius;
  for j=1:k1, xy1=xy(j,:); xy2=xy1/shrink;
    psReg(:,:,:,j)=chnsReg(xy2(2)-ri+1:xy2(2)+ri,xy2(1)-ri+1:xy2(1)+ri,:);
    psSim(:,:,:,j)=chnsSim(xy2(2)-ri+1:xy2(2)+ri,xy2(1)-ri+1:xy2(1)+ri,:);
    t=gt{xy1(3)}.Segmentation(xy1(2)-rg+1:xy1(2)+rg,xy1(1)-rg+1:xy1(1)+rg);
    if(all(t(:)==t(1))), lbls(:,:,j)=1; else [~,~,t]=unique(t);
      lbls(:,:,j)=reshape(t,gtWidth,gtWidth); end
  end
  if(0), figure(1); montage2(squeeze(psReg(:,:,1,:))); drawnow; end
  if(0), figure(2); montage2(lbls(:,:,:)); drawnow; end
  % compute features and store
  ftrs1=[reshape(psReg,[],k1)' stComputeSimFtrs(psSim,opts)];
  ftrs(k+1:k+k1,:)=ftrs1(:,fids); labels(:,:,k+1:k+k1)=lbls;
  k=k+k1; if(k==size(ftrs,1)), tocStatus(tid,1); break; end
  tocStatus(tid,i/nImgs);
end
if(k<size(ftrs,1)), ftrs=ftrs(1:k,:); labels=labels(:,:,1:k); end

% train structured edge classifier (random decision tree)
pTree=struct('minCount',opts.minCount, 'minChild',opts.minChild, ...
  'maxDepth',opts.maxDepth, 'H',opts.nClasses, 'split',opts.split);
t=labels; labels=cell(k,1); for i=1:k, labels{i}=t(:,:,i); end
pTree.discretize=@(hs,H) discretize(hs,H,opts.nSamples,opts.discretize);
tree=forestTrain(ftrs,labels,pTree); tree.hs=cell2array(tree.hs);
tree.fids(tree.child>0) = fids(tree.fids(tree.child>0)+1)-1;
if(~exist(treeDir,'dir')), mkdir(treeDir); end
save([treeFn int2str2(treeInd,3) '.mat'],'tree'); e=etime(clock,tStart);
fprintf('Training of tree %d complete (time=%.1fs).\n',treeInd,e);
RandStream.setGlobalStream( streamOrig );

end

function ftrs = stComputeSimFtrs( chns, opts )
% Compute self-similarity features (order must be compatible w mex file).
w=opts.imWidth/opts.shrink; n=opts.nCells; if(n==0), ftrs=[]; return; end
nSimFtrs=opts.nSimFtrs; nChns=opts.nChns; m=size(chns,4);
inds=round(w/n/2); inds=round((1:n)*(w+2*inds-1)/(n+1)-inds+1);
chns=reshape(chns(inds,inds,:,:),n*n,nChns,m);
ftrs=zeros(nSimFtrs/nChns,nChns,m,'single');
k=0; for i=1:n*n-1, k1=n*n-i; i1=ones(1,k1)*i;
  ftrs(k+1:k+k1,:,:)=chns(i1,:,:)-chns((1:k1)+i,:,:); k=k+k1; end
ftrs = reshape(ftrs,nSimFtrs,m)';
end

function [hs,segs] = discretize( segs, nClasses, nSamples, type )
% Convert a set of segmentations into a set of labels in [1,nClasses].
persistent cache; w=size(segs{1},1); assert(size(segs{1},2)==w);
if(~isempty(cache) && cache{1}==w), [~,is1,is2]=deal(cache{:}); else
  % compute all possible lookup inds for w x w patches
  is=1:w^4; is1=floor((is-1)/w/w); is2=is-is1*w*w; is1=is1+1;
  kp=is2>is1; is1=is1(kp); is2=is2(kp); cache={w,is1,is2};
end
% compute n binary codes zs of length nSamples
nSamples=min(nSamples,length(is1)); kp=randperm(length(is1),nSamples);
n=length(segs); is1=is1(kp); is2=is2(kp); zs=false(n,nSamples);
for i=1:n, zs(i,:)=segs{i}(is1)==segs{i}(is2); end
zs=bsxfun(@minus,zs,sum(zs,1)/n); zs=zs(:,any(zs,1));
if(isempty(zs)), hs=ones(n,1,'uint32'); segs=segs{1}; return; end
% find most representative segs (closest to mean)
[~,ind]=min(sum(zs.*zs,2)); segs=segs{ind};
% apply PCA to reduce dimensionality of zs
U=pca(zs'); d=min(5,size(U,2)); zs=zs*U(:,1:d);
% discretize zs by clustering or discretizing pca dimensions
d=min(d,floor(log2(nClasses))); hs=zeros(n,1);
for i=1:d, hs=hs+(zs(:,i)<0)*2^(i-1); end
[~,~,hs]=unique(hs); hs=uint32(hs);
if(strcmpi(type,'kmeans'))
  nClasses1=max(hs); C=zs(1:nClasses1,:);
  for i=1:nClasses1, C(i,:)=mean(zs(hs==i,:),1); end
  hs=uint32(kmeans2(zs,nClasses,'C0',C,'nIter',1));
end
% optionally display different types of hs
for i=1:0, figure(i); montage2(cell2array(segs(hs==i))); end
end
