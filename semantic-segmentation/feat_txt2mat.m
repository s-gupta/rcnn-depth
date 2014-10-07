function [feats] = feat_txt2mat(featFile, scales, outFile)
% function [feats] = feat_txt2mat(featFile, scales, outFile)
  if(nargin < 3),
    outFile = [];
  end
    
  tmp=strfind(scales,'+');
  nb_scales = numel(tmp)+1;


  fid=fopen(featFile);
  tline = fgetl(fid);
  if ~isequal(tline,'KOEN1'), 
    error('wrong feature file!');
  end
  dim = str2double(fgetl(fid))*nb_scales;
  nb_feats = str2double(fgetl(fid))/nb_scales;

  feats.data = zeros([nb_feats,dim],'uint8');
  feats.pos = zeros([nb_feats,2],'uint16');
  cnt = 1;
  s=1; dt2=[]; 

  while 1
     tline = fgetl(fid);
     
     if ~ischar(tline),   break,   end

     f = strfind(tline,';');
     dt = str2num(tline(f(1)+1:f(2)));
     pp = str2num(tline(8:f(1)-2));
     
     if s<nb_scales,
      dt2 = [dt2 dt];
      s = s + 1;
     else
      feats.data(cnt,:) = [dt2 dt];
      feats.pos(cnt,:) = pp(2:-1:1);
      cnt=cnt+1;
      s=1; dt2=[];
     end
     %pause;
  end

  feats.scale = scales;
  fclose(fid);

  if ~isempty(outFile), 
    save(outFile, 'feats');
  end
end
