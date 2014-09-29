function make_protofiles(protodir, file_name, args), 
% function make_protofiles(protodir, file_name, args) 

% AUTORIGHTS

  in_dir = 'nyud2_finetuning';
  for i = 1:length(file_name),
    copyfile(fullfile(in_dir, file_name{i}), fullfile(protodir, file_name{i}));
    for j = 1:2:length(args),
      helper(fullfile(protodir, file_name{i}), fullfile(protodir, file_name{i}), args{j}, args{j+1});
    end
  end
end


function helper(inFile, outFile, str_a, str_b)
  % Replaces all occurences of string str_a in inFile with str_b
  f = fopen(inFile);
  text = fread(f, inf, '*char')';
  fclose(f);
  text = regexprep(text, str_a, str_b);

  f = fopen(outFile, 'w');
  fwrite(f, text);
  fclose(f);
end
