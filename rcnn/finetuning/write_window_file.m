function write_window_file(imdirs, imexts, imlist, channels, list, window_file)
% function write_window_file(impaths, imexts, imlist, channels, list, window_file)
% Write the window data file
  
  fprintf('Writing to %s\n', window_file);
  fid = fopen(fullfile_ext(window_file, 'txt'), 'wt'); img_dir = 'b';
  for i = 1:length(imlist),
    tic_toc_print('make_window_file [%04d / %04d]\n', i, length(imlist));
    fprintf(fid, '# %d\n', i-1);
    
    for j = 1:length(imdirs),
      fprintf(fid, '%s/%s.%s\n', imdirs{j}, imlist{i}, imexts{j});
    end

    fprintf(fid, '%d\n%d\n%d\n', channels, 0, 0);
    fprintf(fid, '%d', size(list{i},1));
    str = '';
    for j = 1:size(list{i},1),
      str = sprintf('%s\n%s', str, strtrim(sprintf('%g ', list{i}(j,:)))); 
    end
    fprintf(fid, '%s\n', str);
  end
  fclose(fid);
  save(fullfile_ext(window_file, 'mat'), 'imdirs', 'imexts', 'imlist', 'list', 'channels');
end
