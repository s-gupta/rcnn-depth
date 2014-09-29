function s = procid()
% Returns a string identifying the process.

% s1 = mfilename('fullpath');
% [s1, s2] = fileparts(s1);
% [s1, s2] = fileparts(s1);
% s = s1;

d = pwd();
i = strfind(d, filesep);
d = d(i(end)+1:end);
s = d;
