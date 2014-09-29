% Copyright (c) 2012, Jonathan Barron and Jitendra Malik (UC Berkeley)
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of NYU nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY JONATHAN BARRON AND JITENDRA MALIK ``AS IS''
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
% THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


function [X_white, params] = whiten(X, DO_CENTER, V_PAD)
%     [X_white, params] = whiten(X, DO_CENTER)
%
%     X = (X_white * params.inverse) + repmat(params.mean, size(X_white,1),1)
%     X_white = (X - repmat(params.mean, size(X,1),1)) * params.map

if nargin < 2
    DO_CENTER = 1;
end

if nargin < 3
    V_PAD = .1;
end

N_cov = 100000;

X = double(X);
m = mean(X);
if ~DO_CENTER
    m = 0*m;
end
X_zeroed = X - ones(size(X,1),1)*m;

rng('default');
%s = rand('twister');
X_sub = X_zeroed(rand(size(X_zeroed,1),1) <= N_cov/size(X_zeroed,1),:);
C = (X_sub' * X_sub) / size(X_sub,1);
%     C = cov(X_zeroed(rand(size(X_zeroed,1),1) <= N_cov/size(X_zeroed,1),:));
%rand('twister', s);

%     while true
%       [V,D] = eig(C);
%       valid = all(D>=0, 1);% & ~any(imag(D),1);
%       if(all(valid))
%         break
%       end
%       fprintf('WHITEN: warning, poorly behaved covariance, correcting...\n');
%       D(~valid,~valid) = -10*D(~valid,~valid);
%       C = V * D * inv(V);
%       C = (C + C')/2;
%     end

[V,D] = eig(C);

iD = diag(sqrt(1./(diag(D) + V_PAD)));
map = V * iD * V';

inverse = inv(map);

X_white = X_zeroed * map;

params.mean = m;
params.map = map;
params.inverse = inverse;
params.V = V;
params.iD = iD;
params.D = D;
params.C = C;
%params.iC = inv(C);

mag = sqrt(mean(X_white(:).^2));
params.map = params.map / mag;
params.inverse = params.inverse * mag;
X_white = X_white ./ mag;

end