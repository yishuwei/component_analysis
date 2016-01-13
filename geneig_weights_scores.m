function [W,D,A] = geneig_weights_scores(C1,C0,equalize)

% Generalized eigenvectors of inv(C0)*C1 sorted by scores

% C1 = covariance matrix for signals to target; Nchannel x Nchannel
% C0 = covariance matrix for signals to suppress; Nchannel x Nchannel

% W = projection/unmixing matrix (normalized)
% D = scores
% A = back-projection matrix

% To obtain component timecourses use: W'*rawdata

% WATCH OUT!!
% If C0 is rank-deficient there can be suspicious components with Inf scores

if nargin < 3
    equalize = false;
end
if equalize
    C0 = C0/trace(C0);
    C1 = C1/trace(C1);
end


[W,D] = eig(C1,C0);
[D,ord] = sort(diag(D),'descend');
W=W(:,ord);

for j=1:max(ord)
    W(:,j)=W(:,j)/norm(W(:,j));
end

A = inv(W);

end
