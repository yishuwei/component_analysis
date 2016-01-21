function [W,D,A] = geneig_weights_scores(C1,C0,equalize,prune_threshold)

% Generalized eigenvectors of inv(C0)*C1 sorted by scores

% C1 = covariance matrix for signals to target; Nchannel x Nchannel
% C0 = covariance matrix for signals to suppress; Nchannel x Nchannel

% prune_threshold = threshold used to prune the nullspace of C0; default is
%                   1e-6 * trace(C0)

% W = projection/unmixing matrix (normalized)
% D = scores
% A = back-projection matrix

% To obtain component timecourses use: W'*rawdata

if nargin < 3
    equalize = false;
end
if equalize
    C0 = C0/trace(C0);
    C1 = C1/trace(C1);
end

if nargin < 4
    prune_threshold = 1e-6 * trace(C0);
end


%% prune out the nullspace of C0
[U,S] = eig(C0);
[S,ord] = sort(diag(S),'ascend');
U = U(:,ord);

keep = cumsum(S) > prune_threshold;
U = U(:,keep);
C0 = U'* C0 * U;
C1 = U'* C1 * U;

%% generalized eigenvector decomposition
[W,D] = eig(C1,C0);
[D,ord] = sort(diag(D),'descend');
W = W(:,ord);

for j = 1:max(ord)
    W(:,j) = W(:,j)/norm(W(:,j));
end

A = W\U';
W = U * W;

end
