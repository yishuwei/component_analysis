function W = vector_lcmv(L,C1_inverse,normalize_leadfield,normalize_weight,ortho_constraint)

% Vector LCMV adaptive spatial filter for dipolar source AT ONE FIXED LOCATION

% L = lead field for the given source; columns correpond to the spatial responses for
%     orthogonal source orientations (usually 2 or 3 columns); Nchannel x Norientation
% C1_inverse = inverse covariance matrix for sensor signals; Nchannel x Nchannel

% normalize_leadfield = if true, normalize each column of L to reduce spatial bias
% normalize_weight = if true, normalize each column of W to reduce noise interference
% ortho_constraint = if true, enforce the orthogonal constraint W * L = I
%                    cf. Johnson et al. PLoS One 6 (2011): e22251

% W = spatial filter weights for each orientation; Nchannel x Norientation

% To obtain source timecourses use: W'*rawdata


if nargin < 3
    normalize_leadfield = false;
end
if nargin < 4
    normalize_weight = false;
end
if nargin < 5
    ortho_constraint = false;
end


if normalize_leadfield
    for j = 1:size(L,2)
        L(:,j) = L(:,j)/norm(L(:,j));
    end
end


if ortho_constraint
    W = (C1_inverse * L) / (L'* C1_inverse * L);
else
    W = zeros(size(L));
    for j = 1:size(L,2)
        W(:,j) =  (C1_inverse * L(:,j)) / (L(:,j)'* C1_inverse * L(:,j));
    end
end


if normalize_weight
    for j = 1:size(W,2)
        W(:,j) = W(:,j)/norm(W(:,j));
    end
end

end
