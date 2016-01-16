function [W,e,optimum] = scalar_lcmv(L,C1_inverse,normalize_leadfield,normalize_weight,objective,C0_inverse,e)

% Scalar LCMV adaptive spatial filter for dipolar source AT ONE FIXED LOCATION

% L = lead field for the given source; columns correpond to the spatial responses for
%     orthogonal source orientations (usually 2 or 3 columns); Nchannel x Norientation
% C1_inverse = inverse covariance matrix for sensor signals; Nchannel x Nchannel

% normalize_leadfield = if true, normalize combined lead field to reduce spatial bias
% normalize_weight = if true, normalize W to reduce noise interference

% objective = if orientation e is not specified, optimize orientation w.r.t.
%             given objective -- 'power' or 'NAI' (neural activity index)
% C0_inverse = inverse covariance matrix for sensor noise, required for calculating 
%              NAI (if not specified, identity matrix is used); Nchannel x Nchannel

% W = spatial filter weights for the optimal or input orientation
% e = optimal or input orientation (normalized)
% optimum = power or NAI for the optimal or input orientation

% To obtain source timecourse use: W'*rawdata


if nargin < 3
    normalize_leadfield = false;
end
if nargin < 4
    normalize_weight = false;
end
if nargin < 5
    objective = 'power';
end
if nargin < 6 && (strcmpi(objective,'nai') || strcmpi(objective,'nai-1') || strcmpi(objective,'nai-3'))
    C0_inverse = eye(size(C1_inverse));
end


if nargin > 6 && numel(e)==size(L,2) && norm(e,'fro')>eps
    e = e(:)/norm(e,'fro');
    optimum = NaN;

else % optimize orientation
    LC1inv = L'*C1_inverse;
    
    if strcmpi(objective,'power')
        if normalize_weight
            [e,D] = geneig_weights_scores(LC1inv*L, LC1inv*LC1inv');
        elseif normalize_leadfield
            [e,D] = geneig_weights_scores(L'*L, LC1inv*L);
        else
            [e,D] = eig(LC1inv*L);
            [D,ord] = sort(1./diag(D),'descend');
            e = e(:,ord);
        end
    else
        % optimize w.r.t. type 1, type2, or type3 neural activity index (NAI)
        % cf. Huang et al. Brain Topogr 16 (2004): 139-158
        % Note that NAI-2 is the same as power of weight-normalized spatial filter

        switch lower(objective)
            case {'nai', 'nai-1'}
                [e,D] = geneig_weights_scores(L'*C0_inverse*L, LC1inv*L);
            case 'nai-2'
                [e,D] = geneig_weights_scores(LC1inv*L, LC1inv*LC1inv');
            case 'nai-3'
                [e,D] = geneig_weights_scores(LC1inv*L, (LC1inv/C0_inverse)*LC1inv');
            otherwise
                error('objective should be one of ''power'', ''NAI-1'', ''NAI-2'', ''NAI-3''');
        end
    end
    
    e = e(:,1)/norm(e(:,1));
    optimum = D(1);
end

Le = L * e;

if normalize_leadfield
    Le = Le/norm(Le);
end

W = (C1_inverse * Le) / (Le'* C1_inverse * Le);

if normalize_weight
    W = W/norm(W);
end

if ~isfinite(optimum)
    LeC1inv = Le'*C1_inverse;
    switch lower(objective)
        case 'power'
            if normalize_weight
                optimum = (LeC1inv*Le) / (LeC1inv*LeC1inv');
            else
                optimum = 1 / (LeC1inv*Le);
            end
        case {'nai', 'nai-1'}
            optimum = (Le'*C0_inverse*Le) / (LeC1inv*Le);
        case 'nai-2'
            optimum = (LeC1inv*Le) / (LeC1inv*LeC1inv');
        case 'nai-3'
            optimum = (LeC1inv*Le) / ((LeC1inv/C0_inverse)*LeC1inv');
        otherwise
            error('objective should be one of ''power'', ''NAI-1'', ''NAI-2'', ''NAI-3''');
    end
end

end
