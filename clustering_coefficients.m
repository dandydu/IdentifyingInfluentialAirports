function ccfs = clustering_coefficients(A,varargin)

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

options = struct('edge_weight', 'matrix', 'undirected', 0, 'unweighted', 0);
options = merge_options(options,varargin{:});


edge_weights = 0;
edge_weight_opt = 'matrix';

if strcmp(options.edge_weight, 'matrix')
    % do nothing if we are using the matrix weights
else
    edge_weights = 1;
    edge_weight_opt = options.edge_weight;
end

if check
    check_matlab_bgl(A,struct('sym',options.undirected));
    
    if options.unweighted ~= 1 && edge_weights~=1
        try
            check_matlab_bgl(A,struct('nodefault',1,'values',1,'noneg',1));
        catch
            if ~isempty(varargin) && ...
                isfield(varargin{1},'unweighted') && varargin{1}.('unweighted')~=1
                rethrow(lasterror);
            else
                options.unweighted = 1;
            end
        end
    elseif options.unweighted ~= 1 && edge_weights && any(edge_weight_opt < 0)
        error('matlab_bgl:invalidParameter', ...
                'the edge_weight array must be non-negative');
    end
end

if options.undirected && trans, A = A'; end

weight_arg = options.unweighted;
if ~weight_arg
    weight_arg = edge_weight_opt;
else
    weight_arg = 0;
end

ccfs=clustering_coefficients_mex(A,options.undirected,weight_arg);


