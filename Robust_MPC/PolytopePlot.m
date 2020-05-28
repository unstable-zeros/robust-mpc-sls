classdef PolytopePlot < handle
    %POLYTOPEPLOT define functions that are used for plotting.
    
    methods (Static)
       function show_convex(P, varargin)
            P_reduced = projectPolytope2Plane(P);
            switch numel(varargin)
                case 1
                    fill(P_reduced.V(:, 1), P_reduced.V(:, 2), varargin{1});
                case 3
                    fill(P_reduced.V(:, 1), P_reduced.V(:, 2), varargin{1}, varargin{2}, varargin{3});
            end
            hold on;
       end
    end
end

function P_projected = projectPolytope2Plane(P, dim)
    if nargin < 2
        dim = [1 2];
    end
    
    vert = P.V;
    x_vert = round(vert(:, dim(1)), 5);
    y_vert = round(vert(:, dim(2)), 5);
    idx = convhull(x_vert, y_vert);
    P_projected = Polyhedron([x_vert(idx), y_vert(idx)]);
end

