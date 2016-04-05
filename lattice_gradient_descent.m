
function [x, fval] = lattice_gradient_descent(f, x0, max_iters)
%LATTICE_GRADIENT_DESCENT Gradient descent minimisation on a lattice.
%
% [x, fval] = lattice_gradient_descent(f, x0, max)
%
% Return the minimum function value fval at point x of function f starting at
% x0. Only 1D and 2D functions are supported. Max is the maximum number of
% allowed iterations, which defaults to 100.
%
% See also:
% align_images

% Input defaults and validation
d = length(x0);
if d > 2
    error('Only 1 and 2 dimensional functions are supported.');
end
if any(round(x0) ~= x0)
    error('All function inputs must be integers.');
end
if nargin < 3
    max_iters = 100;
end

switch d
    % One-dimensional functions
    case 1
        x = x0;
        bestval = f(x);
        iters = 0;
        while iters < max_iters
            left = f(x-1);
            right = f(x+1);
            grads = bestval - [left, right];
            [best_grad, best_dir] = max(grads);
            if best_grad < 0
                break;
            end
            switch best_dir
                case 1
                    bestval = left;
                    x = x - 1;
                case 2
                    bestval = right;
                    x = x + 1;
            end
            iters = iters + 1;
        end
        if iters >= max_iters
            error('Maximum number of iterations (%d) reached.', max_iters);
        end
        fval = bestval;
    % Two-dimensional functions
    case 2
        x = x0;
        bestval = f(x);
        iters = 0;
        while iters < max_iters
            left = f([x(1)-1,x(2)]);
            right = f([x(1)+1,x(2)]);
            up = f([x(1),x(2)+1]);
            down = f([x(1),x(2)-1]);
            grads = bestval - [left, right, up, down];
            [best_grad, best_dir] = max(grads);
            if best_grad < 0
                break;
            end
            switch best_dir
                case 1
                    bestval = left;
                    x = [x(1)-1,x(2)];
                case 2
                    bestval = right;
                    x = [x(1)+1,x(2)];
                case 3
                    bestval = up;
                    x = [x(1),x(2)+1];
                case 4
                    bestval = down;
                    x = [x(1),x(2)-1];
            end
            iters = iters + 1;
        end
        if iters >= max_iters
            error('Maximum number of iterations (%d) reached.', max_iters);
        end
        fval = bestval;
end
