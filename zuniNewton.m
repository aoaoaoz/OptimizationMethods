tic

n = 1000;eps = 1e-1;
x0 = squeeze(ones(1, n) * 3);
X=sym('X', [1, n]);
iter = 0;
dqdrtic_gradient = gradient(dqdrtic(X), X);
dqdrtic_hessian = hessian(dqdrtic(X), X);
dqdrtic_hessian_inv = inv(dqdrtic_hessian);
dir = (- dqdrtic_hessian_inv * dqdrtic_gradient)';

while true
    syms alpha
    curdir = subs(dir, X, x0);
    iter = iter + 1
    x1 = x0 + alpha * curdir;
    alpha = double( solve(diff(dqdrtic(x1))) );
    x0 = x0 + curdir * alpha;
    grad = double(subs(dqdrtic_gradient, X, x0));
    norm(grad)
    if norm(grad) <= eps
        break
    end
end
x0
dqdrtic(x0)
iter
dir
toc
disp(['运行时间: ',num2str(toc)]);