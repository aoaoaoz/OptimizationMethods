tic

n = 1000;eps = 1e-1;
x0 = squeeze(ones(1, n) * 3);
X=sym('X', [1, n]);
iter = 0;
dqdrtic_gradient = gradient(dqdrtic(X), X);
dir = subs(dqdrtic_gradient, X, x0);
dir = double( dir / norm(dir))';

while true
    syms alpha
    iter = iter + 1
    x1 = x0 + dir * alpha;
    alpha = double( solve(diff(dqdrtic(x1))) );
    x0 = double(x0 + dir * alpha);
    dir = double(subs(dqdrtic_gradient, X, x0));
    grad = dir;
    dir = double( dir / norm(dir))';
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