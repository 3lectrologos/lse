%--------------------------------------------------------------------------
% Author      : Alkis Gotovos <alkisg@inf.ethz.ch>
% Description : LSE demo
%--------------------------------------------------------------------------
xmin = -0.2;
xmax = 0.8;
nsamples = 3000;
% Sample data set from function
x = unifrnd(xmin, xmax, nsamples, 2);
y = fun_cosines(x(:, 1), x(:, 2));
% Set threshold
h = 1;
% Set up GP hyperparameters
hyp.fun.mean = @meanConst;
hyp.fun.cov = @covSEiso;
hyp.fun.lik = @likGauss;
hyp.val.mean = 0;
hyp.val.cov = [-1; 5];
hyp.val.lik = -1;
% Create testcase
tc.samples.x = x;
tc.samples.y = y;
tc.h = h;
tc.hyp = hyp;
tc.name = 'cosines';

% Print true contours
figure('Position', [200, 200, 1400, 670]);
subplot(1, 2, 1);
res = 100;
[x1, x2] = meshgrid(linspace(xmin, xmax, res), linspace(xmin, xmax, res));
f = fun_cosines(x1, x2);
contourf(x1, x2, f, [h, h]);
title('True contour');

% Run LSE on testcase
subplot(1, 2, 2);
epsilon = 0.4;
niter = 200;
obj = Lse(tc, epsilon);
obj.run(niter);
fprintf(1, 'F1-score: %.2f\n', [obj.eval()]);