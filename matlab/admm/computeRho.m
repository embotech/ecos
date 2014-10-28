function rho = computeRho(H)
rows = size(H, 1);
invH = inv(H);
eigenVector = randn(rows, 1);
[eigenVector, eigenMin, iters] = power_method(rows, invH, eigenVector, 10, 1e-4);
eigenMinApprox = 1/eigenMin;
eigenMaxApprox = norm(H, 1);
eigenMin = 1/eigs(invH, 1, 'LA');
eigenMax = eigs(H, 1, 'LA');

rho = sqrt(eigenMaxApprox);
fprintf('eigenMinApprox %g eigenMaxApprox %g rho %g eigenMin %g eigenMax %g\n', eigenMinApprox, eigenMaxApprox, rho, eigenMin, eigenMax);

end
