% ----- helper function (put at end of script or in roelfsema_mod.m) -----
function f = exGauss_mod(p,t)
    mu    = p(1);
    sigma = p(2);
    alpha = p(3);
    c     = p(4);
    d     = p(5);

    % cumulative Gaussian G(t;mu,sigma)
    G  = normcdf((t - mu)./sigma);

    % G(t; mu+alpha*sigma^2, sigma)
    Gs = normcdf((t - (mu + alpha*sigma^2))./sigma);

    % dissipating component solution of dm1/dt = -alpha*m1 + g(t;mu,sigma)
    m1 = exp(alpha*(mu + 0.5*alpha*sigma^2) - alpha*t) .* Gs;

    f = d*m1 + c*G;
end