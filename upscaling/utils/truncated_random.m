function val = truncated_random(mu, sigma, lower, upper, n)
    if nargin < 5, n = 1; end
    pd = makedist('Normal','mu',mu,'sigma',sigma);
    t = truncate(pd, lower, upper);
    sample = random(t,n,1);
    if n == 1
        val = sample(1);
    else
        val = sample;
    end
end
