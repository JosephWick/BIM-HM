% terms in the summation for e13
function y = getE13(x2p, x3p, n, w)
  m=1;
  sum = zeros(size(x2p));

  while true
    sterm = e13Terms(x2p,x3p, m, n, w);
    if anynan(sterm)
      break 
    end
    sum = sum + sterm;
    m = m+1;

  y = (-1/(w*sqrt(n))) .* sum;
end
