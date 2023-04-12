% terms in the summation for e13
function y = getE12(x2p, x3p, n, w)
  m=1;
  sum = zeros(size(x2p));

  while true
    sterm = e12Terms(x2p,x3p, m, n, w);
    if anynan(sterm) == 1
      break
    end
    sum = sum + sterm;
    m = m+1;

  y = 1/(2*w) + (1/w).*sum;
end
