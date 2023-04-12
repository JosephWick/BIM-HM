% terms in the summation for e13
function y = getE12(x2p, x3p, n, w)
  m=1;
  sum = zeros(size(x2p));

  while m <= 3936
    sterm = e12Terms(x2p,x3p, m, n, w);
    sum = sum + sterm;
    m = m+1;

  y = 1/(2*w) + (1/w).*sum;
end
