% terms in the summation for e13
function y = getE12(x2p, x3p, n, w)
  m=1;
  sum = 0.0;

  sterm = e12Terms(x2p,x3p, m, n, w);
  while m <= 5000
    sum = sum + sterm;
    m = m+1;
    sterm = e12Terms(x2p,x3p, m, n, w);

  y = 1/(2*w) + (1/w)*sum;
end
