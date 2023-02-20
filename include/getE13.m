% terms in the summation for e13
function y = getE13(x2p, x3p, n, w)
  m=1;
  sum = 0.0;

  sterm = e13Terms(x2p,x3p, m, n, w);
  while m < 3000
    sum = sum + sterm;
    m = m+1;
    sterm = e13Terms(x2p,x3p, m, n, w);

  y = (-1/(w*n^0.5)) * sum;
end
