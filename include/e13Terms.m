% terms in the summation for e13
function y = e13Terms(x2p, x3p, m, n, w)
  t1 = sinh( (m*pi*(1-x3p)) / (w*sqrt(n)) );
  t2 = sin( (m*pi*x2p) / w );
  t3 = cosh( (m*pi) / (w*sqrt(n)) );

  y = (t1.*t2)./t3;
end
