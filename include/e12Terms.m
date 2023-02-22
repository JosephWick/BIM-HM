% terms in the summation for e12
function y = e12Terms(x2p, x3p, m, n, w)
  t1 = cosh( (m.*.pi.*(1-x3p)) ./ (w.*sqrt(n)) );
  t2 = cos( (m.*pi.*x2p) ./ w );
  t3 = cosh( (m.*pi) ./ (w.*sqrt(n)) );

  y = (t1.*t2)./t3;
end
