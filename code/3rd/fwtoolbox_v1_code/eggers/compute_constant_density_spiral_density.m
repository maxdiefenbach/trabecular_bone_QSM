% Compute constant density spiral density
function [weight] = compute_constant_density_spiral_density( nv, ns, ms, is, lambda, delay, trajectory )

  c1 =  is * sqrt( 1.0 + lambda ) / ns;
  c2 =  lambda / ns;
  c3 =  0.5 * c1 / is;
  c1 =  pi / nv * c1;

  ir =  floor( delay );
  fr =  delay - ir;

  weight =  zeros( ns, 1 );

  for ii=ir:ns-1 

    k  =  fr + ii;

    c4 =  c3 * (1.0 / sqrt( 1.0 + c2 * k ) - 0.5 * c2 * k / (sqrt( 1.0 + c2 * k ))^3);
    
    c5 =  c1 * k / sqrt( 1.0 + c2 * k );

    weight(ii+1) =  c4 * ((cos( c5 ) - sin( c5 )) * trajectory(2*ii+1) + (sin( c5 ) + cos( c5 )) * trajectory(2*ii+2));

  end

  weight(1) =  0.125 * weight(1);

  c =  (pi * ms * ms) / (4.0 * nv * sum( weight ));

  weight =  c * weight;

end
