% Compute constant density spiral trajectory
function [trajectory] = compute_constant_density_spiral_trajectory( nv, ns, ms, lambda, delay )

  c1 =  ms * sqrt( 1.0 + lambda ) / ns;
  c2 =  lambda / ns;
  c3 =  0.5 * c1 / ms;
  c4 =  2.0 * pi / nv;
  c1 =  pi / nv * c1;

  ir =  floor( delay );
  fr =  delay - ir;

  trajectory =  zeros( 2*nv*ns, 1 );

  ik =  1;

  for ii=0:nv-1

    for ij=0:ir-1

      trajectory(ik  ) =  0.0;
      trajectory(ik+1) =  0.0;

      ik = ik + 2;

    end

    ip  =  c4 * ii;

    for ij=ir:ns-1

      n  =  fr + ij;

      c5 =  n / sqrt( 1.0 + c2 * n );

      c6 =  c3 * c5;
      c7 =  c1 * c5 - ip;

      trajectory(ik  ) =  c6 * cos( c7 );
      trajectory(ik+1) =  c6 * sin( c7 );

      ik = ik + 2;
  
    end
  
  end
  
end
