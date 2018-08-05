% Conjugate phase reconstruction
function outParams = cpr( imDataParams, algoParams )
  
  [validParams, algoParams] = checkParamsAndSetDefaults_cpr( imDataParams,algoParams );

  if validParams == 0

    disp('Exiting -- data not processed');
    outParams = [];
    return;

  end
  
  outParams.species = algoParams.species; %DH*: include more info in the output
  
  
  % Derived parameter
  is                   =  size( imDataParams.species(1).amps, 1 );

  nv                   =  algoParams.nv;
  ns                   =  algoParams.ns;

  alpha                =  algoParams.alpha;

  ks                   =  algoParams.ks;
  gti                  =  algoParams.gti;

  tau                  =  algoParams.tau;
  shutter              =  algoParams.shutter;

  trajectory           =  algoParams.trajectory;
  density_compensation =  algoParams.density_compensation;


  bs =  2 * ceil( 0.5 * ks );
  ms =  2 * ceil( alpha * is / 2 );

  alpha =  ms / is;
  mu    =  bs / 2;

  gridding_table      =  compute_gridding_table( ms, is, bs, ks, nv, ns, gti, trajectory );
  gridding_kernel     =  compute_gridding_kernel( gti*ks+1, ks, ms, is );
  gridding_correction =  compute_gridding_correction( 0, shutter, ks, ms, is );

  % Perform transformation to k-space
  
  for n=1:2

    image = imDataParams.species(n).amps;

    fm    = zeros( is*is, 1 );

    for y=0:is-1

      fm( y*is+1:y*is+is ) =  imDataParams.species(n).fieldmap(1:is,y+1);

    end

    for y=0:is-1

      image(1:is,y+1) =  image(1:is,y+1) .* gridding_correction(y*is+1:y*is+is);

    end

    data_ci =  zeros( ms, ms );

    o = (ms-is)/2;

    for y=0:is-1

      data_ci( o+1:o+is, y+o+1 ) = image( 1:is, y+1 );

    end

    data_ck =  fftshift(fft2(fftshift(data_ci)));

    data_cek = zeros( (ms+bs)*(ms+bs), 1 );

    o = bs/2;

    for y=0:ms-1

      data_cek( (y+o)*(ms+bs)+o+1:(y+o)*(ms+bs)+o+ms ) = data_ck( 1:ms, y+1 );

    end

    data_cek =  fold_2d( 1, bs, ms, data_cek );

    data_nck =  zeros( nv*ns, 1 );

    data_nck =  grid_2d( 1, nv*ns, ms+bs, ks, gti, gridding_kernel, gridding_table, data_cek, data_nck );

    % Perform gridding (part I)

    for v=0:nv-1

      data_nck(v*ns+1:(v+1)*ns) =  density_compensation .* data_nck(v*ns+1:(v+1)*ns);

    end

    % Conjugate phase reconstruction

    tmin =  0.0;
    tmax =  (ns - 1) * tau;

    tc =  0.5 * (tmax + tmin);
    tr =  0.5 * (tmax - tmin);

    fmin =  min( fm );
    fmax =  max( fm );

    fc =  0.5 * (fmax + fmin);

    nb =  max( ceil( alpha * (fmax - fmin) * (tmax - tmin) + 2.0 * mu ), 2.0 * mu + 1.0 );

    T =  tr / (0.5 - mu / nb);

    [cpr_table1, cpr_table2, cpr_table3] = compute_cpr_tables( nb, ns, mu, ks, gti, tau, fc );
    
    cpr_correction =  compute_gridding_correction( 1, 0, 2*mu, alpha * round( nb*gti / alpha ), round( nb*gti / alpha ) );

    % Density compensation

    for v=0:nv-1

      data_nck(v*ns+1:(v+1)*ns) =  cpr_table3 .* data_nck(v*ns+1:(v+1)*ns);

    end

    % Data embedding

    data_ncek =  zeros( nb*nv*ns, 1 );

    for s=0:ns-1

      idx1 =  cpr_table1(2*s+1);
      idx2 =  cpr_table1(2*s+2);

      for v=0:nv-1

        for ik=0:ks-1

          data_ncek(((idx1+ik)*nv+v)*ns+s+1) =  gridding_kernel(idx2+ik*gti+1) * data_nck(v*ns+s+1);

        end

      end

    end

    % Transformation to image space

    image =  zeros( is*is, 1 );

    for b=1:nb-1

      idx1 =  cpr_table2(2*b+1);
      idx2 =  cpr_table2(2*b+2);

      if idx2 > 0

          data_cek =  zeros( (ms+bs)*(ms+bs), 1 );

          for v=0:nv-1
              data_cek =  grid_2d( 0, idx2, ms+bs, ks, gti, gridding_kernel, gridding_table(3*(v*ns+idx1)+1:3*(v*ns+idx1+idx2)), data_ncek((b*nv+v)*ns+idx1+1:(b*nv+v)*ns+idx1+idx2), data_cek );
          end

          data_cek =  fold_2d( 0, bs, ms, data_cek );

          data_ck  =  zeros( ms, ms );

          for y=0:ms-1

            data_ck(1:ms,y+1) = data_cek((y+bs/2)*(ms+bs)+bs/2+1:(y+bs/2)*(ms+bs)+bs/2+ms);

          end

          data_ci =  fftshift(ifft2(fftshift(data_ck)));

          o = (ms+bs-is)/2;

          for y=0:is-1

            tf =  2.0 * pi * (T * (fm(y*is+1:(y+1)*is) - fc) * b / nb + fm(y*is+1:(y+1)*is) * tc); 

            image(y*is+1:(y+1)*is) =  image(y*is+1:(y+1)*is) + exp( tf * 1i ) .* data_ci(o+1:o+is, y+o);

          end

      end

    end

    image = image .* gridding_correction;

    for xy=0:is*is-1

      c2 =  (fm(xy+1) - fc) * T;

      c3 =  cpr_correction (round( nb*gti / alpha ) / 2 + round( c2 * gti ));

      image(xy+1) =  c3 * image(xy+1);

    end

    for y=0:is-1

      outParams.species(n).amps(1:is,y+1) = image(y*is+1:(y+1)*is);

    end

  end

end

% Compute gridding table
function [table] = compute_gridding_table( ms, is, bs, ks, nv, ns, gti, trajectory )

  s  =  ms / is;

  c4 =  ms + bs;

  table =  zeros( 3*nv*ns, 1 );

  for ii=0:nv*ns-1
      
    c1 =  trajectory(2*ii+1);
    c2 =  trajectory(2*ii+2);

    if ( c1 <  0  )  c1 =  c1 + is;  end
    if ( c1 >= is )  c1 =  c1 - is;  end
    if ( c2 <  0  )  c2 =  c2 + is;  end
    if ( c2 >= is )  c2 =  c2 - is;  end

    if ( ceil( ks / 2.0 ) > floor( ks / 2.0 ) )

      c5 =  floor( s * c1 );
      c6 =  floor( s * c2 );

      c7 =  round( (0.5 - (s * c1) + (floor( s * c1 ))) * gti );
      c8 =  round( (0.5 - (s * c2) + (floor( s * c2 ))) * gti );

      if ( c7 <= 0 )  c5 =  c5 + 1; c7 =  c7 + gti;  end
      if ( c8 <= 0 )  c6 =  c6 + 1; c8 =  c8 + gti;  end

    else

      c5 =  floor( s * c1 ) + 1;
      c6 =  floor( s * c2 ) + 1;

      c7 =  round( (1.0 - (s * c1) + floor( s * c1 )) * gti );
      c8 =  round( (1.0 - (s * c2) + floor( s * c2 )) * gti );

    end

    table(3*ii+1) =  c4 * c6 + c5;
    table(3*ii+2) =  c7;
    table(3*ii+3) =  c8;

  end

end

% Compute gridding kernel
function [kernel] = compute_gridding_kernel( ns, ks, ms, is )

  b =  pi * (2.0 - is / ms);

  kernel =  zeros( ns, 1 );

  for ii=0:ns-1

    x =  - 0.5 * ks + ii * ks / (ns - 1);

    a =  0.25 * ks * ks - x * x;

    if ( a <= 0.0 )

      kernel(ii+1) =  b / pi;

    else

      kernel(ii+1) =  sinh( b * sqrt( a ) ) / (pi * sqrt( a ));

    end

  end

end

% Compute gridding correction
function [correction] = compute_gridding_correction( m, s, ks, ms, is )

  b =  pi * (2.0 - is / ms);

  c =  besseli( 0, ks / 2.0 * b );

  if ( m == 0 )

    o =  is * is / 2;

  else

    o =  0;

  end

  if ( m == 0 )

    correction =  zeros( is*is, 1 );

  else

    correction =  zeros( is, 1 );

  end

  for ii=0:is-1 

    x =  2.0 * pi * ( - 0.5 * is + ii ) / ms;

    a =  b * b - x * x;

    if ( a < 0.0 )

      correction(ii+o+1) =  0.0;

    else

      correction(ii+o+1) =  c / besseli( 0, 0.5 * ks * sqrt( a ) );

    end

  end

  if ( m == 0 )

    d  =  is / 2;

    ds =  d*d;

    for ii=0:is-1

      dy =  (ii-d) * (ii-d);
 
      for ij=0:is-1 

        dx =  (ij-d) * (ij-d) + dy;

        if ( s )

          if ( dx <= ds )

            correction(ii*is+ij+1) =  correction(ii+o+1) * correction(ij+o+1);

          else

            correction(ii*is+ij+1) =  0.0;

          end

        else

          correction(ii*is+ij+1) =  correction(ii+o+1) * correction(ij+o+1);

        end

      end

    end

    correction =  correction / (c * c);

  end

  if ( m == 1 )

    correction =  correction / c;

  end

end

% Perform gridding
function [r] = grid_2d( m, ns, ms, ks, gti, gt, idx, v, r )

  o1 =  ms-ks;
  o2 =  ks*gti;

  if ( m == 0 )

    % Source-driven processing

    ij =  1;
    ik =  1;

    for ii=0:ns-1 

      im =  idx(ik  ) + 1;

      ix =  idx(ik+1) + 1;
      iy =  idx(ik+2) + 1;

      ik =  ik + 3;

      for y=0:ks-1

        ky =  gt(iy);

        iy =  iy + gti;

        for x=0:ks-1 

          kx =  gt(ix);

          ix =  ix + gti;

          w  =  kx*ky;

          r(im) =  r(im) + w * v(ij);

          im =  im + 1;

        end

        im =  im + o1;

        ix =  ix - o2;

      end

      ij =  ij + 1;

    end

  end

  if ( m == 1 ) 

    % Target-driven processing

    ij =  1;
    ik =  1;

    for ii=0:ns-1

      im =  idx(ik  ) + 1;

      ix =  idx(ik+1) + 1;
      iy =  idx(ik+2) + 1;

      ik =  ik + 3;

      t  =  0.0;

      for y=0:ks-1

        ky =  gt(iy);

        iy =  iy + gti;

        for x=0:ks-1 

          kx =  gt(ix);

          ix =  ix + gti;

          w  =  kx*ky;

          t  =  t + w * v(im);

          im =  im + 1;

        end

        im =  im + o1;

        ix =  ix - o2;

      end

      r(ij) =  t;

      ij = ij + 1;

    end

  end

end

% Perform folding
function [v] = fold_2d( m, bs, ms, v )

  if ( m == 0 )

    % Source-driven processing

    iv1 =  0;
    iv2 =  ms;

    for ii=0:ms+bs-1 

      for ij=0:bs/2-1

        v(iv2+1) =  v(iv2+1) + v(iv1+1);

        iv1 =  iv1 + 1;
        iv2 =  iv2 + 1;

      end

      for ij=0:bs/2-1

        v(iv1+1) =  v(iv1+1) + v(iv2+1);

        iv1 =  iv1 + 1;
        iv2 =  iv2 + 1;

      end

      iv1 =  iv1 + ms;
      iv2 =  iv2 + ms;

    end

    iv1 = bs/2;
    iv2 = bs/2 + ms*(ms+bs);

    for ii=0:bs/2-1 

      for ij=0:ms-1

        v(iv2+1) =  v(iv2+1) + v(iv1+1);

        iv1 =  iv1 + 1;
        iv2 =  iv2 + 1;

      end

      iv1 =  iv1 + bs;
      iv2 =  iv2 + bs;

    end

    for ii=0:bs/2-1 

      for ij=0:ms-1

        v(iv1+1) =  v(iv1+1) + v(iv2+1);

        iv1 =  iv1 + 1;
        iv2 =  iv2 + 1;

      end

      iv1 =  iv1 + bs;
      iv2 =  iv2 + bs;

    end

  end

  if ( m == 1 )

    % Target-driven processing

    iv1 =  bs/2;
    iv2 =  bs/2 + ms*(ms+bs);

    for ii=0:bs/2-1 

      for ij=0:ms-1

        v(iv1+1) =  v(iv2+1);

        iv1 =  iv1 + 1;
        iv2 =  iv2 + 1;

      end

      iv1 =  iv1 + bs;
      iv2 =  iv2 + bs;

    end

    for ii=0:bs/2-1 

      for ij=0:ms-1

        v(iv2+1) =  v(iv1+1);

        iv1 =  iv1 + 1;
        iv2 =  iv2 + 1;

      end

      iv1 =  iv1 + bs;
      iv2 =  iv2 + bs;

    end

    iv1 =  0;
    iv2 =  ms;

    for ii=0:ms+bs-1

      for ij=0:bs/2-1

        v(iv1+1) =  v(iv2+1);

        iv1 =  iv1 + 1;
        iv2 =  iv2 + 1;

      end

      for ij=0:bs/2-1

        v(iv2+1) =  v(iv1+1);

        iv1 =  iv1 + 1;
        iv2 =  iv2 + 1;

      end
      
      iv1 =  iv1 + ms;
      iv2 =  iv2 + ms;

    end

  end

end

% Compute CPR tables
function [table1, table2, table3] = compute_cpr_tables( nb, ns, mu, ks, gti, tau, fc )

  libin =  0;
  n     =  0;

  table1 =  zeros( 2*ns  , 1 );
  table2 =  zeros( 2*nb+2, 1 );
  table3 =  zeros(   ns  , 1 );

  for ii=0:ns-1

    fbin =  ii / (ns - 1) * (nb - 2 * mu);
    ibin =  floor( fbin );

    table1(2*ii+1) =  ibin + 1;
    table1(2*ii+2) =  round( (1.0 - fbin + ibin) * gti );

    if ( ibin > libin )

      table2(2*ibin+1) =  ii - n;
      table2(2*ibin+2) =  n;

      libin =  ibin;
      n     =  1;

    else

      n =  n + 1;

    end

  end

  if ( n > 0 )

    table2(2*(ibin+1)+1) =  ii + 1 - n;
    table2(2*(ibin+1)+2) =  n;

    ibin =  ibin + 1;

  end

  for ii=ibin+1:nb-1

    table2(2*ii+1) =  ns;

  end

  for ii=nb:-1:1

    for ij=1:ks-1

      if ( (ii-ij) >= 0 )

        table2(2*ii+1) =  table2(2*ii+1) - table2(2*(ii-ij)+2);
        table2(2*ii+2) =  table2(2*ii+2) + table2(2*(ii-ij)+2);

      end

    end

  end

  for ii=0:ns-1

    t =  (ii - 0.5 * (ns - 1)) * tau;

    table3(ii+1) = exp( 2.0 * pi * 1i * fc * t );

  end

end