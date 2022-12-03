module tpcore_fvdas_mod

"""
Subroutine do_ymist_pole1_i2d2 sets "dcy" at the Poles.

## Arguments
- `dcy::Matrix{AbstractFloat}` - OUT
- `qqu::Matrix{AbstractFloat}` - IN
- `i1_gl::Integer` - IN
- `i2_gl::Integer` - IN
- `ju1_gl::Integer` - IN
- `j2_gl::Integer` - IN
- `ilo::Integer` - IN
- `ihi::Integer` - IN
- `julo::Integer` - IN
- `jhi::Integer` - IN
- `i1::Integer` - IN
- `i2::Integer` - IN
- `ju1::Integer` - IN
- `j2::Integer` - IN

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function do_ymist_pole1_i2d2!(
	dcy::Matrix{AbstractFloat},
	qqu::Matrix{AbstractFloat},
	i1_gl::Integer,
	i2_gl::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer
)::Nothing
	i2d2 = i2_gl / 2

	r24  = 1.0 / 24.0
	
	if ju1 == ju1_gl
		for il = i1:i2d2
			tmp = ((8.0 * (qqu[il, ju1 + 2] - qqu[il, ju1])) + qqu[il + i2d2, ju1 + 1] - qqu[il, ju1 + 3]) * r24

			pmax = max(qqu[il, ju1], qqu[il, ju1 + 1], qqu[il, ju1 + 2]) - qqu[il, ju1 + 1]

			pmin = qqu(il,ju1+1) - min(qqu[il, ju1], qqu[il, ju1 + 1], qqu[il, ju1 + 2])

			dcy[il, ju1 + 1] = sign(tmp) * min(abs(tmp), pmin, pmax)
		end

		for il = (i1 + i2d2):i2
			tmp = ((8.0 * (qqu[il, ju1 + 2] - qqu[il, ju1])) + qqu[il - i2d2, ju1 + 1] - qqu[il, ju1 + 3]) * r24

			pmax = max(qqu[il, ju1], qqu[il, ju1 + 1], qqu[il, ju1 + 2]) - qqu[il, ju1 + 1]

			pmin = qqu[il, ju1 + 1] - min(qqu[il, ju1], qqu[il, ju1 + 1], qqu[il, ju1 + 2])

			dcy[il, ju1 + 1] = sign(tmp) * min(abs(tmp), pmin, pmax)
		end
	end

	if j2 == j2_gl
		for il = i1:i2d2
			tmp = ((8.0 * (qqu[il, j2] - qqu[il, j2 - 2])) + qqu[il, j2 - 3] - qqu[il+ i2d2, j2 - 1]) * r24

			pmax = max(qqu[il, j2 - 2], qqu[il, j2 - 1], qqu[il, j2]) - qqu[il, j2 - 1]

			pmin = qqu[il, j2 - 1] - min(qqu[il, j2 - 2], qqu[il, j2 - 1], qqu[il, j2])

			dcy[il, j2 - 1] = sign(tmp) * min(abs(tmp), pmin, pmax)
		end

		for il = (i1 + i2d2):i2
			tmp = ((8.0 * (qqu[il, j2] - qqu[il, j2 - 2])) + qqu[il, j2 - 3] - qqu[il - i2d2, j2 - 1]) * r24

			pmax = max(qqu[il, j2 - 2], qqu[il, j2 - 1], qqu[il, j2]) - qqu[il, j2 - 1]

			pmin = qqu[il, j2 - 1] - min(qqu[il, j2 - 2], qqu[il, j2 - 1], qqu[il, j2])

			dcy[il, j2 - 1] = sign(tmp) * min(abs(tmp), pmin, pmax)
		end
	end
end

"""
Subroutine do_ymist_pole2_i2d2 sets "dcy" at the Poles.

## Arguments
- `dcy::Matrix{AbstractFloat}` - OUT
- `qqu::Matrix{AbstractFloat}` - IN
- `i1_gl::Integer` - IN
- `i2_gl::Integer` - IN
- `ju1_gl::Integer` - IN
- `j2_gl::Integer` - IN
- `j1p::Integer` - IN
- `ilo::Integer` - IN
- `ihi::Integer` - IN
- `julo::Integer` - IN
- `jhi::Integer` - IN
- `i1::Integer` - IN
- `i2::Integer` - IN
- `ju1::Integer` - IN
- `j2::Integer` - IN

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function do_ymist_pole2_i2d2!(
	dcy::Matrix{AbstractFloat},
	qqu::Matrix{AbstractFloat},
	i1_gl::Integer,
	i2_gl::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	j1p::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer
)::Nothing
	i2d2 = i2_gl / 2
	
	if ju1 == ju1_gl
		if j1p /= ju1_gl + 1
			dcy[i1:i2, ju1] .= 0.0
		else
			# Determine slope in South Polar cap for scalars.
			for il = i1:i2d2
				tmp = 0.25 * (qqu[il, ju1 + 1] - qqu[il + i2d2, ju1 + 1])

				pmax = max(qqu[il, ju1 + 1], qqu[il, ju1], qqu[il + i2d2, ju1 + 1]) - qqu[il, ju1]

				pmin = qqu[il, ju1] - min(qqu[il, ju1 + 1], qqu[il, ju1], qqu[il + i2d2, ju1 + 1])
							
				dcy[il, ju1] = sign(tmp) * min(abs(tmp), pmax, pmin)
			end

			for il = (i1 + i2d2):i2
				dcy[il, ju1] = -dcy[il - i2d2, ju1]
			end
		end
	end
		
	if j2 == j2_gl

		# TODO: /=
		if j1p /= ju1_gl + 1
			dcy[i1:i2, j2] .= 0.0
		else
			# Determine slope in North Polar cap for scalars.

			for il = i1:i2d2
				tmp  = 0.25 * (qqu[il + i2d2, j2 - 1] - qqu[il, j2 - 1])

				pmax = max(qqu[il + i2d2, j2 - 1], qqu[il, j2], qqu[il, j2 - 1]) - qqu[il, j2]

				pmin = qqu[il, j2] - min(qqu[il + i2d2, j2 - 1], qqu[il, j2], qqu[il, j2 - 1])

				dcy[il, j2] = sign(tmp) * min(abs(tmp), pmax, pmin)
			end

			for il = (i1 + i2d2):i2
				dcy[il, j2] = -dcy[il - i2d2, j2]
			end
			
		end
	end
end

"""
Subroutine fyppm is the 1D "outer" flux form operator based on the Piecewise Parabolic Method (PPM; see also Lin and Rood 1996) for computing the fluxes in the N - S direction.

## Arguments
- `jlmt::Integer` - IN
- `cry::Matrix{AbstractFloat}` - IN
- `dcy::Matrix{AbstractFloat}` - IN
- `qqu::Matrix{AbstractFloat}` - IN
- `qqv::Matrix{AbstractFloat}` - OUT
- `j1p::Integer` - IN
- `j2p::Integer` - IN
- `i1_gl::Integer` - IN
- `i2_gl::Integer` - IN
- `ju1_gl::Integer` - IN
- `j2_gl::Integer` - IN
- `ilong::Integer` - IN
- `ilo::Integer` - IN
- `ihi::Integer` - IN
- `julo::Integer` - IN
- `jhi::Integer` - IN
- `i1::Integer` - IN
- `i2::Integer` - IN
- `ju1::Integer` - IN
- `j2::Integer` - IN

## Author
Original code from Shian - Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S - J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos - chem for complete history.
"""
function fyppm!(
  jlmt::Integer,
  cry::Matrix{AbstractFloat},
  dcy::Matrix{AbstractFloat},
  qqu::Matrix{AbstractFloat},
  qqv::Matrix{AbstractFloat},
  j1p::Integer,
  j2p::Integer,
  i1_gl::Integer,
  i2_gl::Integer,
  ju1_gl::Integer,
  j2_gl::Integer,
  ilong::Integer,
  ilo::Integer,
  ihi::Integer,
  julo::Integer,
  jhi::Integer,
  i1::Integer,
  i2::Integer,
  ju1::Integer,
  j2::Integer
)::Nothing
  a61 = zeros(AbstractFloat, ilong * ((jhi - 1) - (julo + 1) + 1))
  al1 = zeros(AbstractFloat, ilong * ((jhi - 1) - (julo + 1) + 1))
  ar1 = zeros(AbstractFloat, ilong * ((jhi - 1) - (julo + 1) + 1))
  dcy1 = zeros(AbstractFloat, ilong * ((jhi - 1) - (julo + 1) + 1))
  qqu1 = zeros(AbstractFloat, ilong * ((jhi - 1) - (julo + 1) + 1))
  a6 = zeros(AbstractFloat, ilo:ihi, julo:jhi)
  al = zeros(AbstractFloat, ilo:ihi, julo:jhi)
  ar = zeros(AbstractFloat, ilo:ihi, julo:jhi)

  # NOTE: The code was writtein with I1:I2 as the first dimension of AL, AR, A6, AL1, A61, AR1. However, the limits should really should be ILO:IHI. In practice, however, for a global grid (and OpenMP parallelization) ILO=I1 and IHI=I2. Nevertheless, we will change the limits to ILO:IHI to be consistent and to avoid future problems. (bmy, 12/5/08)

  r13 = 1.0 / 3.0
  r23 = 2.0 / 3.0

  for ij = (julo + 1):jhi, il = ilo:ihi
    al[il, ij] = 0.5 * (qqu[il, ij - 1] + qqu[il, ij]) + (dcy[il, ij - 1] - dcy[il, ij]) * r13
    ar[il, ij - 1] = al[il, ij]
  end

  # TODO:
  # call do_fyppm_pole_i2d2(al, ar, i1_gl, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2)

  for ij = (julo + 1):(jhi - 1), il = ilo:ihi
    a6[il, ij] = 3.0 * (qqu[il, ij] + qqu[il, ij] - (al[il, ij] + ar[il, ij]))
  end

  if jlmt <= 2
    lenx = 0

    for ij = (julo + 1):(jhi - 1), il = ilo:ihi
      lenx = lenx + 1

      a61[lenx] = a6[il, ij]
      al1[lenx] = al[il, ij]
      ar1[lenx] = ar[il, ij]
      dcy1[lenx] = dcy[il, ij]
      qqu1[lenx] = qqu[il, ij]
    end

    # TODO:
    # call Lmtppm(lenx, jlmt, a61, al1, ar1, dcy1, qqu1)

    lenx = 0

    for ij = (julo + 1):(jhi - 1), il = ilo:ihi
      lenx = lenx + 1

      a6[il, ij] = a61[lenx]
      al[il, ij] = al1[lenx]
      ar[il, ij] = ar1[lenx]
    end
  end

  for ij = j1p:j2p + 1
    ijm1 = ij - 1

    for il = ilo:ihi
      if cry[il, ij] > 0.0
        qqv[il, ij] = ar[il, ijm1] + 0.5 * cry[il, ij] * (al[il, ijm1] - ar[il, ijm1] + (a6[il, ijm1] * (1.0 - (r23 * cry[il, ij]))))
      else
        qqv[il, ij] = al[il, ij] - 0.5 * cry[il, ij] * (ar[il, ij] - al[il, ij] + (a6[il, ij] * (1.0 + (r23 * cry[il, ij]))))
      end
    end
  end
end

"""
Subroutine do_fyppm_pole_i2d2 sets "al" & "ar" at the Poles.

## Arguments
- `al::Matrix{AbstractFloat}` - INOUT
- `ar::Matrix{AbstractFloat}` - INOUT
- `i1_gl::Integer` - IN
- `i2_gl::Integer` - IN
- `ju1_gl::Integer` - IN
- `j2_gl::Integer` - IN
- `ilo::Integer` - IN
- `ihi::Integer` - IN
- `julo::Integer` - IN
- `jhi::Integer` - IN
- `i1::Integer` - IN
- `i2::Integer` - IN
- `ju1::Integer` - IN
- `j2::Integer` - IN

## Author
Original code from Shian - Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S - J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos - chem for complete history
"""
function do_fyppm_pole_i2d2!(
  al::Matrix{AbstractFloat},
  ar::Matrix{AbstractFloat},
  i1_gl::Integer,
  i2_gl::Integer,
  ju1_gl::Integer,
  j2_gl::Integer,
  ilo::Integer,
  ihi::Integer,
  julo::Integer,
  jhi::Integer,
  i1::Integer,
  i2::Integer,
  ju1::Integer,
  j2::Integer
)::Nothing
  i2d2 = i2_gl / 2

  for il = i1:i2d2
    al[il, ju1] = al[il + i2d2, ju1 + 1]
    al[il + i2d2, ju1] = al[il, ju1 + 1]
    ar[il, j2] = ar[il + i2d2, j2 - 1]
    ar[il + i2d2, j2] = ar[il, j2 - 1]
  end
end

"""
Subroutine do_ytp_pole_sum sets "dq1" at the Poles.

## Arguments
- `geofac_pc::AbstractFloat` - IN
- `dq1::Matrix{AbstractFloat}` - INOUT
- `qqv::Matrix{AbstractFloat}` - IN
- `fy::Matrix{AbstractFloat}` - INOUT
- `i1_gl::Integer` - IN
- `i2_gl::Integer` - IN
- `ju1_gl::Integer` - IN
- `j2_gl::Integer` - IN
- `j1p::Integer` - IN
- `j2p::Integer` - IN
- `ilo::Integer` - IN
- `ihi::Integer` - IN
- `julo::Integer` - IN
- `jhi::Integer` - IN
- `i1::Integer` - IN
- `i2::Integer` - IN
- `ju1::Integer` - IN
- `j2::Integer` - IN

## Author
Original code from Shian - Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S - J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos - chem for complete history.
"""
function do_ytp_pole_sum!(
  geofac_pc::AbstractFloat,
  dq1::Matrix{AbstractFloat},
  qqv::Matrix{AbstractFloat},
  fy::Matrix{AbstractFloat},
  i1_gl::Integer,
  i2_gl::Integer,
  ju1_gl::Integer,
  j2_gl::Integer,
  j1p::Integer,
  j2p::Integer,
  ilo::Integer,
  ihi::Integer,
  julo::Integer,
  jhi::Integer,
  i1::Integer,
  i2::Integer,
  ju1::Integer,
  j2::Integer
)::Nothing
  ri2 = i2_gl

  dqik = zeros(AbstractFloat, 2)

  # ... Integrate N - S flux around polar cap lat circle for each level

  sumsp = 0.0
  sumnp = 0.0
  for il = i1:i2
    sumsp = sumsp + qqv[il, j1p]
    sumnp = sumnp + qqv[il, j2p + 1]
  end

  # ... wrap in E - W direction
  if i1 == i1_gl
    dqik[1] = dq1[i1, ju1]
    dqik[2] = dq1[i1, j2]
  end

  # ... normalize and set inside polar cap

  dq_sp = dqik[1] - (sumsp / ri2 * geofac_pc)
  dq_np = dqik[2] + (sumnp / ri2 * geofac_pc)

  for il = i1:i2
    dq1[il, ju1] = dq_sp
    dq1[il, j2] = dq_np
    # ... save polar flux
    fy[il, ju1] =  - (sumsp / ri2 * geofac_pc)
    fy[il, j2 + 1] = (sumnp / ri2 * geofac_pc)
  end

  # TODO: `/=` works diferently in Julia from Fortran.
  if j1p /= ju1_gl + 1
    for il = i1:i2
      dq1[il, ju1 + 1] = dq_sp
      dq1[il, j2 - 1] = dq_np
      # ... save polar flux
      fy[il, ju1 + 1] =  - (sumsp / ri2 * geofac_pc)
      fy[il, j2] = sumnp / ri2 * geofac_pc
    end
  end
end

"""
Subroutine fzppm is the 1D "outer" flux form operator based on the Piecewise Parabolic Method (PPM; see also Lin and Rood 1996) for computing the fluxes in the vertical direction.

fzppm was modified by S. - J. Lin, 12/14/98, to allow the use of the KORD=7 (klmt=4) option. KORD=7 enforces the 2nd monotonicity constraint of Huynh (1996). Note that in Huynh's original scheme, two constraints are necessary for the preservation of monotonicity. To use Huynh's algorithm, it was modified as follows. The original PPM is still used to obtain the first guesses for the cell edges, and as such Huynh's 1st constraint is no longer needed. Huynh's median function is also replaced by a simpler yet functionally equivalent in - line algorithm.

## Arguments
- `klmt::Integer` - IN
- `delp1::Array{AbstractFloat, 3}` - IN
- `wz::Array{AbstractFloat, 3}` - IN
- `dq1::Array{AbstractFloat, 3}` - INOUT
- `qq1::Array{AbstractFloat, 3}` - IN
- `fz::Array{AbstractFloat, 3}` - OUT
- `j1p::Integer` - IN
- `ju1_gl::Integer` - IN
- `j2_gl::Integer` - IN
- `ilo::Integer` - IN
- `ihi::Integer` - IN
- `julo::Integer` - IN
- `jhi::Integer` - IN
- `ilong::Integer` - IN
- `ivert::Integer` - IN
- `i1::Integer` - IN
- `i2::Integer` - IN
- `ju1::Integer` - IN
- `j2::Integer` - IN
- `k1::Integer` - IN
- `k2::Integer` - IN

## Author
Original code from Shian - Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S - J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos - chem for complete history
"""
function fzppm!(
  klmt::Integer,
  delp1::Array{AbstractFloat,3},
  wz::Array{AbstractFloat,3},
  dq1::Array{AbstractFloat,3},
  qq1::Array{AbstractFloat,3},
  fz::Array{AbstractFloat,3},
  j1p::Integer,
  ju1_gl::Integer,
  j2_gl::Integer,
  ilo::Integer,
  ihi::Integer,
  julo::Integer,
  jhi::Integer,
  ilong::Integer,
  ivert::Integer,
  i1::Integer,
  i2::Integer,
  ju1::Integer,
  j2::Integer,
  k1::Integer,
  k2::Integer
)::Nothing
  a61 = zeros(AbstractFloat, ilong * (ivert - 4))
  al1 = zeros(AbstractFloat, ilong * (ivert - 4))
  ar1 = zeros(AbstractFloat, ilong * (ivert - 4))
  dca1 = zeros(AbstractFloat, ilong * (ivert - 4))
  qq1a1 = zeros(AbstractFloat, ilong * (ivert - 4))

  a6 = zeros(AbstractFloat, i1:i2, k1:k2)
  al = zeros(AbstractFloat, i1:i2, k1:k2)
  ar = zeros(AbstractFloat, i1:i2, k1:k2)
  dca = zeros(AbstractFloat, i1:i2, k1:k2)
  dlp1a = zeros(AbstractFloat, i1:i2, k1:k2)
  qq1a = zeros(AbstractFloat, i1:i2, k1:k2)
  wza = zeros(AbstractFloat, i1:i2, k1:k2)

  dc = zeros(AbstractFloat, i1:i2, ju1:j2, k1:k2)

  # Work array
  dp = zeros(AbstractFloat, i1:i2, ju1:j2, k1:k2)

  # !.sds... diagnostic vertical flux for species - set top to 0.0
  fz[:, :, :] .= zeros(AbstractFloat, ilo:ihi, julo:jhi, k1:k2)

  k1p1 = k1 + 1
  k1p2 = k1 + 2

  k2m1 = k2 - 1
  k2m2 = k2 - 2

  r13 = 1.0 / 3.0
  r23 = 2.0 / 3.0

  # Compute dc for PPM.

  for ik = k1:k2m1
    dpi[:, :, ik] .= qq1[i1:i2, ju1:j2, ik + 1] - qq1[i1:i2, ju1:j2, ik]
  end

  for ik = k1p1:k2m1, ij = ju1:j2, il = i1:i2
    c0 = delp1[il, ij, ik] / (delp1[il, ij, ik - 1] + delp1[il, ij, ik] + delp1[il, ij, ik + 1])
    c1 = (delp1[il, ij, ik - 1] + (0.5 * delp1[il, ij, ik])) / (delp1[il, ij, ik + 1] + delp1[il, ij, ik])
    c2 = (delp1[il, ij, ik + 1] + (0.5 * delp1[il, ij, ik])) / (delp1[il, ij, ik - 1] + delp1[il, ij, ik])

    tmp = c0 * ((c1 * dpi[il, ij, ik]) + (c2 * dpi[il, ij, ik - 1]))

    qmin = qq1[il, ij, ik] - min(qq1[il, ij, ik - 1], qq1[il, ij, ik], qq1[il, ij, ik + 1])
    qmax = max(qq1[il, ij, ik - 1], qq1[il, ij, ik], qq1[il, ij, ik + 1]) - qq1[il, ij, ik]

    dc[il, ij, ik] = min(abs(tmp), qmax, qmin) * sign(tmp)
  end

  # c?
  # Loop over latitudes (to save memory).

  @label ijloop
  for ij = ju1:j2
    # TODO: `/=` has a different meaning in Julia than it does in Fortran 
    if (ij == ju1_gl + 1 || ij == j2_gl - 1) && (j1p /= ju1_gl + 1)
      @goto ijloop
    end

    for ik = k1:k2, il = i1:i2
      dca[il, ik] = dc[il, ij, ik] # the monotone slope
      wza[il, ik] = wz[il, ij, ik]
      dlp1a[il, ik] = delp1[il, ij, ik]
      qq1a[il, ik] = qq1[il, ij, ik]
    end

    # Compute first guesses at cell interfaces. First guesses are required to be continuous. Three - cell parabolic subgrid distribution at model top; two - cell parabolic with zero gradient subgrid distribution at the surface.

    # First guess top edge value.

    for il = i1:i2
      # Three - cell PPM; compute a, b, & c of q = aP^2 + bP + c using cell averages and dlp1a.

      fac1 = dpi[il, ij, k1p1] - dpi[il, ij, k1] * (dlp1a[il, k1p1] + dlp1a[il, k1p2]) / (dlp1a[il, k1] + dlp1a[il, k1p1])

      fac2 = (dlp1a[il, k1p1] + dlp1a[il, k1p2]) * (dlp1a[il, k1] + dlp1a[il, k1p1] + dlp1a[il, k1p2])

      aa = 3.0 * fac1 / fac2

      bb = 2.0 * dpi[il, ij, k1] / (dlp1a[il, k1] + dlp1a[il, k1p1]) - r23 * aa * (2.0 * dlp1a[il, k1] + dlp1a[il, k1p1])

      al[il, k1] = qq1a[il, k1] - dlp1a[il, k1] * (r13 * aa * dlp1a[il, k1] + 0.5 * bb)

      al[il, k1p1] = dlp1a[il, k1] * (aa * dlp1a[il, k1] + bb) + al[il, k1]

      # Check if change sign.

      if qq1a[il, k1] * al[il, k1] <= 0.0
        al[il, k1] = 0.0
        dca[il, k1] = 0.0
      else
        dca[il, k1] = qq1a[il, k1] - al[il, k1]
      end
    end

    # Bottom.
    for il = i1:i2
      # 2 - cell PPM with zero gradient right at the surface.

      fac1 = dpi[il, ij, k2m1] * (dlp1a[il, k2] * dlp1a[il, k2]) / ((dlp1a[il, k2] + dlp1a[il, k2m1]) * (2.0 * dlp1a[il, k2] + dlp1a[il, k2m1]))

      ar[il, k2] = qq1a[il, k2] + fac1
      al[il, k2] = qq1a[il, k2] - (fac1 + fac1)

      if qq1a[il, k2] * ar[il, k2] <= 0.0
        ar[il, k2] = 0.0
      end

      dca[il, k2] = ar[il, k2] - qq1a[il, k2]
    end

    # 4th order interpolation in the interior.

    for ik = k1p2:k2m1, il = i1:i2
      c1 = (dpi[il, ij, ik - 1] * dlp1a[il, ik - 1]) / (dlp1a[il, ik - 1] + dlp1a[il, ik])
      c2 = 2.0 / (dlp1a[il, ik - 2] + dlp1a[il, ik - 1] + dlp1a[il, ik] + dlp1a[il, ik + 1])

      a1 = (dlp1a[il, ik - 2] + dlp1a[il, ik - 1]) / (2.0 * dlp1a[il, ik - 1] + dlp1a[il, ik])
      a2 = (dlp1a[il, ik] + dlp1a[il, ik + 1]) / (2.0 * dlp1a[il, ik] + dlp1a[il, ik - 1])

      al[il, ik] = qq1a[il, ik - 1] + c1 + c2 * (dlp1a[il, ik] * (c1 * (a1 - a2) + a2 * dca[il, ik - 1]) - dlp1a[il, ik - 1] * a1 * dca[il, ik])
    end

    for ik = k1:k2m1, il = i1:i2
      ar[il, ik] = al[il, ik + 1]
    end

    # Top & Bottom 2 layers always monotonic.

    lenx = i2 - i1 + 1

    for ik = k1:k1p1
      for il = i1:i2
        a6[il, ik] = 3.0 * (qq1a[il, ik] + qq1a[il, ik] - (al[il, ik] + ar[il, ik]))
      end

      # TODO:
      # call Lmtppm(lenx, 0, a6[i1, ik], al[i1, ik], ar[i1, ik], dca[i1, ik], qq1a[i1, ik])
    end

    for ik = k2m1:k2
      for il = i1:i2
        a6[il, ik] = 3.0 * (qq1a[il, ik] + qq1a[il, ik] - (al[il, ik] + ar[il, ik]))
      end

      # TODO:
      # call Lmtppm(lenx, 0, a6[i1, ik], al[i1, ik], ar[i1, ik], dca[i1, ik], qq1a[i1, ik])
    end

    # Interior depending on klmt.

    if klmt == 4
      # KORD=7, Huynh's 2nd constraint.

      for ik = k1p1:k2m1, il = i1:i2
        dca[il, ik] = dpi[il, ij, ik] - dpi[il, ij, ik - 1]
      end

      for ik = k1p2:k2m2, il = i1:i2
        # Right edges.

        qmp = qq1a[il, ik] + (2.0 * dpi[il, ij, ik - 1])
        lac = qq1a[il, ik] + (1.5 * dca[il, ik - 1]) + (0.5 * dpi[il, ij, ik - 1])
        qmin = min(qq1a[il, ik], qmp, lac)
        qmax = max(qq1a[il, ik], qmp, lac)

        ar[il, ik] = min(max(ar[il, ik], qmin), qmax)

        # Left edges.

        qmp = qq1a[il, ik] - (2.0 * dpi[il, ij, ik])
        lac = qq1a[il, ik] + (1.5 * dca[il, ik + 1]) - (0.5 * dpi[il, ij, ik])
        qmin = min(qq1a[il, ik], qmp, lac)
        qmax = max(qq1a[il, ik], qmp, lac)

        al[il, ik] = min(max(al[il, ik], qmin), qmax)

        # Recompute a6.

        a6[il, ik] = 3.0 * (qq1a[il, ik] + qq1a[il, ik] - (ar[il, ik] + al[il, ik]))
      end
    elseif klmt <= 2
      lenx = 0

      for ik = k1p2:k2m2, il = i1:i2
        lenx = lenx + 1

        al1[lenx] = al[il, ik]
        ar1[lenx] = ar[il, ik]
        dca1[lenx] = dca[il, ik]
        qq1a1[lenx] = qq1a[il, ik]

        a61[lenx] = 3.0 * (qq1a1[lenx] + qq1a1[lenx] - (al1[lenx] + ar1[lenx]))
      end

      # TODO:
      # call Lmtppm(lenx, klmt, a61, al1, ar1, dca1, qq1a1)

      lenx = 0
      for ik = k1p2:k2m2, il = i1:i2
        lenx = lenx + 1

        a6[il, ik] = a61[lenx]
        al[il, ik] = al1[lenx]
        ar[il, ik] = ar1[lenx]
        dca[il, ik] = dca1[lenx]
        qq1a[il, ik] = qq1a1[lenx]
      end
    end

    for ik = k1:k2m1, il = i1:i2
      if wza[il, ik] > 0.0
        cm = wza[il, ik] / dlp1a[il, ik]

        dca[il, ik + 1] = ar[il, ik] + 0.5 * cm * (al[il, ik] - ar[il, ik] + a6[il, ik] * (1.0 - r23 * cm))
      else
        cp = wza[il, ik] / dlp1a[il, ik + 1]

        dca[il, ik + 1] = al[il, ik + 1] + 0.5 * cp * (al[il, ik + 1] - ar[il, ik + 1] - a6[il, ik + 1] * [1.0 + r23 * cp])

      end
    end

    for ik = k1:k2m1, il = i1:i2
      dca[il, ik + 1] = wza[il, ik] * dca[il, ik + 1]
      # .sds.. save vertical flux for species as diagnostic
      fz[il, ij, ik + 1] = dca[il, ik + 1]
    end

    for il = i1:i2
      dq1[il, ij, k1] = dq1[il, ij, k1] - dca[il, k1p1]
      dq1[il, ij, k2] = dq1[il, ij, k2] + dca[il, k2]
    end

    for ik = k1p1:k2m1, il = i1:i2
      dq1[il, ij, ik] = dq1[il, ij, ik] + dca[il, ik] - dca[il, ik + 1]
    end
  end
end

"""
Subroutine average_press_poles averages pressure at the Poles when the Polar cap is enlarged. It makes the last two latitudes equal.

## Arguments
- `area_1d::Array{AbstractFloat}` - IN
- `press::Matrix{AbstractFloat}` - INOUT
- `i1::Integer` - IN
- `i2::Integer` - IN
- `ju1::Integer` - IN
- `j2::Integer` - IN
- `ilo::Integer` - IN
- `ihi::Integer` - IN
- `julo::Integer` - IN
- `jhi::Integer` - IN

## Author
Philip Cameron - Smith and John Tannahill, GMI project @ LLNL (2003).
Implemented into GEOS - Chem by Claire Carouge (ccarouge@seas.harvard.edu).

## Remarks
Subroutine from pjc_pfix. Call this one once everything is working fine.

## Revision History
05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S - J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos - chem for complete history.
"""
function average_press_poles!(
  area_1d::Array{AbstractFloat},
  press::Matrix{AbstractFloat},
  i1::Integer,
  i2::Integer,
  ju1::Integer,
  j2::Integer,
  ilo::Integer,
  ihi::Integer,
  julo::Integer,
  jhi::Integer
)::Nothing
  # Compute the sum of surface area

  # TODO: `dble` is not an excisting function, it's probably imported from somewhere else

  rel_area = zeros(AbstractFloat, ju1:j2)
  # TODO: looks like a map
  sum_area = sum[area_1d] * dble(i2)
  # calculate rel_area for each lat. (ccc, 11/20/08)
  for j = ju1:j2
    rel_area[j] = area_1d[j] / sum_area
  end

  # South Pole

  # Surface area of the S. Polar cap
  sum_area = sum(rel_area[ju1:(ju1 + 1)]) * dble(i2)

  # TODO: This looks like a reduce
  # Zero
  meanp = 0.0
  # Sum pressure * surface area over the S. Polar cap
  for j = ju1:(ju1 + 1), i = i1:i2
    meanp = meanp + (rel_area[j] * press[i, j])
  end

  # Normalize pressure in all grid boxes w/in the S. Polar cap
  press[:, ju1:(ju1 + 1)] .= meanp / sum_area

  # North Pole

  # Surface area of the N. Polar cap
  sum_area = sum(rel_area[(j2 - 1):j2]) * dble(i2)

  # TODO: This looks like a reduce
  # Zero
  meanp = 0.0
  # ! Sum pressure * surface area over the N. Polar cap
  for j = (j2 - 1):j2, i = i1:i2
    meanp = meanp + (rel_area[j] * press[i, j])
  end

  # ! Normalize pressure in all grid boxes w/in the N. Polar cap
  press[:, (j2 - 1):j2] .= meanp / sum_area
end

end # module tpcore_fvdas_mod
