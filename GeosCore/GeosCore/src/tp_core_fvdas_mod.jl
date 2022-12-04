module tpcore_fvdas_mod

"""
Subroutine Do_Cross_Terms_Pole_I2d2 sets "va" at the Poles.

## Arguments
- `cry::Integer` - IN
- `va::Integer` - OUT
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
function do_cross_terms_pole_i2d2!(
	cry::Matrix{AbstractFloat},
	va::Matrix{AbstractFloat},
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
	
	if j1p == ju1_gl + 1
		# Polar Cap NOT Enlarged: Get cross terms for N-S horizontal advection.
		
		if ju1 == ju1_gl
			for il = i1:i2d2
				va[il, ju1] = 0.5 * (cry[il, ju1 + 1] - cry[il + i2d2, ju1 + 1])
				va[il + i2d2, ju1] = -va[il, ju1]
			end
		end

		if j2 == j2_gl
			for il = i1:i2d2
				va[il, j2] = 0.5 * (cry[il, j2] - cry[il + i2d2, j2 - 1])
				va[il + i2d2, j2] = -va[il, j2]
			end
		end
	end
end

"""
Subroutine Xadv_Dao2 is the advective form E-W operator for computing the adx (E-W) cross term.

## Arguments
- `iad::Integer` - IN
- `jn::Integer` - IN
- `js::Integer` - IN
- `adx::Matrix{AbstractFloat}` - OUT
- `qqv::Matrix{AbstractFloat}` - IN
- `ua::Matrix{AbstractFloat}` - IN
- `ilo::Integer` - IN
- `ihi::Integer` - IN
- `julo::Integer` - IN
- `jhi::Integer` - IN
- `ju1_gl::Integer` - IN
- `j2_gl::Integer` - IN
- `j1p::Integer` - IN
- `j2p::Integer` - IN
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
function xadv_dao2!(
	iad::Integer,
	jn::Integer,
	js::Integer,
	adx::Matrix{AbstractFloat},
	qqv::Matrix{AbstractFloat},
	ua::Matrix{AbstractFloat},
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	j1p::Integer,
	j2p::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer
)::Nothing
	qtmp = zeros(AbstractFloat, (i2 / 3):(i2 + i2 / 3), julo:jhi)

	# Zero output array
	adx = 0
	for ij = julo:jhi
		for il = 1:i2
			qtmp[il, ij] = qqv[il, ij]
		end

		for il = (-i2 / 3):0
			qtmp[il, ij] = qqv[i2 + il, ij]
		end

		for il = (i2 + 1):(i2 + i2 / 3)
			qtmp[il, ij] = qqv[il - i2, ij]
		end
	end
			
	if iad == 1
		# 1st order.
		
		for ij = j1p:j2p
			if ij <= js || ij >= jn
				# In Polar area.

				for il = i1:i2
					iu = ua[il, ij]
					riu = iu
					ru = ua[il, ij] - riu
					iu = il - iu

					if ua[il, ij] >= 0.0
						rdiff = qtmp[iu - 1, ij] - qtmp[iu, ij]
					else
						rdiff = qtmp[iu, ij] - qtmp[iu + 1, ij]
					end

					adx[il, ij] = (qtmp[iu, ij] - qtmp[il, ij]) + (ru * rdiff)
				end
			else # js < ij < jn
				# Eulerian upwind.
				
				for il = i1:i2
					ril = il
					iu = ril - ua[il, ij]
					
					adx[il, ij] = ua[il, ij] * (qtmp[iu, ij] - qtmp[iu + 1, ij])
				end
			end
		end
	elseif iad == 2
		for ij = j1p:j2p
			if ij <= js || ij >= jn
				# In Polar area.
				
				for il = i1:i2
					iu = round(ua[il, ij])
					riu = iu
					ru = riu - ua[il, ij]
					iu = il - iu

					a1 = 0.5 * (qtmp[iu + 1, ij] + qtmp[iu - 1, ij]) - qtmp[iu, ij]

					b1 = 0.5 * (qtmp[iu + 1, ij] - qtmp[iu - 1, ij])

					c1 = qtmp[iu, ij] - qtmp[il, ij]

					adx[il, ij] = (ru * ((a1 * ru) + b1)) + c1
				end
			else # js < ij < jn
				# Eulerian upwind.

				for il = i1:i2
					iu = round(ua[il, ij])
					riu = iu
					ru = riu - ua[il, ij]
					iu = il - iu

					a1 = 0.5 * (qtmp[iu + 1, ij] + qtmp[iu - 1, ij]) - qtmp[iu, ij]

					b1 = 0.5 * (qtmp[iu + 1, ij] - qtmp[iu - 1, ij])

					c1 = qtmp[iu, ij] - qtmp[il, ij]

					adx[il, ij] = (ru * ((a1 * ru) + b1)) + c1
				end
			end
		end
	end
	
	if ju1 == ju1_gl
		adx[i1:i2, ju1] .= 0.0

		if j1p != ju1_gl + 1
			adx[i1:i2, ju1 + 1] .= 0.0
		end
	end

	if j2 == j2_gl
		adx[i1:i2, j2] .= 0.0

		if j1p != ju1_gl + 1
			adx[i1:i2, j2 - 1] .= 0.0
		end
	end
end

"""
Subroutine Yadv_Dao2 is the advective form N-S operator for computing the ady (N-S) cross term.

## Arguments
- `iad::Integer` - IN
- `ady::Matrix{AbstractFloat}` - OUT
- `qqu::Matrix{AbstractFloat}` - IN
- `va::Matrix{AbstractFloat}` - IN
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
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function yadv_dao2!(
	iad::Integer,
	ady::Matrix{AbstractFloat},
	qqu::Matrix{AbstractFloat},
	va::Matrix{AbstractFloat},
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
	# We may need a small ghost zone depending on the polar cap used
	qquwk = zeros(AbstractFloat, ilo:ihi, (julo - 2):(jhi + 2))

	# Zero output array
	ady = 0

	# Make work array
	for ij = julo:jhi
		qquwk[:, ij] .= qqu[:, ij]
	end

	# This routine creates a ghost zone in latitude in case of not enlarged polar cap (ccc, 11/20/08)
	# call do_yadv_pole_i2d2(qqu, qquwk, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
	
	if iad == 1
		# 1st order.
		for ij = (j1p - 1):(j2p + 1), il = i1:i2
			# !c?
			rij = ij
			jv = rij - va[il, ij]

			ady[il, ij] = va[il, ij] * (qquwk[il, jv] - qquwk[il, jv + 1])
		end
	elseif iad == 2
		for ij = (j1p - 1):(j2p + 1), il = i1:i2
			# c?
			jv  = round(va[il, ij])
			rjv = jv
			rv  = rjv - va[il, ij]
			jv  = ij - jv

			a1 = 0.5 * (qquwk[il, jv + 1] + qquwk[il, jv - 1]) - qquwk[il, jv]

			b1 = 0.5 * (qquwk[il, jv + 1] - qquwk[il, jv - 1])

			c1 = qquwk[il, jv] - qquwk[il, ij]

			ady[il, ij] = (rv * ((a1 * rv) + b1)) + c1
		end
	end
	
	# call do_yadv_pole_sum( ady, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
end

"""
Subroutine Do_Yadv_Pole_I2d2 sets "qquwk" at the Poles.

## Arguments
- `qqu::Integer` - IN
- `qquwk::Integer` - IN
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
# 05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function do_yadv_pole_i2d2!(
	qqu::Matrix{AbstractFloat},
	qquwk::Matrix{AbstractFloat},
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
	
	if j1p == ju1_gl + 1
		# Polar Cap NOT Enlarged.

		if ju1 == ju1_gl
			for il = i1:i2d2, inb = 1:2
				qquwk[il, ju1 - inb] = qqu[il + i2d2, ju1 + inb]
				qquwk[il + i2d2, ju1 - inb] = qqu[il, ju1 + inb]
			end
		end

		if j2 == j2_gl
			for il = i1:i2d2, inb = 1:2
				qquwk[il, j2 + inb] = qqu[il + i2d2, j2 - inb]
				qquwk[il + i2d2, j2 + inb] = qqu[il, j2 - inb]
			end
		end
	end
end

"""
Subroutine Do_Yadv_Pole_Sum sets the cross term due to N-S advection at the Poles.

## Arguments
- `ady::Matrix{AbstractFloat}` - INOUT
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
function do_yadv_pole_sum!(
	ady::Matrix{AbstractFloat},
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
	# Test if we are using extended polar caps (i.e. the S pole and next N latitude and N. Pole and next S latitude).  Do this outside the loops. (bmy, 12/11/08)
	is_ext_polar_cap = j1p == ju1_gl + 2

	# South Pole

	sumsp = 0.0
	sumnp = 0.0

	if is_ext_polar_cap
		# For a 2-latitude polar cap (S. Pole + next Northward latitude)
		for il = i1:i2
			sumsp = sumsp + ady[il, ju1 + 1]
			sumnp = sumnp + ady[il, j2 - 1]
		end
	else
		# For a 1-latitude polar cap (S. Pole only)
		for il = i1:i2
			sumsp = sumsp + ady[il, ju1]
			sumnp = sumnp + ady[il, j2]
		end
	end

	sumsp = sumsp / i2_gl
	sumnp = sumnp / i2_gl
				
	if is_ext_polar_cap 
		# For a 2-latitude polar cap (S. Pole + next Northward latitude)
		for il = i1:i2
			ady[il, ju1 + 1] = sumsp
			ady[il, ju1] = sumsp
			ady[il, j2 - 1] = sumnp
			ady[il, j2] = sumnp
		end
	else
		# For a 1-latitude polar cap (S. Pole only)
		for il = i1:i2
			ady[il, ju1] = sumsp
			ady[il, j2] = sumnp
		end
	end
end

"""
Subroutine Xtp does horizontal advection in the E-W direction.

## Arguments
- `ilmt::Integer` - IN
- `jn::Integer` - IN
- `js::Integer` - IN
- `pu::Matrix{AbstractFloat}` - IN
- `crx::Matrix{AbstractFloat}` - IN
- `dq1::Matrix{AbstractFloat}` - INOUT
- `qqv::Matrix{AbstractFloat}` - INOUT
- `xmass::Matrix{AbstractFloat}` - IN
- `fx::Matrix{AbstractFloat}` - OUT
- `j1p::Integer` - IN
- `j2p::Integer` - IN
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
- `iord::Integer` - IN

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function xtp!(
	ilmt::Integer,
	jn::Integer,
	js::Integer,
	pu::Matrix{AbstractFloat},
	crx::Matrix{AbstractFloat},
	dq1::Matrix{AbstractFloat},
	qqv::Matrix{AbstractFloat},
	xmass::Matrix{AbstractFloat},
	fx::Matrix{AbstractFloat},
	j1p::Integer,
	j2p::Integer,
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
	j2::Integer,
	iord::Integer
)::Nothing
	isav = zeros(Integer, i1:i2)
	dcx = zeros(AbstractFloat, (-i2 / 3):(i2 + i2 / 3), julo:jhi)
	qtmp = zeros(AbstractFloat, (-i2 / 3):(i2 + i2 / 3), julo:jhi)
	fx[:, :] .= 0.0

	imp = i2 + 1

	# NOTE: these loops do not parallelize well (bmy, 12/5/08)

	# Populate qtmp
	for il = i1:i2
		qtmp[il, :] .= qqv[il, :]
	end

	for il = (-i2 / 3):0
		qtmp[il, :] = qqv[i2 + il, :]
	end

	for il = (i2 + 1):(i2 + i2 / 3)
		qtmp[il, :] = qqv[il - i2, :]
	end

	if iord != 1
		qtmp[i1 - 1, :] .= qqv[i2, :]
		qtmp[i1 - 2, :] .= qqv[i2 - 1, :]
		qtmp[i2 + 1, :] .= qqv[i1, :]
		qtmp[i2 + 2, :] .= qqv[i1 + 1, :]
		
		# call Xmist(dcx, qtmp, j1p, j2p, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
	end
	
	jvan = max(1, j2_gl / 18)
	
	for ij = j1p:j2p
		if ij > js && ij < jn
			# Do horizontal Eulerian advection in the E-W direction.
			
			if iord == 1 || ij == j1p || ij == j2p
				for il = i1:i2
					ril = il
					iu = ril - crx[il, ij]

					fx[il, ij] = qtmp[iu, ij]
				end
			else
				if iord == 2 || ij <= j1p + jvan || ij >= j2p - jvan
					for il = i1:i2
						ril = il
						iu = ril - crx[il, ij]

						fx[il, ij] = qtmp[iu, ij] + (dcx[iu, ij] * (sign(crx[il, ij]) - crx[il, ij]))
					end
				else	
					# call Fxppm(ij, ilmt, crx, dcx, fx, qtmp, -i2/3, i2+i2/3, julo, jhi, i1, i2)
					# qtmp (inout) - can be updated
				end
			end
			
			for il = i1:i2
				fx[il, ij] = fx[il, ij] * xmass[il, ij]
			end
		else
			# Do horizontal Conservative (flux-form) Semi-Lagrangian advection in the E-W direction (van Leer at high latitudes).
			if iord == 1 || ij == j1p  ij == j2p
				for il = i1:i2
					ic = crx[il, ij]
					isav[il] = il - ic
					ril = il
					iu = ril - crx[il, ij]
					ric = ic
					rc = crx[il, ij] - ric

					fx[il, ij] = rc * qtmp[iu, ij]
				end
			else
				for il = i1:i2
					ic = crx[il, ij]
					isav[il] = il - ic
					ril = il
					iu = ril - crx[il, ij]
					ric = ic
					rc = crx[il, ij] - ric

					fx[il, ij] = rc * (qtmp[iu, ij] + (dcx[iu, ij] * (sign(rc) - rc)))
				end
			end

			for il = i1:i2
				if crx[il, ij] > 1.0
					for ix = isav[il]:(il - 1)
						fx[il, ij] = fx[il, ij] + qtmp[ix, ij]
					end
				elseif crx[il, ij] < -1.0
					for ix = il:(isav[il] - 1)
						fx[il, ij] = fx[il, ij] - qtmp[ix, ij]
					end
				end
			end

			for il = i1:i2
				fx[il, ij] = pu[il, ij] * fx[il, ij]
			end
		end
	end
			
	# NOTE: This loop does not parallelize well (bmy, 12/5/08)
	for ij = j1p:j2p
		for il = i1:(i2 - 1)
			dq1[il, ij] = dq1[il, ij] + (fx[il, ij] - fx[il + 1, ij])
		end
		dq1[i2, ij] = dq1[i2, ij] + (fx[i2, ij] - fx[i1, ij])
	end
end

"""
Subroutine Xmist computes the linear tracer slope in the E-W direction. It uses the Lin et. al. 1994 algorithm.

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function xmist!(
	dcx::Matrix{AbstractFloat},
	qqv::Matrix{AbstractFloat},
	j1p::Integer,
	j2p::Integer,
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
	r24 = 1.0 / 24.0

	for ij = (j1p + 1):(j2p - 1), il = i1:i2
		tmp = ((8.0 * (qqv[il + 1, ij] - qqv[il - 1, ij])) + qqv[il - 2, ij] - qqv[il + 2, ij]) * r24

		pmax = max(qqv[il - 1, ij], qqv[il, ij], qqv[il + 1, ij]) - qqv[il, ij]

		pmin = qqv[il, ij] - min(qqv[il - 1, ij], qqv[il, ij], qqv[il + 1, ij])

		dcx[il, ij] = sign(tmp) * min(abs(tmp), pmax, pmin)
	end

	# Populate ghost zones of dcx (ccc, 11/20/08)
	for ij = julo:jhi
		for il = (-i2 / 3):0
			dcx[il, ij] = dcx[i2 + il, ij]
		end

		for il = (i2 + 1):(i2 + i2 / 3)
			dcx[il, ij] = dcx[il - i2, ij]
		end
	end
end

"""
Subroutine Fxppm is the 1D "outer" flux form operator based on the Piecewise Parabolic Method (PPM; see also Lin and Rood 1996) for computing the fluxes in the E-W direction.

## Arguments
- `ij::Integer` - IN
- `ilmt::Integer` - IN
- `crx::Matrix{AbstractFloat}` - IN
- `dcx::Matrix{AbstractFloat}` - OUT
- `fx::Matrix{AbstractFloat}` - OUT
- `qqv::Matrix{AbstractFloat}` - INOUT
- `ilo::Integer` - IN
- `ihi::Integer` - IN
- `julo::Integer` - IN
- `jhi::Integer` - IN
- `i1::Integer` - IN
- `i2::Integer` - IN

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Remarks
This routine is called from w/in a OpenMP parallel loop fro

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function fxppm!(
	ij::Integer,
	ilmt::Integer,
	crx::Matrix{AbstractFloat},
	dcx::Matrix{AbstractFloat},
	fx::Matrix{AbstractFloat},
	qqv::Matrix{AbstractFloat},
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer
)::Nothing
	# Zero arrays (bmy, 12/5/08)
	a6 = zeros(AbstractFloat, ilo:ihi)
	al = zeros(AbstractFloat, ilo:ihi)
	ar = zeros(AbstractFloat, ilo:ihi)
	a61 = zeros(AbstractFloat, (ihi-1) - (ilo+1) + 1)
	al1 = zeros(AbstractFloat, (ihi-1) - (ilo+1) + 1)
	ar1 = zeros(AbstractFloat, (ihi-1) - (ilo+1) + 1)
	dcxi1 = zeros(AbstractFloat, (ihi-1) - (ilo+1) + 1)
	qqvi1 = zeros(AbstractFloat, (ihi-1) - (ilo+1) + 1)

	r13 = 1.0 / 3.0
	r23 = 2.0 / 3.0

	for il = (ilo + 1):ihi
		rval = 0.5 * (qqv[il - 1, ij] + qqv[il, ij]) + (dcx[il - 1, ij] - dcx[il, ij]) * r13
		al[il] = rval
		ar[il - 1] = rval
	end

	for il = (ilo + 1):(ihi - 1)
		a6[il] = 3.0 * (qqv[il, ij] + qqv[il, ij] - (al[il] + ar[il]))
	end
	
	if ilmt <= 2
		a61[:] .= 0.0
		al1[:] .= 0.0
		ar1[:] .= 0.0

		dcxi1[:] .= 0.0
		qqvi1[:] .= 0.0

		lenx = 0
		for il = (ilo + 1):(ihi - 1)
			lenx = lenx + 1

			a61[lenx] = a6[il]
			al1[lenx] = al[il]
			ar1[lenx] = ar[il]

			dcxi1[lenx] = dcx[il, ij]
			qqvi1[lenx] = qqv[il, ij]
		end
		
		# call lmtppm(lenx, ilmt, a61, al1, ar1, dcxi1, qqvi1)

		lenx = 0
		for il = (ilo + 1):(ihi - 1)
			lenx = lenx + 1

			a6[il] = a61[lenx]
			al[il] = al1[lenx]
			ar[il] = ar1[lenx]

			dcx[il, ij] = dcxi1[lenx]
			qqv[il, ij] = qqvi1[lenx]
		end

		# Populate ghost zones of qqv and dcx with new values (ccc, 11/20/08)
		for il = (-i2 / 3):0
			dcx[il, ij] = dcx[i2 + il, ij]
			qqv[il, ij] = qqv[i2 + il, ij]
		end

		for il = (i2 + 1):(i2 + i2 / 3)
			dcx[il, ij] = dcx[il - i2, ij]
			qqv[il, ij] = qqv[il - i2, ij]
		end
	end

	for il = (i1 + 1):i2
		if crx[il, ij] > 0.0
			ilm1 = il - 1
			fx[il, ij] = ar[ilm1] + 0.5 * crx[il, ij] * (al[ilm1] - ar[ilm1] + (a6[ilm1] * (1.0 - (r23 * crx[il, ij]))))
		else
			fx[il, ij] = al[il] - 0.5 * crx[il, ij] * (ar[il] - al[il] + (a6[il] * (1.0 + (r23 * crx[il, ij]))))
		end
	end

	# First box case (ccc, 11/20/08)
	if crx[i1, ij] > 0.0
			ilm1 = i2
			fx[i1, ij] = ar[ilm1] + 0.5 * crx[i1, ij] * (al[ilm1] - ar[ilm1] + (a6[ilm1] * (1.0 - (r23 * crx[i1, ij]))))
	else
		fx[i1, ij] = al[i1] - 0.5 * crx[i1, ij] * (ar[i1] - al[i1] + (a6[i1] * (1.0 + (r23 * crx[i1, ij]))))
	end
end

"""
Subroutine Lmtppm enforces the full monotonic, semi-monotonic, or the positive-definite constraint to the sub-grid parabolic distribution of the Piecewise Parabolic Method (PPM).

## Arguments
- `lenx::Integer` - IN
- `lmt::Integer` - IN
- `a6::Array{AbstractFloat}` - INOUT
- `al::Array{AbstractFloat}` - INOUT
- `ar::Array{AbstractFloat}` - INOUT
- `dc::Array{AbstractFloat}` - INOUT
- `qa::Array{AbstractFloat}` - INOUT

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function lmtppm!(
	lenx::Integer,
	lmt::Integer,
	a6::Array{AbstractFloat},
	al::Array{AbstractFloat},
	ar::Array{AbstractFloat},
	dc::Array{AbstractFloat},
	qa::Array{AbstractFloat}
)::Nothing
	r12 = 1.0 / 12.0
	
	if lmt == 0
		# Full constraint.
		for il = 1:lenx
			if dc[il] == 0.0
					a6[il] = 0.0
					al[il] = qa[il]
					ar[il] = qa[il]
			else
				da1 = ar[il] - al[il]
				da2 = da1 * da1
				a6da = a6[il] * da1

				if a6da < -da2
					a6[il] = 3.0 * (al[il] - qa[il])
					ar[il] = al[il] - a6[il]
				elseif a6da > da2
					a6[il] = 3.0 * (ar[il] - qa[il])
					al[il] = ar[il] - a6[il]
				end
			end
		end
	elseif lmt == 1
		# Semi-monotonic constraint.
		for il = 1:lenx
			if abs(ar[il] - al[il]) < -a6[il]
				if qa[il] < ar[il] && qa[il] < al[il]
					a6[il] = 0.0
					al[il] = qa[il]
					ar[il] = qa[il]
				elseif ar[il] > al[il]
					a6[il] = 3.0 * (al[il] - qa[il])
					ar[il] = al[il] - a6[il]
				else
					a6[il] = 3.0 * (ar[il] - qa[il])
					al[il] = ar[il] - a6[il]
				end
			end
		end
	elseif lmt == 2
		for il = 1:lenx
			if abs(ar[il] - al[il]) < -a6[il]

				ftmp = ar[il] - al[il]

				fmin = qa[il] + 0.25 * (ftmp * ftmp) / a6[il] + a6[il] * r12

				if fmin < 0.0
					if qa[il] < ar[il] && qa[il] < al[il]
						a6[il] = 0.0
						al[il] = qa[il]
						ar[il] = qa[il]
					elseif ar[il] > al[il]
						a6[il] = 3.0 * (al[il] - qa[il])
						ar[il] = al[il] - a6[il]
					else
						a6[il] = 3.0 * (ar[il] - qa[il])
						al[il] = ar[il] - a6[il]
					end
				end
			end
		end
	end
end

"""
Subroutine Ytp does horizontal advection in the N-S direction.

## Arguments
- `jlmt::Integer` - IN
- `geofac_pc::AbstractFloat` - IN
- `geofac::Matrix{AbstractFloat}` - IN
- `cry::Matrix{AbstractFloat}` - IN
- `dq1::Matrix{AbstractFloat}` - INOUT
- `qqu::Matrix{AbstractFloat}` - IN
- `qqv::Matrix{AbstractFloat}` - INOUT
- `ymass::Matrix{AbstractFloat}` - IN
- `fy::Matrix{AbstractFloat}` - OUT
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
- `jord::Integer` - IN

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function ytp!(
	jlmt::Integer,
	geofac_pc::AbstractFloat,
	geofac::Matrix{AbstractFloat},
	cry::Matrix{AbstractFloat},
	dq1::Matrix{AbstractFloat},
	qqu::Matrix{AbstractFloat},
	qqv::Matrix{AbstractFloat},
	ymass::Matrix{AbstractFloat},
	fy::Matrix{AbstractFloat},
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
	j2::Integer,
	jord::Integer
)::Nothing
	dcy = zeros(AbstractFloat, ilo:ihi, julo:jhi)
	fy[:, :] .= 0.0

	rj1p = j1p
	
	if jord == 1
		for ij = j1p:(j2p + 1), il = i1:i2
			# c?
			jv = rj1p - cry[il, ij]
			qqv[il, ij] = qqu[il, jv]
		end
	else
			# TODO:
			# call ymist(4, dcy, qqu, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
		if jord <= 0 || jord >= 3
			# TODO:
			# call fyppm(jlmt, cry, dcy, qqu, qqv, j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl, ilong, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
		else
			for ij = j1p:(j2p + 1), il = i1:i2
				# c?
				jv = rj1p - cry[il, ij]
				qqv[il, ij] = qqu[il, jv] + ((sign(cry[il, ij]) - cry[il, ij]) * dcy[il, jv])
			end
		end
	end

	for ij = j1p:(j2p + 1)
		qqv[i1:i2, ij] = qqv[i1:i2, ij] * ymass[i1:i2, ij]
	end

	# .sds.. save N-S species flux as diagnostic
	for ij = i1:i2
		fy[ij, j1p:(j2p + 1)] = qqv[ij, j1p:(j2p + 1)] * geofac[j1p:(j2p + 1)]
	end

	# ... meridional flux update
	for ij = j1p:j2p
		dq1[i1:i2, ij] = dq1[i1:i2, ij] + (qqv[i1:i2, ij] - qqv[i1:i2, ij + 1]) * geofac[ij]
	end
	
	# TODO:
	# call do_ytp_pole_sum(geofac_pc, dq1, qqv, fy, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, j2p, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
end

"""
Subroutine ymist computes the linear tracer slope in the N-S direction.  It uses the Lin et. al. 1994 algorithm.

## Arguments
- `id::Integer` - IN
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
function ymist!(
	id::Integer,
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
	# I suppose the values for these indexes are 0. It should work as the pole values are re-calculated in the pole functions. (ccc)
	qtmp = zeros(AbstractFloat, ilo:ihi, (julo - 2):(jhi + 2))

	r24  = 1.0 / 24.0

	# Populate qtmp
	qtmp = 0.0
	for ij = ju1:j2
			qtmp[:, ij] .= qqu[:, ij]
	end
	
	if id == 2
		for ij = (ju1 - 1):(j2 - 1), il = i1:i2
			tmp  = 0.25 * (qtmp[il, ij + 2] - qtmp[il, ij])

			pmax = max(qtmp[il, ij], qtmp[il, ij + 1], qtmp[il, ij + 2]) - qtmp[il, ij + 1]

			pmin = qtmp[il, ij + 1] - min(qtmp[il, ij], qtmp[il, ij + 1], qtmp[il, ij + 2])

			dcy[il, ij + 1] = sign(tmp) * min(abs(tmp), pmin, pmax)
		end
	else
		# TODO:
		# call do_ymist_pole1_i2d2(dcy, qtmp, i1_gl, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2)

		for ij = (ju1 - 2):(j2 - 2), il = i1:i2
			tmp = ((8.0 * (qtmp[il, ij + 3] - qtmp[il, ij + 1])) + qtmp[il, ij] - qtmp[il, ij + 4]) * r24

			pmax = max(qtmp[il, ij + 1], qtmp[il, ij + 2], qtmp[il, ij + 3]) - qtmp[il, ij + 2]

			pmin = qtmp[il, ij + 2] - min(qtmp[il, ij + 1], qtmp[il, ij + 2], qtmp[il, ij + 3])

			dcy[il, ij + 2] = sign(tmp) * min(abs(tmp), pmin, pmax)
		end
	end
	# TODO:
	# call do_ymist_pole2_i2d2(dcy, qtmp, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
end

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
		if j1p != ju1_gl + 1
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
		if j1p != ju1_gl + 1
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

  if j1p != ju1_gl + 1
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
    if (ij == ju1_gl + 1 || ij == j2_gl - 1) && (j1p != ju1_gl + 1)
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
