c
      subroutine wrapup(pot_id)
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'charge.h'
      include 'elecpot.h'
      include 'energy.h'
      include 'compat.h'
      include 'version.h'
      include 'pseudowave.h'
      include 'corecorr.h'
c
c     ecuts now set in compat_params...
c
      double precision zero, tpfive, one
      parameter (zero=0.d0,tpfive=2.5d0, one=1.d0)
c
      integer i, icore, j, jcut, lp, noi, npotd, npotu, ifull
      integer nops(norbmx), position, iunit
      character ray(6)*10, title*70, pot_id*40, id*1
      character id_original*1
c
      double precision zval, zratio, zion, ac, bc, cdcp, tanb, rbold,
     &                 rbnew, pi, ecut, vp2z, fcut, zot, vpsdm,
     &                 vps, rmind, vpsum, rminu, zelu, zeld, zelt,
     &                 viodj, viouj, cc
      double precision rcut(10), v(nrmax)
      double precision orb_charge(nrmax)
c
      double precision fourier_area(5), qc(5)
      double precision fourier_eps
      parameter (fourier_eps = 1.0d-2)

      integer n_shells_down(5), n_shells_up(5)
      integer n_channels, lun

      double precision wfn(nrmax)

      double precision cutoff_function, force_underflow
      external cutoff_function, force_underflow

      logical new_cc_scheme, defined, down_channel
c
      external logder, defined, down_channel
c
      pi = 4*atan(one)
c
c     Do not use relativity for what follows
c
      id_original = id
      if (polarized) then
         id = 's'
      else
         id = ' '
      endif
c
c  Reset the n quantum numbers to include all valence orbitals.
c  Compute the ratio between the valence charge present and the
c  valence charge of a neutral atom.
c  Transfer pseudo valence charge to charge array
c
      zval = zero
      zratio = zero
c
      do 5 i = 1, norb
         nops(i) = 0
    5 continue
c
!     nops is the "effective principal quantum number", which
!     sets the number of "nodes"     

      do lp = 1, lmax
         n_shells_down(lp) = 0
         n_shells_up(lp) = 0
      enddo
      do 10 i = ncp, norb
         lp = lo(i) + 1
         if (down(i)) then
            n_shells_down(lp) = n_shells_down(lp) + 1
            nops(i) = lo(i) + n_shells_down(lp)
         else
            n_shells_up(lp) = n_shells_up(lp) + 1
            nops(i) = lo(i) + n_shells_up(lp)
         endif
!         nops(i) = lo(i) + 1
         zval = zval + zo(i)
   10 continue
      zion = zval + znuc - zel
      if (zval .ne. zero) zratio = zion/zval
c

!     New method able to deal with non-pseudized states in the valence.

!     We could use vio{u,d} and vi{u,d} coming out of the
!     generation routine to test the whole configuration
!     and generate the total valence charge density 
!     (including any non-pseudized states) by calling:

!         call dsolv2(0,1,id,ncp,norb,0,nops)


!     On entry now: vio{u,d}: ionic ps (de-screened with the total AE charge)
!                   vi{u,d} : V_HXC(total AE charge)
!
!     The net V_ext would be then the screened pseudo that correctly
!     generates the pseudo-wavefunction for each channel, as constructed
!     by the appropriate pseudization scheme.
!     When used to generate any upper levels, it will give "pseudo-wfs"
!     with nodes, e.g., one node for Ba 6s if we have 5s as semicore and
!     pseudized. This is precisely what we need to generate the full
!     pseudo-valence charge for later de-screening.

!     Notes
!       For pseudized levels, this approach does not give exactly the same
!       charge density as the synthetic method used in each "ps-generator".
!
!       There is indeed a difference between the charge density generated
!       by the pseudizer and the output of the call to dsolv2: the integration
!       is done non-relativistically, whereas the original charge (outside
!       rc) might be relativistic.
!       On the other hand, the inversion of the Schrodinger equation in 
!       the pseudizer is non-relativistic, so if we use 'r' in 'id' we might 
!       impact the r<r_c shape of the charge density.

!       It is then be better to get the "pseudo" versions of non-pseudized
!       orbitals only, by managing the implicit loop in dsolv2 explicitly
!       here and adding the charge to the work arrays vo{u,d} which 
!       contain the 'pseudo-valence' charge constructed in the 'ps_generator'
!       routine by explicit squaring of the pseudo-wavefunction of
!       the pseudized levels.

         do 880 i = ncp, norb
            if (.not. pseudized(i)) then
               call dsolv2_single(0,id,i,nops,orb_charge)
               if (down(i)) then
                  do 580 j = 1, nr
                     vod(j) = vod(j) + orb_charge(j)
 580              continue
               else
                  do 590 j = 1, nr
                     vou(j) = vou(j) + orb_charge(j)
 590              continue
               endif
            endif
 880     continue

!     Now re-load the charge arrays with the total
!     pseudo-valence charge

         do 20 i = 1, nr
            cdd(i) = vod(i)
            cdu(i) = vou(i)
 20      continue

c
c=====================================================================
c  If a core correction is indicated construct pseudo core charge
c  if cfac < 0 or the valence charge is zero the full core is used
c
      ifcore = job - 1
      if (ifcore .ne. 0) then
         ac = zero
         bc = zero
         cc = zero
         icore = 1
         if (cfac .le. zero .or. zratio .eq. zero) then
            write(6,9000) r(icore), ac, bc, cc
            write(6,'(a)') '(Full core used)'
            call coreq
         else
            if (rcfac .le. zero) then
               do 30 i = nr, 2, -1
                  if (cdc(i) .gt. cfac*zratio*
     &                (cdd(i)+cdu(i))) go to 50
   30          continue
            else
               do 40 i = nr, 2, -1
                  if (r(i) .le. rcfac) go to 50
   40          continue
            end if
   50       continue
            icore = i
C
C--------------------------------------------------------------------
C Choice of methods to compute the core correction:
C
C 1. Traditional 'Froyen-Louie-Cohen' with cdc(r) = Arsin(Br)
C    and value and first-derivative matching. It is the default
c    for LDA calculations. See compat_params.f
C    and input.f for info on how to force its use from the input file.
C
C 2. New by Jose Luis Martins' group, using a Kerker-like exp( ) function
C    and matching also the second derivative. This is the default with
C    the 'mons' compatibility mode for GGA calculations.
C

            new_cc_scheme = .false.
            if (is_gga) then
               write(6,'(a)') 'Note: GGA calculation ==> New CC scheme'
               new_cc_scheme = .true.
            endif

            if ( (use_old_cc) .and. (is_gga)) then

               write(6,'(/,2a,/a)')
     $              'WARNING: Using old-style core corrections',
     $            ' despite this being a GGA calculation.',
     $              'I hope you know what you are doing...'
               new_cc_scheme = .false.
            endif

            if ( (use_new_cc) .and. (.not. is_gga)) then

               write(6,'(/,2a,/a)')
     $              'WARNING: Using new core corrections',
     $            ' despite this being an LDA calculation.',
     $         'Results will not be compatible with older versions.'
               new_cc_scheme = .true.
            endif

            rc_core = r(icore)

            if (.not. new_cc_scheme) then

               n_of_continuous_derivs = 1

c           Fit to  cdc(r) = ac*r * sin(bc*r) inside r(icore)
c
c           Find derivative at core radius. Use a five-point formula
c           instead of the old two-point method:
c              cdcp = (cdc(icore+1)-cdc(icore))/(r(icore+1)-r(icore))
c           Since the r abscissae are not equally spaced, we perform
c           the derivative using the chain rule:
c
c           r(i) = a [ exp(b*(i-1)) - 1 ]
c           r(x) = a [ exp(b*x) - 1 ]
c
c           f'(r) = f'(x) / r'(x)
c
c           r'(x) = a*b*exp(b*x) = (r(x)+a)*b
c           r'(i) = a*b*exp[b*(i-1)] = rab(i)
c
c           To compute f'(x), use eq. 25.3.6 
c           (p. 883) in Abramowitz-Stegun, with h=1
c
c           f'(0) = 1/12 f(-2) - 2/3 f(-1) + 2/3 f(1) - 1/12 f(2)
c
cold            cdcp = (cdc(icore+1)-cdc(icore))/(r(icore+1)-r(icore))

            cdcp =   cdc(icore-2) - 8*cdc(icore-1) +
     $             8*cdc(icore+1) -   cdc(icore+2)
            cdcp = cdcp / 12.d0 / rab(icore)
c
c           Now fit ac and bc using the function and derivative
c           information. 
c
c           g(r) = Arsin(Br) ==>  g / (rg'-g) = tanBr/(Br)
c
c           Use an iterative method to find Br and from that ac and bc.
c           Start near Br=2.5 so that we are in an invertible region 
c           of the tanx/x function.
c 
c
            tanb = cdc(icore)/(r(icore)*cdcp-cdc(icore))
            rbold = tpfive
            do 70 i = 1, 50
               rbnew = pi + atan(tanb*rbold)
               if (abs(rbnew-rbold) .lt. .00001D0) then
                  bc = rbnew/r(icore)
                  ac = cdc(icore)/(r(icore)*sin(rbnew))
                  do 60 j = 1, icore
                     cdc(j) = ac*r(j)*sin(bc*r(j))
   60             continue
                  write(6,9000) r(icore), ac, bc, cc
c
                  call coreq
                  go to 80
c
               else
                  rbold = rbnew
               end if
   70       continue
            write(6,9010)
            call ext(830)

            else

c                 Use subroutine provided by JLM to fit
c                 cdc(r) = r^2*exp(ac+bc*r^2+cc*r^4) inside r(icore)
c
                  n_of_continuous_derivs = 2

                  CALL PCC_EXP(NR,ICORE,AC,BC,CC,R,CDC)
                  write(6,9000) r(icore), ac, bc, cc

c
           endif

           call coreq

         end if
      end if
C---------------------------------------------------------------------
 9000 format(//' Core correction used',/' Pseudo core inside r =',f6.3,
     &      /' ac =',f6.3,' bc =',f6.3,' cc =',f6.3,/)
 9010 format(//' Error in pseudo - nonconvergence in finding ',
     &      /'pseudo-core values')
c
c  End the pseudo core charge.
c======================================================================

   80 continue
!
c
c  Compute the potential due to pseudo valence charge.
!  (Total configuration in the new method)
c
c  njtj  ***  NOTE  ***
c  Spin-polarized potentials should be unscreened with
c  spin-polarized valence charge.  This was not
c  done in pseudo and pseudok in earlier versions
c  of this program.
c  njtj  ***  NOTE  ***
c
!     Computes vo{u,d}, which is V_HXC for the current charge
!     density in cd{u,d} (plus cdc if core-corrections are used)

!     On entry now, cd{u,d} was traditionally synthetically generated from
!     the pseudo-wavefunctions. In the new method, it is the complete
!     valence pseudo-charge, including any upper states.

      call Velect(0,1,id,zval)
c
c  Construct the ionic pseudopotential and find the cutoff.
c  ecut should be adjusted to give a reasonable ionic cutoff
c  radius, but should not alter the pseudopotential, ie.,
c  the ionic cutoff radius should not be inside the pseudopotential
c  cutoff radius.
c
c  Note that the cutting off of the pseudopotentials (making
c  them approach -2*Zion/r faster) is not strictly necessary.
c  It might even be argued that it should be left to "client"
c  programs to decide what to do.

      write(6,9020)
 9020 format(/)
      ecut = ecuts


      do 150 lp = 1, lmax
         if (indd(lp) .ne. 0) then
            i = indd(lp)
            do 90 j = 2, nr
               ! This is the screened pseudopotential
               ! vid is still the AE V_HXC
               v(j) = viod(lp,j)/r(j) + vid(j)
               ! De-screen with the potential from the pseudo-charge
               ! (including all the valence and the pseudo-core)
               viod(lp,j) = viod(lp,j) + (vid(j)-vod(j))*r(j)
               vp2z = viod(lp,j) + 2*zion
               if (abs(vp2z) .gt. ecut) jcut = j
   90       continue

c           Default cutoff function: f(r)=exp(-5*(r-r_cut)). It damps
c           down the residual of rV+2*Zion.
c           Should be made smoother... Vps ends up with a kink at rcut.
c           Maybe use one of the Vanderbilt generalized gaussians.

            rcut(i-ncore) = r(jcut)
            if (rcut(i-ncore) .lt. rc(lp)) then
               write(6,'(a,2f8.4)') 'Vps rcut point moved out to rc: ',
     $              rcut(i-ncore), rc(lp)
               rcut(i-ncore) = rc(lp)
            endif
            do 100 j = jcut, nr
cag               fcut = exp(-5*(r(j)-r(jcut)))
               fcut = cutoff_function(r(j)-r(jcut))
               viod(lp,j) = -2*zion + fcut*(viod(lp,j)+2*zion)
  100       continue

         endif

         if (indu(lp) .ne. 0) then
            i = indu(lp)
            do 120 j = 2, nr
               v(j) = viou(lp,j)/r(j) + viu(j)
               viou(lp,j) = viou(lp,j) + (viu(j)-vou(j))*r(j)
               vp2z = viou(lp,j) + 2*zion
               if (abs(vp2z) .gt. ecut) jcut = j
  120       continue

            rcut(i-ncore) = r(jcut)
            if (rcut(i-ncore) .lt. rc(lp)) then
               write(6,'(a,2f8.4)') 'Vps rcut point moved out to rc: ',
     $              rcut(i-ncore), rc(lp)
               rcut(i-ncore) = rc(lp)
            endif
            do 130 j = jcut, nr
cag               fcut = exp(-5*(r(j)-r(jcut)))
               fcut = cutoff_function(r(j)-r(jcut))
               viou(lp,j) = -2*zion + fcut*(viou(lp,j)+2*zion)
  130       continue

         end if
c
  150 continue
c
c   Convert spin-polarized potentials back to nonspin-polarized
c   by occupation weight(zo).  Assumes core polarization is
c   zero, ie. polarization is only a valence effect.
c
      if (polarized) then
         do 180 lp=1,lmax
            if ((indd(lp) .eq. 0) .or. (indu(lp) .eq. 0)) then
               goto 180
            endif
            i = indd(lp)
!!         do 180 i = ncp, norb, 2
!!            lp = lo(i) + 1
            zot = zo(i) + zo(i+1)
            if (zot .ne. zero) then
               do 160 j = 2, nr
                  viod(lp,j) = (viod(lp,j)*zo(i)+viou(lp,j)*zo(i+1))/zot
                  viou(lp,j) = viod(lp,j)
  160          continue
            else
               do 170 j = 2, nr
                  viod(lp,j) = viod(lp,j)/2 + viou(lp,j)/2
                  viou(lp,j) = viod(lp,j)
  170          continue
            end if
  180    continue
      end if
c
!     Reset V_HXC(AE) for further calculations to 
!     V_HXC(pseudo-valence + pseudo-core)
!
      do 190 i = 2, nr
         vid(i) = vod(i)
         viu(i) = vou(i)
  190 continue
c
c   Test the pseudopotential self consistency.  Spin-polarized
c   is tested as spin-polarized(since up/down potentials are
c   now the same)
c
!     This uses implicitly (via common blocks):
!         -- The un-screened ionic ps
!         -- The (pseudo-only) vi{d,u}
!
      call dsolv2(0,1,id,ncp,norb,0,nops)
c
c  Printout the pseudo eigenvalues after cutoff.
c
      write(6,9030)
 9030 format(//' test of eigenvalues',//' rcut :')
 9035 format((2x,i1,a1,f7.2,2x,f8.5))
      do i = ncp, norb
         lp = lo(i) + 1
         if ((indd(lp) .eq. i) .or. (indu(lp) .eq. i)) then
            write(6,9035) no(i), il(lo(i)+1),rcut(i-ncore), ev(i)
         endif
      enddo
c
c  Printout the data for potentials.
c
      write(6,9050)
 9050 format(///' l    vps(0)    vpsmin      at r',/)
c
      do 220 i = 1, lmax
         if (indd(i)+indu(i) .eq. 0) go to 220
         if (indd(i) .ne. 0) then
            vpsdm = zero
            do 200 j = 2, nr
               if (r(j) .lt. .00001D0) go to 200
               vps = viod(i,j)/r(j)
               if (vps .lt. vpsdm) then
                  vpsdm = vps
                  rmind = r(j)
               end if
  200       continue
            write(6,9060) il(i), viod(i,2)/r(2), vpsdm, rmind
         end if
         if (indu(i) .ne. 0) then
            vpsum = zero
            do 210 j = 2, nr
               if (r(j) .lt. .00001D0) go to 210
               vps = viou(i,j)/r(j)
               if (vps .lt. vpsum) then
                  vpsum = vps
                  rminu = r(j)
               end if
  210       continue
            write(6,9060) il(i), viou(i,2)/r(2), vpsum, rminu
         end if
 9060    format(1x,a1,3f10.3)
  220 continue
c
c   Print out the energies from etotal. (Valence only...)
c
      call etotal(ncp,norb)
c
c
c  Find the jobname and date.
c
      ray(1) = atom_id
      call cal_date(ray(2))
c  
      read(pot_id,'(4a10)') (ray(i),i=3,6)
c
c  Encode the title array.
c
      title = ' '
      position = 1
      do 240 i = 1, lmax
         if (indd(i) .eq. 0 .and. indu(i) .eq. 0) go to 240
         zelu = zero
         zeld = zero
         if (indd(i) .ne. 0) then
            noi = no(indd(i))
            zeld = zo(indd(i))
         end if
         if (indu(i) .ne. 0) then
            noi = no(indu(i))
            zelu = zo(indu(i))
         end if
         zelt = zeld + zelu
         if ( .not. polarized) then
            write(title(position:),9070) noi, il(i), zelt, ispp, rc(i)
 9070       format(i1,a1,f5.2,a1,' r=',f5.2,'/')
            position = position + 17
         else
            write(title(position:),9090)
     $                             noi, il(i), zeld, zelu, ispp, rc(i)
 9090       format(i1,a1,f4.2,',',f4.2,a1,f4.2,'/')
            position = position + 17
         end if
  240 continue
c
c  Construct relativistic sum and difference potentials.
c
      if (relativistic) then
c
c   ***  The s potential is from now on considered as "down", even
c        though s=0.5 in the relativistic case.
c
         if (indu(1) .eq. 0) go to 260
         indd(1) = indu(1)
         indu(1) = 0
         do 250 j = 2, nr
            viod(1,j) = viou(1,j)
            viou(1,j) = zero
  250    continue
  260    continue
         do 280 i = 2, lmax
            if (indd(i) .eq. 0 .or. indu(i) .eq. 0) go to 280
            do 270 j = 2, nr
               viodj = viod(i,j)
               viouj = viou(i,j)
               viod(i,j) = ((i-1)*viodj+i*viouj)/(2*i-1)
               viou(i,j) = 2*(viouj-viodj)/(2*i-1)
  270       continue
  280    continue
      end if

c  Determine the number of  potentials.  Coded them as
c  two digits, where the first digit is the number
c  of down or sum potentials and the second the number of
c  up or difference potentials.
c
      npotd = 0
      npotu = 0
      do 290 i = 1, lmax
         if (indd(i) .ne. 0) npotd = npotd + 1
         if (indu(i) .ne. 0) npotu = npotu + 1
  290 continue
!
!     Plotting and Fourier analysis should properly be
!     done here, for the "down" versions only.
!     The screened potential is no longer plotted

      n_channels = 0
      do 500 lp = 1, lmax
         if (indd(lp) .ne. 0) then
            n_channels = n_channels + 1
            do 510 j = 2, nr
               v(j) = viod(lp,j)/r(j)
  510       continue
            call potran(lp,v,r,nr,zion,fourier_area(n_channels),
     $                  fourier_eps,qc(n_channels))
            call potrv(v,r,nr-120,lp-1,zion)
         endif
 500  enddo

!     Compute the logarithmic derivative as a function of energy 
!     An energy range appropriate for the full valence complex
!     is selected, and only the "down" potentials are used.
!
      if (logder_radius .gt. 0.d0) call logder('PS')
c
c     Write out the Fourier area for each pseudo channel
c
      call get_unit(lun)
      open(unit=lun,file="FOURIER_AREA",form="formatted",
     $     status="unknown")
      rewind(lun)
      write(lun,"(i4)") n_channels
      write(lun,"(5f10.5)") (fourier_area(j),j=1,n_channels)
      close(lun)
c
c     Write out the Fourier threshold for each channel
c
      call get_unit(lun)
      open(unit=lun,file="FOURIER_QMAX",form="formatted",
     $     status="unknown")
      rewind(lun)
      write(lun,"(i4)") n_channels
      write(lun,"(5f14.5)") (qc(j),j=1,n_channels)
      close(lun)
c
c     Generate and save all the valence pseudowavefunctions
c     as generated by the "down" (average, major) pseudos.
c
c     Store in arrays in pseudowave.h
c
      nshells_stored = 0
      do 885 i = ncp, norb
         if (down_channel(i)) then
            call dsolv2_save_wf(0,id,i,nops,wfn)
            nshells_stored = nshells_stored + 1
            n_pswf(nshells_stored) = no(i)
            l_pswf(nshells_stored) = lo(i)
            e_pswf(nshells_stored) = ev(i)
            do j = 2, nr
               pswf(nshells_stored,j) = wfn(j)
            enddo
         endif
 885  continue
c
c  Write the heading to the current pseudo.dat
c  file (unit=1).
c
      ifull = 0
      if (cfac .le. zero .or. zratio .eq. zero) ifull = 1
      if (ifcore .eq. 1) then
         if (ifull .eq. 0) then
            nicore = 'pcec'
         else
            nicore = 'fcec'
         end if
      else if (ifcore .eq. 2) then
         if (ifull .eq. 0) then
            nicore = 'pche'
         else
            nicore = 'fche'
         end if
      else
         nicore = 'nc  '
      end if
      if (polarized) then
         irel = 'isp'
      else if (relativistic) then
         irel = 'rel'
      else
         irel = 'nrl'
      end if

      open(unit=1,file='VPSOUT',status='unknown',form='unformatted')
      rewind 1
      open(unit=2,file='VPSFMT',status='unknown',form='formatted')
      rewind 2
      open(unit=3,file='PSWFFMT',status='unknown',form='formatted')
      rewind 3

      write(1) nameat, icorr, irel, nicore, (ray(i),i=1,6), title,
     &         npotd, npotu, nr - 1, a, b, zion
      write(1) (r(i),i=2,nr)
c
      do 700 iunit=2,3
         write(iunit,8005) nameat, icorr, irel, nicore
         write(iunit,8010) (ray(j),j=1,6), title
         write(iunit,8015) npotd, npotu, nr-1, a, b, zion
         write(iunit,8040) 'Radial grid follows' 
         write(iunit,8030) (r(j),j=2,nr)
 700  continue
c
 8000 format(1x,i2)
 8005 format(1x,a2,1x,a2,1x,a3,1x,a4)
 8010 format(1x,6a10,/,1x,a70)
 8015 format(1x,2i3,i5,3g20.12)
 8030 format(4(g20.12))
 8040 format(1x,a)
c
c  Write the potentials to file (unit=1, and 2).
c
      do 300 i = 1, lmax
         if (indd(i) .eq. 0) go to 300
         write(1) i - 1, (viod(i,j),j=2,nr)
         write(2,8040) 'Down Pseudopotential follows (l on next line)'
         write(2,8000) i-1
         write(2,8030) (force_underflow(viod(i,j)),j=2,nr)
  300 continue
      do 310 i = 1, lmax
         if (indu(i) .eq. 0) go to 310
         write(1) i - 1, (viou(i,j),j=2,nr)
         write(2,8040) 'Up Pseudopotential follows (l on next line)'
         write(2,8000) i-1
         write(2,8030) (force_underflow(viou(i,j)),j=2,nr)
 310  continue
c
c  Write the charge densities      
c  Note that this charge density is the "pseudo" one.
c
      write(2,8040) 'Core charge follows'

      if (ifcore .ne. 1) then
         write(1) (zero,i=2,nr)
         write(2,8030) (zero,i=2,nr)
      else
         write(1) (cdc(i),i=2,nr)
         write(2,8030) (force_underflow(cdc(i)),i=2,nr)
      end if
      write(1) (zratio*(cdd(i)+cdu(i)),i=2,nr)
      write(2,8040) 'Valence charge follows'
      write(2,8030) (force_underflow(zratio*(cdd(i)+cdu(i))),i=2,nr)

      close(1)
      close(2)
c
c  Write the pseudo-wavefunctions (only in formatted form, in
c  auxiliary file) 'u_n,l (r) = 1/r R_n,l (r)'
c
 8045 format(1x,i2,i2)
      do 600 i = 1, lmax
         if (indd(i) .eq. 0) go to 600
         write(3,8040)
     $      'Down Pswavefunction (R/r) follows (l,n on next line)'
         write(3,8045) i-1, no(indd(i))
         write(3,8030) (force_underflow(pswfnrd(i,j)),j=2,nr)
  600 continue
      do 620 i = 1, lmax
         if (indu(i) .eq. 0) go to 620
         write(3,8040)
     $      'Up Pswavefunction (R/r) follows (l,n on next line)'
         write(3,8045) i-1, no(indu(i))
         write(3,8030) (force_underflow(pswfnru(i,j)),j=2,nr)
  620 continue
c
      close(3)
c
      open(unit=3,file='PSCHARGE',form='formatted',status='unknown')
c
c     NOTE: We no longer put "zratio" here!!!
c     (We still do in the ps file for compatibility with PW and
c      SIESTA) 
c     (Only affects plots for ionic configurations)
c     BEAR THIS IN MIND IF YOU ARE USING THE HEURISTIC CORE CORRECTION
c     CRITERION: If you specify a given "pc_weight" in the input file,
c     do not be surprised if the plot does not show rcore
c     in the place you expect it to be.
c
      if (ifcore .ne. 1) then
         do 400 j = 2, nr
            write(3,9900) r(j), cdu(j), cdd(j), zero
 400     continue
      else
         do 410 j = 2, nr
            write(3,9900) r(j), cdu(j), cdd(j), cdc(j)
 410     continue
      endif
c
 9900 format(1x,f15.10,3x,3f15.8)
c
      close(unit=3)
c
c     Write the pseudopotential in PSML format
c
      call write_psml( ray, npotd, npotu, zion )
c
      return
c
      CONTAINS

      function extrapol(f,r) result(val)
      ! extrapolate quadratically to zero
      double precision, intent(in) :: f(2:3), r(2:3)
      double precision             :: val

      double precision :: r2

      r2 = r(2)/(r(3)-r(2))
      val = f(2) - (f(3)-f(2))*r2

      end function extrapol

      end subroutine wrapup

      logical function down_channel(i)
      integer, intent(in) :: i

      include 'orbital.h'
      include 'param.h'

      if (relativistic) then
         down_channel = ((lo(i) .eq. 0) .or. down(i))
      else
         down_channel = down(i)
      endif

      end function down_channel
          
      double precision function cutoff_function(r)
      implicit none
c
c     Generalized cutoff function
c
      double precision r

      logical defined
      double precision vander
      external defined
      external vander
c
c     Standard cutoff
c
      if (defined('DO_NOT_CUT_TAILS')) then
         cutoff_function = 1.d0
      else if (defined('SMOOTH_TAIL_CUTOFF')) then
         cutoff_function = vander(1.d0,3.d0*r)
      else
         cutoff_function = exp(-5.d0*r)
      endif

      end

!      The famous "Vanderbilt generalized cutoff"
!
       function vander(a,x) result(f)
       integer, parameter :: dp = kind(1.d0)
       real(dp), intent(in) :: a    ! Generalized gaussian shape
       real(dp), intent(in) :: x    
       real(dp)             :: f

       real(dp), parameter :: log10_e = 0.4343
       real(dp), parameter :: exp_range = (range(1.0_dp)-1)/log10_e

!!     real(dp), parameter :: exp_range = 40.0_dp
       real(dp)   :: gexp

       gexp = sinh(a*x)/sinh(a)
       gexp = gexp*gexp

       if (gexp .lt. exp_range) then
          f=exp(-gexp)  
       else
          f = 0.0_dp
       endif

       end function vander
