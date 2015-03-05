program build_ff

  !! Build force fields for DL_POLY
 
  implicit none

  integer, parameter :: Ndim=3, Nmax=100000, icount=1
  integer, parameter :: izero=0, ione=1, itwo=2, ithree=3
  integer, parameter :: iconfig=0, icell=3

  real(8), parameter :: zero=0.d0

  !-----------------------------------------------------------------------

  character(20) :: pdb, topology, output

  !------------------------------------------------------------------------

  character(4) :: wff_type, uname1, uname2
  character(10) :: ff_key
  character(10), allocatable :: vdwtype(:), atname(:), ff_type(:)
  character(10), allocatable :: uname(:), utype(:)
  character(10), allocatable :: atype(:), atype1(:), atype2(:)
  character(10), allocatable :: atype3(:), atype4(:)

  integer :: iatom, Natom, iatom1, iatom2, iatom3, iatom4
  integer :: ibond, iangle, idihedral, ivdw, icheck
  integer :: itype, Ntype, itmp, iuni, Nuni, itype1, itype2
  integer :: isum, idiff, Nfac, wifac, wn, wifacok
  integer :: Nbond, Nangle, Ndihe, Nvdw, Nnew, Ncheck
  integer :: Nbond_fit, Nangle_fit, Ndihe_fit, Nimp_fit
  integer, allocatable :: iat1(:), iat4(:), n(:), ifac(:)
  
  real(8) :: charge_tot, pi, dot,ANG_BELOW,ANG_ABOVE,ANG
  real(8) :: box_x, box_y, box_z
  real(8) :: wk, k1, k1_new, k2_new, weq, eq1, eq1_new, eq2_new
  real(8) :: walpha, alpha1_new, alpha2_new, wphi
  real(8) :: wcc1, wcc2, wcc3
  real(8) :: ee, ss, e1, e2, s1, s2
  real(8) :: r2, r, cutoff, cutoff_OZnO
  real(8) :: ssx, ssy, ssz, xss, yss, zss, det
  real(8) :: dx(3), box(3), angle(3), cell(9), rcell(9)
  real(8) :: vec1(3),vec2(3),rA,rB
  real(8), allocatable :: pos(:,:), sigma(:), eps(:), charge(:)
  real(8), allocatable :: mass(:), e14(:), vdw14(:)
  real(8), allocatable :: ucharge(:), umass(:)
  real(8), allocatable :: k(:), eq(:), alpha(:), phi(:)
  real(8), allocatable :: cc1(:), cc2(:), cc3(:)

  character(10), allocatable :: atomtype(:), resname(:)
  integer :: ierr
  integer, allocatable :: resnum(:)

  !-------------------------------------------------------------------------

  character(10) :: junk
  integer :: ijunk, ncount, i
  real(8) :: rjunk
  
  !-------------------------------------------------------------------------

  namelist/input/Nuni, Natom, Nbond, Nangle, Ndihe, Nbond_fit, &
                 Nangle_fit, Ndihe_fit, Nimp_fit, Nvdw, Nfac, wifacok
  
  !-------------------------------------------------------------------------

  !! Input files.
  open(unit=1, file='input', status="old", iostat=ierr)
  if (ierr>0) stop "*** Can't open input ***"
  read(1,nml=input)
  close(1)

  open(unit=2, file='force_field.dat', status="old", iostat=ierr)
  if (ierr>0) stop "*** Can't open force_field.dat ***"
  
  open(unit=3, file='14scaling.dat', status="old", iostat=ierr)
  if (ierr>0) stop "*** Can't open 14scaling.dat ***"

  open(unit=11, file='reference.pdb', status="old", iostat=ierr)   
  if (ierr>0) stop "*** Can't open reference.pdb ***" 

  open(unit=21, file='new_params.dat', status="old", iostat=ierr)
  if (ierr>0) stop "*** Can't open new_params.dat ***"

  !! Output files.
  open(unit=101, file='FIELD_TMP')
  open(unit=102, file='CONFIG_TMP')

  !-------------------------------------------------------------------------

  !! Allocation.
  allocate(e14(0:Nfac-1))
  allocate(vdw14(0:Nfac-1))

  allocate(iat1(Nmax))
  allocate(iat4(Nmax))
  allocate(ff_type(Nmax))
  allocate(k(Nmax))
  allocate(eq(Nmax))
  allocate(alpha(Nmax))
  allocate(cc1(Nmax))
  allocate(cc2(Nmax))
  allocate(cc3(Nmax))
  allocate(phi(Nmax))
  allocate(n(Nmax))
  allocate(ifac(Nmax))
  allocate(uname(Nuni))
  allocate(utype(Nuni))
  allocate(atname(Natom))
  allocate(atype(Nmax))  
  allocate(atype1(Nmax))
  allocate(atype2(Nmax))
  allocate(atype3(Nmax))
  allocate(atype4(Nmax))
  allocate(umass(Nuni))
  allocate(mass(Natom))
  allocate(ucharge(Nuni))
  allocate(charge(Natom))
  allocate(pos(Ndim,Natom))

  !-------------------------------------------------------------------------

  !! Get some values.
  do i = 0, Nfac-1
     read(3,*) itmp, e14(i), vdw14(i)
  end do

  cell(:) = 0.d0
  pi = acos(-1.d0)

  !!===========================================
  !! Read input from Gaussian pdb.
  !!===========================================

  read(11,*) junk, box(1), box(2), box(3), &
             angle(1), angle(2), angle(3)

  call calculate_cell(box,angle,cell)
   
  call invert(cell,rcell,det)

  do iatom = 1, Natom
     !read(11,*) junk, ijunk, atname(iatom), &
     !           junk, ijunk, pos(:,iatom)
     read(11,*) junk, ijunk, atname(iatom), &
                ijunk, pos(:,iatom)
  end do
  close(11)
  
  ! Write CONFIG.
  write(102,'(a14)') "MOF SIMULATION"
  write(102,'(3i10,f20.10)') iconfig, icell, natom, zero
  write(102,'(3f20.10)') cell(1:3)
  write(102,'(3f20.10)') cell(4:6)
  write(102,'(3f20.10)') cell(7:9)
  
  do iatom = 1, Natom
     write(102,'(a8,i10)') atname(iatom), iatom
     write(102,'(3f20.10)') pos(:,iatom)
  end do


  !!===========================================
  !! Loop over atoms.
  !!===========================================

  read(2,*) ff_key, Nuni  
  do iuni = 1, Nuni
     read(2,*) uname(iuni), utype(iuni), &
               ucharge(iuni), umass(iuni)
  end do

  charge_tot = 0.d0
  ncount = 0
  write(101,'(a6)') 'DLPOLY'
  write(101,'(a10)') 'UNITS kcal'
  write(101,'(a11)') 'MOLECULES 1'
  write(101,'(a3)') 'MOF'
  write(101,'(a8)') 'NUMMOL 1'
  write(101,'(a10,i8)') ff_key, Natom

  do iatom = 1, Natom
     do iuni = 1, Nuni
        if (atname(iatom).eq.uname(iuni)) then
           atype(iatom) = utype(iuni) 
           charge(iatom) = ucharge(iuni) 
           mass(iatom) = umass(iuni) 
           ncount = ncount + 1
           exit
        end if
     end do
     charge_tot = charge_tot + charge(iatom)
     write(101,'(a8,2f12.7,i8)') &
        atype(iatom), mass(iatom), charge(iatom), icount
  end do
  read(2,*)

  print*, 'number of atom = ', ncount
  print*, 'total charge = ', charge_tot


  !!===========================================
  !! Loop over bonds.
  !!===========================================

  ibond = 0
  read(2,*) ff_key, Ntype
  !print*, ff_key, Ntype
  do itype = 1, Ntype
     read(2,*) atype1(itype), atype2(itype), ff_type(itype), &
               k(itype), eq(itype), cc1(itype), cc2(itype)
  end do
  read(2,*)

  do itype = 1, Nbond_fit
     read(21,*) k(itype), eq(itype), cc1(itype), cc2(itype)
  end do

  write(101,'(a10,i10)') ff_key, Nbond

  !! Identify 1st atom.
  do iatom1 = 1, Natom

     !! Identify 2nd atom.
     do iatom2 = iatom1+1, Natom

        dx(:) = pos(:,iatom2) - pos(:,iatom1)
        ssx = rcell(1)*dx(1) + rcell(4)*dx(2) + rcell(7)*dx(3)
        ssy = rcell(2)*dx(1) + rcell(5)*dx(2) + rcell(8)*dx(3)
        ssz = rcell(3)*dx(1) + rcell(6)*dx(2) + rcell(9)*dx(3)
        xss = ssx - nint(ssx)
        yss = ssy - nint(ssy)
        zss = ssz - nint(ssz)
        dx(1) = cell(1)*xss + cell(4)*yss + cell(7)*zss
        dx(2) = cell(2)*xss + cell(5)*yss + cell(8)*zss
        dx(3) = cell(3)*xss + cell(6)*yss + cell(9)*zss

        r2 = sum(dx(:)**2)
        r = sqrt(r2)

        do itype = 1, Ntype

           if ((atype(iatom1).eq.atype1(itype) .and. &
                atype(iatom2).eq.atype2(itype)) .or. &
               (atype(iatom1).eq.atype2(itype) .and. &
                atype(iatom2).eq.atype1(itype))) then

               call cut(atype(iatom1), atype(iatom2), cutoff)

               if (r.lt.cutoff) then
                  wff_type = ff_type(itype)
                  wk = k(itype)
                  weq = eq(itype)
                  wcc1 = cc1(itype)
                  wcc2 = cc2(itype)
                  ibond = ibond + 1
                  write(101,'(a4,2i8,4f15.6)') wff_type, iatom1, iatom2, &
                                               wk, weq, wcc1, wcc2
               end if

           end if

        end do

     end do

  end do

  print*, 'bonds = ', ibond
  !write(101,*) 'bonds  ',ibond


  !!===========================================
  !! Loop over angles.
  !!===========================================

  iangle = 0
  read(2,*) ff_key, Ntype
  !print*, ff_key, Ntype
  
  do itype = 1, Ntype
     read(2,*) atype1(itype), atype2(itype), atype3(itype), &
               ff_type(itype), k(itype), eq(itype), &
               cc1(itype), cc2(itype)
  end do
  read(2,*)
!override normal setting and tell it to look for one additional angle
  do itype = 1, Nangle_fit+1
     read(21,*) k(itype), eq(itype), cc1(itype), cc2(itype)
  end do
  write(101,'(a10,i10)') ff_key, Nangle

  !! Identify central atom.
  do iatom2 = 1, Natom

     !! Identify 1st atom.
     do iatom1 = 1, Natom
        if (iatom1.eq.iatom2) cycle
        dx(:) = pos(:,iatom1) - pos(:,iatom2)
        ssx = rcell(1)*dx(1) + rcell(4)*dx(2) + rcell(7)*dx(3)
        ssy = rcell(2)*dx(1) + rcell(5)*dx(2) + rcell(8)*dx(3)
        ssz = rcell(3)*dx(1) + rcell(6)*dx(2) + rcell(9)*dx(3)
        xss = ssx - nint(ssx)
        yss = ssy - nint(ssy)
        zss = ssz - nint(ssz)
        dx(1) = cell(1)*xss + cell(4)*yss + cell(7)*zss
        dx(2) = cell(2)*xss + cell(5)*yss + cell(8)*zss
        dx(3) = cell(3)*xss + cell(6)*yss + cell(9)*zss
        r2 = sum(dx(:)**2)
        r = sqrt(r2)
        rA=r
        vec1(:)=dx(:)
        call cut(atype(iatom1), atype(iatom2), cutoff)
        if (r.lt.cutoff) then

           !! Identify 3rd atom.
           do iatom3 = iatom1+1, Natom
              if (iatom3.eq.iatom1 .or. iatom3.eq.iatom2) cycle
              dx(:) = pos(:,iatom3) - pos(:,iatom2)
              ssx = rcell(1)*dx(1) + rcell(4)*dx(2) + rcell(7)*dx(3)
              ssy = rcell(2)*dx(1) + rcell(5)*dx(2) + rcell(8)*dx(3)
              ssz = rcell(3)*dx(1) + rcell(6)*dx(2) + rcell(9)*dx(3)
              xss = ssx - nint(ssx)
              yss = ssy - nint(ssy)
              zss = ssz - nint(ssz)
              dx(1) = cell(1)*xss + cell(4)*yss + cell(7)*zss
              dx(2) = cell(2)*xss + cell(5)*yss + cell(8)*zss
              dx(3) = cell(3)*xss + cell(6)*yss + cell(9)*zss
              r2 = sum(dx(:)**2)
              r = sqrt(r2)
              rB=r
              vec2(:)=dx(:)
              call cut(atype(iatom3), atype(iatom2), cutoff)
              if (r.lt.cutoff) then
                 dot=vec1(1)*vec2(1)+vec1(2)*vec2(2)+vec1(3)*vec2(3)
                 !! Loop over types.
                 do itype = 1, Ntype

                    if ((atype(iatom1).eq.atype1(itype) .and.  &
                         atype(iatom2).eq.atype2(itype) .and.  &
                         atype(iatom3).eq.atype3(itype)) .or.  &
                        (atype(iatom1).eq.atype3(itype) .and.  &
                         atype(iatom2).eq.atype2(itype) .and.  &
                         atype(iatom3).eq.atype1(itype))) then
                        
                       wff_type = ff_type(itype)
                       wk = k(itype)
                       weq = eq(itype)
                       wcc1 = cc1(itype)
                       wcc2 = cc2(itype)
! Enter a routine specific for your system.  I am trying to eliminate a o1-Ni-o1 triad that has two angles: 91 and 160.
                       ANG=ACOS(dot/(rA*rB))*180.0/pi
!print *, iatom1, iatom2, iatom3, ACOS(dot/(rA*rB))*180.0/pi
!print *, atype(iatom1), atype(iatom2), atype(iatom3), ANG
                       ANG_BELOW=70.0
                       ANG_ABOVE=100.0
                       if ((atype(iatom1).eq.'o1').and.(atype(iatom2).eq.'Ni').and.(atype(iatom3).eq.'o1')) then
                           if((ANG.gt.ANG_BELOW).and.(ANG.lt.ANG_ABOVE)) then
                               write(101,'(a4,3i8,4f15.6)') wff_type, iatom1, iatom2, iatom3, wk, weq, wcc1, wcc2!,ANG
                              iangle = iangle + 1
                           else
!Add in the special angle
                               wk = k(Nangle_fit+1)
                               weq = eq(Nangle_fit+1)
                               wcc1 = cc1(Nangle_fit+1)
                               wcc2 = cc2(Nangle_fit+1)
                               write(101,'(a4,3i8,5f15.6)') wff_type, iatom1, iatom2, iatom3, wk, weq, wcc1, wcc2, ANG
                               iangle = iangle + 1
                           endif
                       else
                       write(101,'(a4,3i8,4f15.6)') wff_type, iatom1, iatom2, iatom3, wk, weq, wcc1, wcc2
                       iangle = iangle + 1
                      end if 
                    end if

                 end do

              end if
           end do

        end if
     end do

  end do

  print*, 'angles = ', iangle
  !write(101,*) 'angles  ', iangle


  !!===========================================
  !! Loop over dihedrals.
  !!===========================================

  idihedral = 0
  read(2,*) ff_key, Ntype
  !print*, ff_key, Ntype
  do itype = 1, Ntype
     read(2,*) atype1(itype), atype2(itype), atype3(itype), atype4(itype), &
               ff_type(itype), k(itype), phi(itype), n(itype), ifac(itype)
  end do
  read(2,*)

  do itype = 1, Ndihe_fit
     read(21,*) k(itype), phi(itype), n(itype), ifac(itype)
     ifac(itype) = wifacok
  end do
  write(101,'(a10,i10)') ff_key, Ndihe

  !! Identify 2nd atom.
  do iatom2 = 1, Natom

     !! Identify 3rd atom.
     do iatom3 = iatom2+1, Natom
        if (iatom3.eq.iatom2) cycle
        dx(:) = pos(:,iatom3) - pos(:,iatom2)
        ssx = rcell(1)*dx(1) + rcell(4)*dx(2) + rcell(7)*dx(3)
        ssy = rcell(2)*dx(1) + rcell(5)*dx(2) + rcell(8)*dx(3)
        ssz = rcell(3)*dx(1) + rcell(6)*dx(2) + rcell(9)*dx(3)
        xss = ssx - nint(ssx)
        yss = ssy - nint(ssy)
        zss = ssz - nint(ssz)
        dx(1) = cell(1)*xss + cell(4)*yss + cell(7)*zss
        dx(2) = cell(2)*xss + cell(5)*yss + cell(8)*zss
        dx(3) = cell(3)*xss + cell(6)*yss + cell(9)*zss
        r2 = sum(dx(:)**2)
        r = sqrt(r2)
        call cut(atype(iatom3), atype(iatom2), cutoff)
        if (r.lt.cutoff) then

           !! Identify 1st atom.
           do iatom1 = 1, Natom
              if (iatom1.eq.iatom2 .or. iatom1.eq.iatom3) cycle
              dx(:) = pos(:,iatom1) - pos(:,iatom2)
              ssx = rcell(1)*dx(1) + rcell(4)*dx(2) + rcell(7)*dx(3)
              ssy = rcell(2)*dx(1) + rcell(5)*dx(2) + rcell(8)*dx(3)
              ssz = rcell(3)*dx(1) + rcell(6)*dx(2) + rcell(9)*dx(3)
              xss = ssx - nint(ssx)
              yss = ssy - nint(ssy)
              zss = ssz - nint(ssz)
              dx(1) = cell(1)*xss + cell(4)*yss + cell(7)*zss
              dx(2) = cell(2)*xss + cell(5)*yss + cell(8)*zss
              dx(3) = cell(3)*xss + cell(6)*yss + cell(9)*zss
              r2 = sum(dx(:)**2)
              r = sqrt(r2)
              call cut(atype(iatom1), atype(iatom2), cutoff)
              if (r.lt.cutoff) then

                 !! Identify 4th atom.
                 do iatom4 = 1, Natom
                    if (iatom4.eq.iatom2 .or. iatom4.eq.iatom3) cycle
                    dx(:) = pos(:,iatom4)-pos(:,iatom3)
                    ssx = rcell(1)*dx(1) + rcell(4)*dx(2) + rcell(7)*dx(3)
                    ssy = rcell(2)*dx(1) + rcell(5)*dx(2) + rcell(8)*dx(3)
                    ssz = rcell(3)*dx(1) + rcell(6)*dx(2) + rcell(9)*dx(3)
                    xss = ssx - nint(ssx)
                    yss = ssy - nint(ssy)
                    zss = ssz - nint(ssz)
                    dx(1) = cell(1)*xss + cell(4)*yss + cell(7)*zss
                    dx(2) = cell(2)*xss + cell(5)*yss + cell(8)*zss
                    dx(3) = cell(3)*xss + cell(6)*yss + cell(9)*zss
                    r2 = sum(dx(:)**2)
                    r = sqrt(r2)
                    call cut(atype(iatom4), atype(iatom3), cutoff)
                    if (r.lt.cutoff) then

                       !! Loop over type.
                       do itype = 1, Ntype

                          if ((atype(iatom1).eq.atype1(itype) .and.  &
                               atype(iatom2).eq.atype2(itype) .and.  &
                               atype(iatom3).eq.atype3(itype) .and.  &
                               atype(iatom4).eq.atype4(itype)) .or.  &
                              (atype(iatom1).eq.atype4(itype) .and.  &
                               atype(iatom2).eq.atype3(itype) .and.  &
                               atype(iatom3).eq.atype2(itype) .and.  &
                               atype(iatom4).eq.atype1(itype))) then

                             idihedral = idihedral + 1
                             iat1(idihedral) = iatom1
                             iat4(idihedral) = iatom4
                             wff_type = ff_type(itype)
                             wk = k(itype)
                             wphi = phi(itype)
                             wn = n(itype)
                             wifac = ifac(itype)
                             do icheck = 1, idihedral-1
                                if ((iat1(icheck).eq.iatom1 .and. &
                                     iat4(icheck).eq.iatom4) .or. &
                                    (iat1(icheck).eq.iatom4 .and. &
                                     iat4(icheck).eq.iatom1)) then
                                   wifac = izero
                                end if
                             end do
                             write(101,'(a4,4i8,f15.6,f15.6,i8,2f15.6)') &
                                wff_type, iatom1, iatom2, iatom3, iatom4, &
                                wk, wphi, abs(wn), e14(wifac), vdw14(wifac)
                          end if

                       end do

                    end if
                 end do

              end if
           end do

        end if
     end do

  end do

  print*, 'dihedrals = ', idihedral
  !write(101,*) 'dihedral  ', idihedral

  !!===========================================
  !! Impropers.
  !!===========================================

  !! Read parameters.
  read(2,*) ff_key, Ntype
  !print*, ff_key, Ntype
  do itype = 1, Ntype
     read(2,*) atype1(itype), atype2(itype), atype3(itype), atype4(itype), &
               ff_type(itype), k(itype), phi(itype), n(itype), ifac(itype)
  end do
  read(2,*)

  do itype = 1, Nimp_fit
     read(21,*) k(itype), phi(itype), n(itype), ifac(itype)
  end do

  !! Loop over central atoms.
  do iatom3 = 1, Natom

     !! Identify 1st atom.
     do iatom1 = 1, Natom
        if (iatom1 .eq. iatom3) cycle
        dx(:) = pos(:,iatom1)-pos(:,iatom3)
        dx(:) = dx(:) - box(:)*nint(dx(:)/box(:))
        r2 = sum(dx(:)**2)
        r = sqrt(r2)
        call cut(atype(iatom1), atype(iatom3), cutoff)
        if (r .lt. cutoff) then

           !! Identify 2nd atom.
           do iatom2 = iatom1+1, Natom
              if (iatom2.eq.iatom3 .or. iatom2.eq.iatom1) cycle
              dx(:) = pos(:,iatom2)-pos(:,iatom3)
              dx(:) = dx(:) - box(:)*nint(dx(:)/box(:))
              r2 = sum(dx(:)**2)
              r = sqrt(r2)
              call cut(atype(iatom2), atype(iatom3), cutoff)
              if (r .lt. cutoff) then

                 !! Identify 3rd atom
                 do iatom4 = iatom2+1, Natom
                    if (iatom4.eq.iatom3 .or. iatom4.eq.iatom1 .or. iatom4.eq.iatom2) cycle
                    dx(:) = pos(:,iatom4)-pos(:,iatom3)
                    dx(:) = dx(:) - box(:)*nint(dx(:)/box(:))
                    r2 = sum(dx(:)**2)
                    r = sqrt(r2)
                    call cut(atype(iatom4), atype(iatom3), cutoff)
                    if (r .lt. cutoff) then
 
                       !! Loop over impropers.
                       do itype = 1, Ntype

                          if ((atype(iatom1).eq.atype1(itype) .and.  &
                               atype(iatom2).eq.atype2(itype) .and.  &
                               atype(iatom3).eq.atype3(itype) .and.  &
                               atype(iatom4).eq.atype4(itype)) .or.  &
                              (atype(iatom1).eq.atype2(itype) .and.  &
                               atype(iatom2).eq.atype4(itype) .and.  &
                               atype(iatom3).eq.atype3(itype) .and.  &
                               atype(iatom4).eq.atype1(itype)) .or.  &
                              (atype(iatom1).eq.atype4(itype) .and.  &
                               atype(iatom2).eq.atype1(itype) .and.  &
                               atype(iatom3).eq.atype3(itype) .and.  &
                               atype(iatom4).eq.atype2(itype))) then 

                             idihedral = idihedral + 1
                             wff_type = ff_type(itype)
                             wk = k(itype)
                             wphi = phi(itype)
                             wn = n(itype)
                             wifac = izero          
                             write(101,'(a4,4i8,f15.6,f15.6,i8,2f15.6)') &
                                wff_type, iatom3, iatom1, iatom2, iatom4, &
                                wk, wphi, abs(wn), e14(wifac), vdw14(wifac)
                          end if

                       end do

                    end if
                 end do

              end if
           end do

        end if
     end do

  end do

  print*, 'impropers = ', idihedral
  write(101,'(a6)') 'FINISH'

  !! VdW.
  read(2,*) ff_key, Ntype
  write(101,'(a10,i10)') ff_key, Nvdw
  allocate(vdwtype(Ntype))
  allocate(sigma(Ntype))
  allocate(eps(Ntype))
  do itype = 1, Ntype
     read(2,*) vdwtype(itype), wff_type, sigma(itype), eps(itype)
  end do
  ivdw = 0
  !do itype1 = 1, Ntype
  !   do itype2 = itype1, Ntype
  !      ss = sigma(itype1)+sigma(itype2)
  !      ee = sqrt(eps(itype1)*eps(itype2))
  !      if (ee.eq.0.d0) cycle
  !      ivdw = ivdw + 1
  !      write(101,'(2a8,a4,2f12.6)') vdwtype(itype1), vdwtype(itype2), wff_type, ee, ss
  !   end do
  !end do
  !print*, 'vdw = ', ivdw
  do itype1 = 1, Ntype
     do itype2 = itype1, Ntype
        if ((vdwtype(itype1).eq."Ni".and.vdwtype(itype2).eq."oz").or.(vdwtype(itype1).eq."oz".and.vdwtype(itype2).eq."Ni"))  then
           write(101,'(2a8,a4,4f18.6)') vdwtype(itype1), vdwtype(itype2), 'bucd', 80851.8, 0.261444, 18881.4, 0.68789
        else
        ss = sigma(itype1)+sigma(itype2)
        ee = sqrt(eps(itype1)*eps(itype2))
        if (ee.eq.0.d0) cycle
        ivdw = ivdw + 1
        write(101,'(2a8,a4,2f12.6)') vdwtype(itype1), vdwtype(itype2), wff_type, ee, ss
     end if
    enddo
  end do
  print*, 'vdw = ', ivdw+1
  write(101,'(a5)') 'CLOSE'

end program build_ff


!=======================================================================================

subroutine cut(atype1, atype2, cutoff)

  character(2) :: atype1, atype2
  real(8) :: cutoff

  cutoff = 1.65d0

  if (atype1.eq."Ni".or.atype2.eq."Ni") cutoff = 2.4d0

  !if ((atype1.eq."c1".and.atype2.eq."o") .or. &
  !    (atype1.eq."o".and.atype2.eq."c1"))  cutoff = 1.3d0


end subroutine cut

!======================================================================================

subroutine invert(a,b,d)

   implicit none

   real(8) a(9),b(9),d,r

   !! calculate adjoint matrix
   b(1)=a(5)*a(9)-a(6)*a(8)
   b(2)=a(3)*a(8)-a(2)*a(9)
   b(3)=a(2)*a(6)-a(3)*a(5)
   b(4)=a(6)*a(7)-a(4)*a(9)
   b(5)=a(1)*a(9)-a(3)*a(7)
   b(6)=a(3)*a(4)-a(1)*a(6)
   b(7)=a(4)*a(8)-a(5)*a(7)
   b(8)=a(2)*a(7)-a(1)*a(8)
   b(9)=a(1)*a(5)-a(2)*a(4)

   !! calculate determinant
   d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
   r=0.d0
   if(abs(d).gt.0.d0)r=1.d0/d

   !! complete inverse matrix
   b(1)=r*b(1)
   b(2)=r*b(2)
   b(3)=r*b(3)
   b(4)=r*b(4)
   b(5)=r*b(5)
   b(6)=r*b(6)
   b(7)=r*b(7)
   b(8)=r*b(8)
   b(9)=r*b(9)

end

!=======================================================================

subroutine calculate_cell(box,angle,cell)

  implicit none

  real(8) :: box(3), angle(3), cell(9)

  real(8) :: pi, vv, deg2rad
  real(8) :: sin_alpha, cos_alpha
  real(8) :: sin_beta, cos_beta
  real(8) :: sin_gamma, cos_gamma

  pi = acos(-1.d0)

  deg2rad = pi / 180.d0

  sin_alpha = sin(angle(1) * deg2rad)
  cos_alpha = cos(angle(1) * deg2rad)

  sin_beta = sin(angle(2) * deg2rad)
  cos_beta = cos(angle(2) * deg2rad)

  sin_gamma = sin(angle(3) * deg2rad)
  cos_gamma = cos(angle(3) * deg2rad)
  
  vv = sqrt(1.d0 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 &
           + 2*cos_alpha*cos_beta*cos_gamma)

  cell(1) = box(1) 
  if (abs(cell(1)).lt.1d-6) cell(1) = 0.d0
  cell(2) = 0.d0   
  cell(3) = 0.d0   

  cell(4) = box(2) * cos_gamma
  if (abs(cell(4)).lt.1d-6) cell(4) = 0.d0
  cell(5) = box(2) * sin_gamma
  if (abs(cell(5)).lt.1d-6) cell(5) = 0.d0
  cell(6) = 0.d0   

  cell(7) = box(3) * cos_beta
  if (abs(cell(7)).lt.1d-6) cell(7) = 0.d0
  cell(8) = box(3) * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
  if (abs(cell(8)).lt.1d-6) cell(8) = 0.d0
  cell(9) = box(3) * vv / sin_gamma
  if (abs(cell(9)).lt.1d-6) cell(9) = 0.d0

end

!=======================================================================

