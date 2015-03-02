program build_ff

  !! Build MOF force field for DL_POLY
 
  implicit none

  integer, parameter :: Ndim=3, Nmax=100000, icount=1, izero=0, ione=1, itwo=2, ithree=3

  real(8), parameter :: zero=0.d0

  !-----------------------------------------------------------------------

  character(20) :: pdb, topology, output

  !------------------------------------------------------------------------

  character(4) :: wff_type, uname1, uname2
  character(10) :: ff_key
  character(10), allocatable :: vdwtype(:), atname(:), ff_type(:)
  character(10), allocatable :: uname(:), utype(:)
  character(10), allocatable :: atype(:), atype1(:), atype2(:), atype3(:), atype4(:)

  integer :: iatom, Natom, iatom1, iatom2, iatom3, iatom4, ibond, iangle, idihedral, ivdw, icheck
  integer :: itype, Ntype, itmp, iuni, Nuni, itype1, itype2, isum, idiff, Nfac, wifac, wn, wifacok
  integer :: Nbond, Nangle, Ndihe, Nvdw, Nnew, Ncheck,nb
  integer :: Nbond_fit, Nangle_fit, Ndihe_fit, Nimp_fit
  integer, allocatable :: iat1(:), iat4(:), n(:), ifac(:)
  
  real(8) :: charge_tot
  real(8) :: box_x, box_y, box_z
  real(8) :: wk, k1, k1_new, k2_new, weq, eq1, eq1_new, eq2_new, walpha, alpha1_new, alpha2_new, wphi
  real(8) :: wcc1, wcc2, wcc3
  real(8) :: ee, ss, e1, e2, s1, s2
  real(8) :: r2, r, cutoff, cutoff_OZnO
  real(8) :: dx(Ndim), box(Ndim)
  real(8), allocatable :: pos(:,:), sigma(:), eps(:), charge(:), mass(:), e14(:), vdw14(:)
  real(8), allocatable :: ucharge(:), umass(:),refcut(:)
  real(8), allocatable :: k(:), eq(:), alpha(:), phi(:), cc1(:), cc2(:), cc3(:)
  
  character(10), allocatable :: atomtype(:), resname(:), ref1(:),ref2(:)
  integer :: ierr
  integer, allocatable :: resnum(:)

  !-------------------------------------------------------------------------

  character(10) :: junk
  integer :: ijunk, ncount, i
  real(8) :: rjunk
  
  !-------------------------------------------------------------------------

  namelist/input/Nuni, Natom, Nbond, Nangle, Ndihe, Nbond_fit, Nangle_fit, Ndihe_fit, Nimp_fit, Nvdw, Nfac, wifacok
  
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
  !open(unit=102, file='CONFIG_TMP')

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
  allocate(ref1(Nuni))
  allocate(ref2(Nuni))
  allocate(refcut(Nuni))
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

  !-------------------------------------------------------------------------

  !! Read input from Gaussian pdb.
  read(11,*) junk, box(1), box(2), box(3)
  do iatom = 1, Natom
     !read(11,*) junk, ijunk, atname(iatom), junk, ijunk, pos(:,iatom)
     read(11,*) junk, ijunk, atname(iatom), ijunk, pos(:,iatom)
     !write(6,*) junk, ijunk, atname(iatom), ijunk, pos(:,iatom)
  end do
  close(11)
  
  !! Atoms

  read(2,*) ff_key, Nuni  
  do iuni = 1, Nuni
     read(2,*) uname(iuni), utype(iuni), ucharge(iuni), umass(iuni)
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
           charge(iatom) = ucharge(iatom) 
           mass(iatom) = umass(iuni) 
           ncount = ncount + 1
           exit
        end if
     end do
     charge_tot = charge_tot + charge(iatom)
     write(101,'(a8,2f12.7,i8)') atype(iatom), mass(iatom), charge(iatom), icount
     !write(102,'(a8,i10,2f12.6)') atype(iatom), iatom, mass(iatom), charge(iatom)
     !write(102,'(3f15.6)') pos(:,iatom)
  end do
  read(2,*)
  print*, 'atom = ', ncount
  !print*, 'total charge = ', charge_tot

  !! Bonds.
  ibond = 0
  read(2,*) ff_key, Ntype

  do itype = 1, Ntype
     read(2,*) atype1(itype), atype2(itype), ff_type(itype), k(itype), eq(itype), cc1(itype), cc2(itype)
  end do
  nb=Ntype

  read(2,*)
  do itype = 1, Nbond_fit
     read(21,*) k(itype), eq(itype), cc1(itype), cc2(itype)
  end do
  write(101,'(a10,i10)') ff_key, Nbond
  do itype = 1, Ntype
     !! Identify 1st atom.
     do iatom1 = 1, Natom
        !! Identify 2nd atom.
        do iatom2 = iatom1+1, Natom
           dx(:) = pos(:,iatom2)-pos(:,iatom1)
           dx(:) = dx(:) - box(:)*nint(dx(:)/box(:))
           r2 = sum(dx(:)**2)
           r = sqrt(r2)
           if ((atype(iatom1).eq.atype1(itype) .and. atype(iatom2).eq.atype2(itype)) .or. &
               (atype(iatom1).eq.atype2(itype) .and. atype(iatom2).eq.atype1(itype))) then
               call cut(atype(iatom1), atype(iatom2),nb,cutoff)
               if (r.lt.cutoff) then
                  wff_type = ff_type(itype)
                  wk = k(itype)
                  weq = eq(itype)
                  wcc1 = cc1(itype)
                  wcc2 = cc2(itype)
                  ibond = ibond + 1
                  write(101,'(a4,2i8,4f15.6)') wff_type, iatom1, iatom2, wk, weq, wcc1, wcc2
                !  print *, r,atype(iatom1),'  ',atype(iatom2),iatom1,iatom2
               end if
           end if
        end do
     end do
  end do
  print*, 'bonds = ', ibond
  !write(101,*) 'bonds  ',ibond


  !! Angles.
  iangle = 0
  read(2,*) ff_key, Ntype
  !print*, ff_key, Ntype
  do itype = 1, Ntype
     read(2,*) atype1(itype), atype2(itype), atype3(itype), ff_type(itype), k(itype), eq(itype), cc1(itype), cc2(itype)
  end do
  read(2,*)
  do itype = 1, Nangle_fit
     read(21,*) k(itype), eq(itype), cc1(itype), cc2(itype)
  end do
  write(101,'(a10,i10)') ff_key, Nangle

  
  do itype = 1, Ntype
     !! End atoms equal.
     if (atype1(itype).eq.atype3(itype)) then
        !! Identify central (2nd) atom.
        do iatom2 = 1, Natom
           if (atype(iatom2).eq.atype2(itype)) then
              !! Identify 1st atom.
              do iatom1 = 1, Natom
                 if (atype(iatom1).eq.atype1(itype) .and. iatom1.ne.iatom2) then
                    dx(:) = pos(:,iatom2)-pos(:,iatom1)
                    dx(:) = dx(:) - box(:)*nint(dx(:)/box(:))
                    r2 = sum(dx(:)**2)
                    r = sqrt(r2)
                   call cut(atype(iatom1), atype(iatom2),nb,cutoff)
                    if (r.lt.cutoff) then
                       !! Identify 3rd atom.
                       do iatom3 = iatom1+1, Natom
                          if (atype(iatom3).eq.atype3(itype) .and. iatom3.ne.iatom2 .and. iatom3.ne.iatom1) then
                             dx(:) = pos(:,iatom3)-pos(:,iatom2)
                             dx(:) = dx(:) - box(:)*nint(dx(:)/box(:))
                             r2 = sum(dx(:)**2)
                             r = sqrt(r2)
                             call cut(atype(iatom1), atype(iatom2),nb,cutoff)
                             if (r.lt.cutoff) then
                                wff_type = ff_type(itype)
                                wk = k(itype)
                                weq = eq(itype)
                                wcc1 = cc1(itype)
                                wcc2 = cc2(itype)
                                iangle = iangle + 1
                                write(101,'(a4,3i8,4f15.6)') wff_type, iatom1, iatom2, iatom3, wk, weq, wcc1, wcc2
                             end if
                          end if
                       end do
                    end if
                 end if
              end do
           end if
        end do
     else
        do iatom2 = 1, Natom
           if (atype(iatom2).eq.atype2(itype)) then
              !! Identify 1st atom.
              do iatom1 = 1, Natom
                 if (atype(iatom1).eq.atype1(itype) .and. iatom1.ne.iatom2) then
                    dx(:) = pos(:,iatom2)-pos(:,iatom1)
                    dx(:) = dx(:) - box(:)*nint(dx(:)/box(:))
                    r2 = sum(dx(:)**2)
                    r = sqrt(r2)
                    call cut(atype(iatom1), atype(iatom2),nb,cutoff)
                    if (r.lt.cutoff) then
                       !! Identify 3rd atom.
                       do iatom3 = 1, Natom
                          if (atype(iatom3).eq.atype3(itype) .and. iatom3.ne.iatom2 .and. iatom3.ne.iatom1) then
                             dx(:) = pos(:,iatom3)-pos(:,iatom2)
                             dx(:) = dx(:) - box(:)*nint(dx(:)/box(:))
                             r2 = sum(dx(:)**2)
                             r = sqrt(r2)
                             call cut(atype(iatom1), atype(iatom2),nb,cutoff)
                             if (r.lt.cutoff) then
                                iangle = iangle + 1
                                wff_type = ff_type(itype)
                                wk = k(itype)
                                weq = eq(itype)
                                wcc1 = cc1(itype)
                                wcc2 = cc2(itype)
                                write(101,'(a4,3i8,4f15.6)') wff_type, iatom1, iatom2, iatom3, wk, weq, wcc1, wcc2
                             end if
                          end if
                       end do
                    end if
                 end if
              end do
           end if
        end do
     end if
  end do
  print*, 'angles = ', iangle
  !write(101,*) 'angles  ', iangle


  !! Dihedrals.
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
  do itype = 1, Ntype
     !! End atoms equal.
     if (atype1(itype).eq.atype4(itype)) then
        !! Identify 2nd atom.
        do iatom2 = 1, Natom
           if (atype(iatom2).eq.atype2(itype)) then
              !! Identify 3rd atom.
              do iatom3 = 1, Natom
                 if (atype(iatom3).eq.atype3(itype) .and. iatom3.ne.iatom2) then
                    dx(:) = pos(:,iatom3)-pos(:,iatom2)
                    dx(:) = dx(:) - box(:)*nint(dx(:)/box(:))
                    r2 = sum(dx(:)**2)
                    r = sqrt(r2)
                    call cut(atype(iatom1), atype(iatom2),nb,cutoff)
                    if (r<cutoff) then
                       !! Identify 1st atom.
                       do iatom1 = 1, Natom
                          if (atype(iatom1).eq.atype1(itype) .and. iatom1.ne.iatom2 .and. iatom1.ne.iatom3) then
                             dx(:) = pos(:,iatom2)-pos(:,iatom1)
                             dx(:) = dx(:) - box(:)*nint(dx(:)/box(:))
                             r2 = sum(dx(:)**2)
                             r = sqrt(r2)
                            call cut(atype(iatom1), atype(iatom2),nb,cutoff)
                             if (r<cutoff) then
                                !! Identify 4th atom.
                                do iatom4 = iatom1+1, Natom
                                   if (atype(iatom4).eq.atype4(itype) .and. iatom4.ne.iatom1 .and. iatom4.ne.iatom2 .and. iatom4.ne.iatom3) then
                                      dx(:) = pos(:,iatom4)-pos(:,iatom3)
                                      dx(:) = dx(:) - box(:)*nint(dx(:)/box(:))
                                      r2 = sum(dx(:)**2)
                                      r = sqrt(r2)
                                      call cut(atype(iatom1), atype(iatom2),nb,cutoff)
                                      if (r<cutoff) then
                                         idihedral = idihedral + 1
                                         iat1(idihedral) = iatom1
                                         iat4(idihedral) = iatom4
                                         wff_type = ff_type(itype)
                                         wk = k(itype)
                                         wphi = phi(itype)
                                         wn = n(itype)
                                         wifac = ifac(itype)
                                         do icheck = 1, idihedral-1
                                            if ((iat1(icheck).eq.iatom1 .and. iat4(icheck).eq.iatom4) .or. (iat1(icheck).eq.iatom4 .and. iat4(icheck).eq.iatom1)) then
                                               wifac = izero
                                            end if
                                         end do
                                         write(101,'(a4,4i8,f15.6,f15.6,i8,2f15.6)') wff_type, iatom1, iatom2, iatom3, iatom4, wk, wphi, abs(wn), e14(wifac), vdw14(wifac)
                                      end if
                                   end if
                                end do
                             end if
                          end if
                       end do
                    end if
                 end if
              end do
           end if
        end do
     else
        !! Identify 2nd atom.
        do iatom2 = 1, Natom
           if (atype(iatom2).eq.atype2(itype)) then
              !! Identify 3rd atom.
              do iatom3 = 1, Natom
                 if (atype(iatom3).eq.atype3(itype) .and. iatom3.ne.iatom2) then
                    dx(:) = pos(:,iatom3)-pos(:,iatom2)
                    dx(:) = dx(:) - box(:)*nint(dx(:)/box(:))
                    r2 = sum(dx(:)**2)
                    r = sqrt(r2)
                    call cut(atype(iatom1), atype(iatom2),nb,cutoff)
                    if (r<cutoff) then
                       !! Identify 1st atom.
                       do iatom1 = 1, Natom
                          if (atype(iatom1).eq.atype1(itype) .and. iatom1.ne.iatom2.and. iatom1.ne.iatom3) then
                             dx(:) = pos(:,iatom2)-pos(:,iatom1)
                             dx(:) = dx(:) - box(:)*nint(dx(:)/box(:))
                             r2 = sum(dx(:)**2)
                             r = sqrt(r2)
                             call cut(atype(iatom1), atype(iatom2),nb,cutoff)
                             if (r<cutoff) then
                                !! Identify 4th atom.
                                do iatom4 = 1, Natom
                                   if (atype(iatom4).eq.atype4(itype) .and.iatom4.ne.iatom1 .and. iatom4.ne.iatom2 .and. iatom4.ne.iatom3) then
                                      dx(:) = pos(:,iatom4)-pos(:,iatom3)
                                      dx(:) = dx(:) - box(:)*nint(dx(:)/box(:))
                                      r2 = sum(dx(:)**2)
                                      r = sqrt(r2)
                                      call cut(atype(iatom1), atype(iatom2),nb,cutoff)
                                      if (r<cutoff) then
                                         idihedral = idihedral + 1
                                         iat1(idihedral) = iatom1
                                         iat4(idihedral) = iatom4
                                         wff_type = ff_type(itype)
                                         wk = k(itype)
                                         wphi = phi(itype)
                                         wn = n(itype)
                                         wifac = ifac(itype)
                                         do icheck = 1, idihedral-1
                                            if ((iat1(icheck).eq.iatom1 .and. iat4(icheck).eq.iatom4) .or. (iat1(icheck).eq.iatom4 .and. iat4(icheck).eq.iatom1)) then
                                               wifac = izero
                                            end if
                                         end do
                                         write(101,'(a4,4i8,f15.6,f15.6,i8,2f15.6)') wff_type, iatom1, iatom2, iatom3, iatom4, wk, wphi, abs(wn), e14(wifac), vdw14(wifac)
                                      end if
                                   end if
                                end do
                             end if
                          end if
                       end do
                    end if
                 end if
              end do
           end if
        end do
     end if
  end do
  print*, 'dihedrals = ', idihedral
  !write(101,*) 'dihedral  ', idihedral

  !! Improper.
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
  do itype = 1, Ntype
     !! Identify central atom.
     do iatom3 = 1, Natom
        if (atype(iatom3).eq.atype3(itype)) then
           !! Identify 1st atom.
           do iatom1 = 1, Natom
              if (atype(iatom1).eq.atype1(itype) .and. iatom1.ne.iatom3) then
                 dx(:) = pos(:,iatom1)-pos(:,iatom3)
                 dx(:) = dx(:) - box(:)*nint(dx(:)/box(:))
                 r2 = sum(dx(:)**2)
                 r = sqrt(r2)
                 call cut(atype(iatom1), atype(iatom2),nb,cutoff)
                 if (r<cutoff) then
                    !! Identify 2nd atom.
                    do iatom2 = 1, Natom
                       if (atype(iatom2).eq.atype2(itype) .and. iatom2.ne.iatom3 .and. iatom2.ne.iatom1) then
                          dx(:) = pos(:,iatom2)-pos(:,iatom3)
                          dx(:) = dx(:) - box(:)*nint(dx(:)/box(:))
                          r2 = sum(dx(:)**2)
                          r = sqrt(r2)
                         call cut(atype(iatom1), atype(iatom2),nb,cutoff)
                          if (r<cutoff) then
                             !! Identify 4th atom.
                             do iatom4 = 1, Natom
                                if (atype(iatom4).eq.atype4(itype) .and. iatom4.ne.iatom3 .and. iatom4.ne.iatom1 .and. iatom4.ne.iatom2) then
                                   dx(:) = pos(:,iatom4)-pos(:,iatom3)
                                   dx(:) = dx(:) - box(:)*nint(dx(:)/box(:))
                                   r2 = sum(dx(:)**2)
                                   r = sqrt(r2)
                                   call cut(atype(iatom1), atype(iatom2),nb,cutoff)
                                   if (r<cutoff) then
                                      idihedral = idihedral + 1
                                      wff_type = ff_type(itype)
                                      wk = k(itype)
                                      wphi = phi(itype)
                                      wn = n(itype)
                                      wifac = izero          
                                      write(101,'(a4,4i8,f15.6,f15.6,i8,2f15.6)') wff_type, iatom1, iatom2, iatom3, iatom4, wk, wphi, abs(wn), e14(wifac), vdw14(wifac)
                                   end if
                                end if
                             end do
                          end if
                       end if
                    end do
                 end if
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
  do itype1 = 1, Ntype
     do itype2 = itype1, Ntype
        ss = sigma(itype1)+sigma(itype2)
        ee = sqrt(eps(itype1)*eps(itype2))
        if (ee.eq.0.d0) cycle
        ivdw = ivdw + 1
        write(101,'(2a8,a4,2f12.6)') vdwtype(itype1), vdwtype(itype2), wff_type, ee, ss
     end do
  end do
  print*, 'vdw = ', ivdw

  write(101,'(a5)') 'CLOSE'

end program build_ff


!=======================================================================================
subroutine cut(atype1, atype2,ntypes,cutoff)

  character(2) :: atype1, atype2
  character(2),allocatable :: ref_1(:), ref_2(:)
  real(8), allocatable :: ref_cut(:)
  integer :: ntypes
  real(8) :: cutoff

  allocate(ref_1(ntypes))
  allocate(ref_2(ntypes))
  allocate(ref_cut(ntypes))

  open(11,file='cutoffs.dat')
  do itype = 1, ntypes
     read(11,*) ref_1(itype),ref_2(itype),ref_cut(itype)
  end do
  close(11)
  
  do i=1,ntypes
     if ((atype1.eq.ref_1(i).and.atype2.eq.ref_2(i)) .or. (atype1.eq.ref_2(i).and.atype2.eq.ref_1(i)))  cutoff = ref_cut(i)
  enddo


end subroutine cut
!======================================================================================
