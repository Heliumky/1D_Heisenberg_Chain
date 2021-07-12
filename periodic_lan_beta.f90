!This is a beta version program for 1D Heisenberg's model by lanczos method.
!Author: Jerry Chen (Qiancan Chen).

!abbreviation words
!###############################################################
!de_ : decimal.
!coeff : coefficent of wavefunction.
!bi_ : binary.
!ba_:backup of. ex: ba_bi_uparr: as same as bi_uparr, it is a chache of bi_uparr.
!config : configurations.
!arr: array.
!tran_ : transformate.
!lan_: lanczos
!hamil: hamiltonian
!ladd: ladder, [(s+)i * (s-)i+1 + (s-)i * (s+)i+1]
!z: sz, szi * szi+1
!pos: position in array or index in array.
!dic_: dictionary of index of de_uparr. 
!it_ : iterated
!
!###############################################################


!The introduction of parameter.
!###############################################################
!n: sites of system.
!n_up: how many up states of the sites.
!de_up: the up state is in decimal form. 
!com_up: it is a number that mean how many configurations of the C n takes n_up.
!de_n: any state in decimal form.
!up_pos: the position of up state is in de_uparr.
!de_up_ladd:
!de_up_z:
!jp: the coeff of  ladder part.
!jz: the coeff of sz part. 
!de_up_ladd: ladder part of de_up.
!de_up_z: sz part of de_up.
!norm_arr: normalized coeff_arr
!out_coeff_arr: the array output by Hamiltonian act on coeff_arr. 
!beta_arr: contains the lan_vec diagonal elements of the tridiagonal lanczos matrix.
!alpha_arr: contains the n-1 subdiagonal elements of the tridiagonal lanczos matrix.
!###############################################################

!The introduction of array.
!###############################################################
!bi_arr: any state is saved as binary form.
!bi_uparr: the up state is saved as binary form.
!coeff_arr: the coeff of up state which is used to superposition to get wf.
!de_uparr: the decimal form state by up states.
!###############################################################

!The introdution of subroutine.
!###############################################################
!config_up(n,n_up,com_up): to calculate how many config of n_up.
!tran_de_bi(n,de_n,bi_arr): transform de_n to bi_arr.
!tran_bi_de(n,bi_arr,de_n): transform bi_arr to de_n.
!filter_uparr(n,n_up,bi_arr,bi_uparr): choose bi_uparr from bi_arr.
!hamilt(n,com_up,jp,jz,de_uparr,dic_arr,coeff_arr): hamiltonian of 1d heisenberg's model.
!dstev('N', lan_vec, alpha_arr, beta_arr, lan_eigvec, ldz, work, info): sovling tridiagonal matrix eigenvalue problems by intel mkl. more info: https://software.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-fortran/top/lapack-routines/lapack-least-squares-and-eigenvalue-problem-routines/lapack-least-squares-and-eigenvalue-problem-driver-routines/symmetric-eigenvalue-problems-lapack-driver-routines/stev.html
!###############################################################


program lan_method
  implicit none
  integer(4),parameter :: n = 30, n_up = n/2, max_lan_vec = 100
  integer(4) :: de_up, com_up, de_n, up_pos, lan_vec, i, ldz, info
  real(8),parameter :: jp = 0.5d0, jz = 1.d0
  real(8) :: lmin_eigval = 5.d0, min_eigval = 3.d0, diff_e = 3.d0, invrs_norm
  real(8),dimension(max_lan_vec,max_lan_vec) :: lan_eigvec 
  real(8),dimension(2*max_lan_vec-2) :: work
  integer(4),dimension(n) :: bi_arr, bi_uparr
  integer(4),dimension(2**n) :: dic_arr
  real(8),allocatable :: coeff_arr(:), out_coeff_arr(:), tuda_coeff_arr(:)
  real(8),allocatable :: beta_arr(:), alpha_arr(:), eigvec_arr(:),it_norm_arr(:), norm_arr(:),eigval_arr(:),ba_beta_arr(:)
  integer(4),allocatable :: de_uparr(:)
  real(8), external :: dnrm2, ddot
  call config_up(n,n_up,com_up)
  allocate(coeff_arr(com_up))
  call random_number(coeff_arr)
  allocate(beta_arr(max_lan_vec-1))
  allocate(ba_beta_arr(max_lan_vec-1))
  allocate(alpha_arr(max_lan_vec))
  allocate(eigval_arr(max_lan_vec))
  allocate(tuda_coeff_arr(com_up))
  allocate(out_coeff_arr(com_up))
  allocate(it_norm_arr(com_up))
  allocate(eigvec_arr(com_up))
  allocate(norm_arr(com_up))
  allocate(de_uparr(com_up))
  de_n  = 2**n -1 
  up_pos = 0
!built dictionary and de_uparr
!###############################################################
  do i = 0, de_n
    call tran_de_bi(n,i,bi_arr)
    call filter_uparr(n,n_up,bi_arr,bi_uparr)
    call tran_bi_de(n,bi_uparr,de_up)
    if (i == de_up .and. de_up /= 0) then
      up_pos = up_pos + 1
      dic_arr(i) = up_pos
      de_uparr(up_pos) = de_up
    end if
  end do
!write(*,*) dic_arr
!###############################################################
  
!save alpah_arr, beta_arr.
!###############################################################
!  open(UNIT=10, FILE='alpha.txt')
!  open(UNIT=11, FILE='beta.txt')
!###############################################################

!lanczos method
!###############################################################
  lan_vec = 2  
  invrs_norm  =dnrm2(com_up,coeff_arr,1)
  invrs_norm =1.d0/ invrs_norm
  it_norm_arr=coeff_arr* invrs_norm
  call hamilt(n,com_up,jp,jz,de_uparr,dic_arr,it_norm_arr,out_coeff_arr)
  alpha_arr(1) = ddot(com_up,it_norm_arr,1,out_coeff_arr,1)
  tuda_coeff_arr = out_coeff_arr - alpha_arr(1) * it_norm_arr
  do while (.true.)
    beta_arr(lan_vec-1) = dnrm2(com_up,tuda_coeff_arr,1)
    invrs_norm=1.d0/beta_arr(lan_vec-1)
    norm_arr = tuda_coeff_arr* invrs_norm
    call hamilt(n,com_up,jp,jz,de_uparr,dic_arr,norm_arr,out_coeff_arr)
    alpha_arr(lan_vec) = ddot(com_up,norm_arr,1,out_coeff_arr,1)
    tuda_coeff_arr = out_coeff_arr - alpha_arr(lan_vec)*norm_arr - beta_arr(lan_vec-1)*it_norm_arr
    it_norm_arr = norm_arr
    lmin_eigval = min_eigval
  if (lan_vec > 39) then
    eigval_arr = alpha_arr 
    ba_beta_arr = beta_arr
    ldz = lan_vec
    call dstev('V', lan_vec, eigval_arr, ba_beta_arr, lan_eigvec, ldz, work, info)
    min_eigval = eigval_arr(1)
    diff_e = abs(abs(min_eigval) - abs(lmin_eigval))
   if (diff_e < 1.d-16 ) then
     exit
     end if
   end if
   lan_vec = lan_vec + 1
  end do 
write(*,*) lan_vec, min_eigval 
!###############################################################
end program lan_method

subroutine config_up(n,n_up,com_up)
  implicit none
  integer(4) :: n, n_up, com_up, i
  real(8) :: c
  real(8),dimension(n_up) :: p_arr, q_arr, r_arr
  do i = 1, n_up
    p_arr(i) = n - i + 1
  end do 
  do i = 1, n_up
    q_arr(i) = n_up - i + 1
  end do
  r_arr = p_arr/q_arr
  c = product(r_arr)
  if (abs(c - int(c,kind = 4)) >=0.5) then
    c = int(c,kind = 4) + 1
  else
    c = c
  end if
  com_up = int(c,kind = 4)
  return 
end subroutine config_up

subroutine tran_de_bi(n,de_n,bi_arr)
  implicit none
  integer(4) :: de_n,ba_de_n
  integer(4) :: n, j
  integer(4),dimension(n) :: bi_arr
  bi_arr(:) = 0
  ba_de_n = de_n
  j = 0 
  do while (ba_de_n>0)
    j = j +1
    bi_arr(n + 1 - j) = mod(ba_de_n, 2) 
    ba_de_n = ba_de_n / 2
  end do
end subroutine tran_de_bi

subroutine tran_bi_de(n,bi_arr,de_n)
  implicit none
  integer(4) :: n, de_n, i, a 
  integer(4),dimension(n) :: bi_arr
  de_n = 0
  a = 1
  do i = 1, n
   ! de_n = de_n + bi_arr(n-i+1)*(2**(i-1))
   de_n = de_n + bi_arr(n-i+1) * a 
   a = a + a 
  end do
  return
end subroutine tran_bi_de

subroutine filter_uparr(n,n_up,bi_arr,bi_uparr)
  implicit none
  integer(4) :: n, n_up 
  integer(4),dimension(n) :: bi_arr, bi_uparr
  if (sum(bi_arr) == n_up) then
    bi_uparr = bi_arr
  end if
  return
end subroutine filter_uparr

subroutine hamilt(n,com_up,jp,jz,de_uparr,dic_arr,coeff_arr,out_coeff_arr)
  implicit none
  integer(4) :: n, com_up, de_up_ladd, de_up_z, i, j, k
  real(8) :: jp, jz
  integer(4),dimension(n) :: bi_uparr, ba_bi_uparr
  integer(4),dimension(2**n) :: dic_arr
  integer(4),dimension(com_up) :: de_uparr
  real(8), dimension(com_up) :: coeff_arr,out_coeff_arr! , ladd_arr!, pauliz_arr
!  ladd_arr(:) = 0.d0
!  pauliz_arr(:) = 0.d0
  out_coeff_arr(:)=0.d0
  do i = 1, com_up
    call tran_de_bi(n,de_uparr(i),bi_uparr)
   ! write(*,*)  n
   do j = 1, n
      k =mod((j + 1),n + 1)
      if (k ==0) then
        k = 1
      end if
     ! write(*,*) k
!##############################################################
!     if (bi_uparr(j)==0 .and. bi_uparr(k)==0) then
!       call tran_bi_de(n,bi_uparr,de_up_z)
!       pauliz_arr(dic_arr(de_up_z)) = 0.25d0 * coeff_arr(i) + pauliz_arr(dic_arr(de_up_z))
!     
!     else if (bi_uparr(j)==1 .and. bi_uparr(k)==1) then
!       call tran_bi_de(n,bi_uparr,de_up_z)
!       pauliz_arr(dic_arr(de_up_z)) = 0.25d0 * coeff_arr(i) + pauliz_arr(dic_arr(de_up_z))
!###!###########################################################
     if (bi_uparr(j)==0 .and. bi_uparr(k)==0) then
       call tran_bi_de(n,bi_uparr,de_up_z)
       out_coeff_arr(dic_arr(de_up_z)) = 0.25d0*jz * coeff_arr(i) + out_coeff_arr(dic_arr(de_up_z))
     
     else if (bi_uparr(j)==1 .and. bi_uparr(k)==1) then
       call tran_bi_de(n,bi_uparr,de_up_z)
       out_coeff_arr(dic_arr(de_up_z)) = 0.25d0*jz * coeff_arr(i) + out_coeff_arr(dic_arr(de_up_z))
!###!###########################################################

     else if (bi_uparr(j)==0 .and. bi_uparr(k)==1) then
       ba_bi_uparr = bi_uparr
       ba_bi_uparr(j) = 1
       ba_bi_uparr(k) = 0
      ! write(*,*) de_up, ba_bi_uparr
       call tran_bi_de(n,bi_uparr,de_up_z)
       call tran_bi_de(n,ba_bi_uparr,de_up_ladd)
       !ladd_arr(dic_arr(de_up_ladd)) = coeff_arr(i) + ladd_arr(dic_arr(de_up_ladd))
       out_coeff_arr(dic_arr(de_up_ladd)) = jp* coeff_arr(i) + out_coeff_arr(dic_arr(de_up_ladd))
       !pauliz_arr(dic_arr(de_up_z)) = (-0.25d0) * coeff_arr(i) + pauliz_arr(dic_arr(de_up_z))
       out_coeff_arr(dic_arr(de_up_z)) = (-0.25d0)*jz * coeff_arr(i) + out_coeff_arr(dic_arr(de_up_z))

     else if (bi_uparr(j)==1 .and. bi_uparr(k)==0) then
       ba_bi_uparr = bi_uparr
       ba_bi_uparr(j) = 0
       ba_bi_uparr(k) = 1
       call tran_bi_de(n,bi_uparr,de_up_z)
       call tran_bi_de(n,ba_bi_uparr,de_up_ladd)
      ! write(*,*) de_up, ba_bi_uparr
       !ladd_arr(dic_arr(de_up_ladd)) = coeff_arr(i) + ladd_arr(dic_arr(de_up_ladd))
       out_coeff_arr(dic_arr(de_up_ladd)) = jp*coeff_arr(i) + out_coeff_arr(dic_arr(de_up_ladd))
       !pauliz_arr(dic_arr(de_up_z)) = (-0.25d0) * coeff_arr(i) + pauliz_arr(dic_arr(de_up_z))
       out_coeff_arr(dic_arr(de_up_z)) = (-0.25d0)*jz * coeff_arr(i) + out_coeff_arr(dic_arr(de_up_z))
     end if

!     write(*,*) de_up_ladd, de_up_z
!##############################################################
   end do
 end do
!out_coeff_arr = jz * pauliz_arr + jp * ladd_arr
 ! write(*,*) coeff_arr
end subroutine hamilt
