
c-----------------------------------------------------------------------
c
c Mar 14, 2016, Jay
c
c * I have fixed some bugs that have mostly do with with print-outs
c   of diagnostic items when the items don't get initialized properly.
c  
c
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c
c Nov 2, 2012: Jay
c
c I looked over all the changes that had been made since I last gave Jeff
c the code (pretty much those made on Oct 30 and Oct 31) and have the
c following comments about things that were done and things that still
c need to be done
c
c GENERAL CHANGES
c
c 1) renamed all instances of "scl" to "rat" for consistency
c    perhaps I should have stayed with "scl" since it's more
c    descriptive... we can change that back if we want
c 2) the warm start doesn't work anymore; need to read in meteor 
c    info from the save_mod() file ; will that be enough
c 3) fixed a bug in read_mod() ; it didn't know ahead of time 
c    what NDIMs was ; it should use _NDIM_ to reset unim to .false. 
c    (start off flagging all exposures as unusable)
c 4) note any "print" statement is temporary  ; write(6,...) is
c    what should be used for permanent code
c
c
c LEINS/CONCERNS FOR THE FUTURE
c
c 1) find_rati_bar has a hardcoded size for the pixel array ; 
c    I can do better...
c 2) find_lsqtem  has harcoded limits for pixel index at the 
c    beginning
c 3) I still normalize each spectrum to have 50,000 counts average 
c    in the middle regions in the read-in routtine ; no good reason 
c    to do this
c 4) I still need to be more deliberate about how to use the flags:  
c    rej/uuu/met
c 5) in the template smoothing, I use hard numbers for the value 
c    of the template to determine how much smoothing to do; since
c    there is no longer any continuum normalization, these numbers
c    are not anchored to anything real.  This isn't really a huge 
c    deal now, since the continuum tends to be at 1.05 for TauCeti, 
c    but we may want to do some kind of rough continuum normalization
c    at the very least to define where to and not to do smoothing.
c    Of course, if we go the SVD route, then smoothing may be obviated
c 6) in the template normalization, I am currently just using the 
c    wavelengths between pixel 501 and 3501 in the first (?) observation; 
c    we may need a better way to specify this, perhaps at the very 
c    least with compiler defined variables and not hardcoded numbers
c
c------------------------------------------------------------------------






c-----------------------------------------------------------
c
c Oct 31, 2012:  Jay
c
c CHANGE: include meteor:  met_imn() contains
c         0 if no meteor impact at all
c         1 if the meteor is closer to this order than others
c         2 if the meteor crosses within 12 pixels
c         9 if this pixel is the center of the meteor crossing
c
c CHANGE: do better job wih smoother rati...
c         I now don't do a lsq spline, but rather
c         I just do a piecewise average; looks better;
c         I also read-in info for the meteor and
c         flag as bad the +/- 12 pixels about the
c         center ; the new piecewise avg does a better
c         job of identifying outliers for rejection
c
c CHANGE: the bulk velocities now ignore the pixels
c         that have been flagged-bad pixels, plus all
c         the pixels possibly impacted by the meteor 
c         (wider swath), and only fit between pixels
c         i=500 and i=3500.
c
c-----------------------------------------------------------
c
c Oct 30, 2012:  Jay
c
c CHANGE#1  
c
c    For some reason I have two arrays that do largely
c    the same thing, uuu_imn and rej_imn... in some
c    places they're set to the same thing.  I should
c    probably unify them, but need to be careful in
c    doing so that they may not be exactly the same.
c    the code is:
c        uuu_imn(i,m,n) = 0  <---- good pixel
c        uuu_imn(i,m,n) = 1  <---- rejected by sigma clipping
c        uuu_imn(i,m,n) = 2  <---- within the first 20 good pixels 
c                                  (_RATA_<= i <=_RATB_)
c        uuu_imn(i,m,n) = 3  <---- within the outer 20 good pixels 
c                                  (_RATY_<= i <=_RATZ_)
c        uuu_imn(i,m,n) = 4  <---- bad value entirely
c        uuu_imn(i,m,n) = 9  <---- rejected by meteor (not really anymore)
c        ---> still need to differentiate uuu/rat/met
c
c
c CHANGE#2
c
c    Fixed the lsqtem routine to go all the way to the end 
c    of the template; added a routine find_loctem_oe() do 
c    to it for the even/odd nodes separately, but this lowers 
c    the resolution
c
c CHANGE#3
c
c    New normalization; no longer try to find the continuum,
c    but rather just normalize such that the template over
c    the region of interest (central 3000 pixels) has an
c    average of 1.00000
c    ---> new routine cnorm_tem_bar() to do this
c    ---> old routine cnorm_tem() not used at moment
c
c CHANGE#4
c
c    When reading things in rdsk()
c    1) I now divide by the average 4th order blaze function 
c       from the B-star observations, order by order
c  ? 2) I then normalize every observation to have an 
c       average pixel value of 50000 in the middle 3000 pixels
c    3) just in case pixel i=4021 ever gets used, I set
c       its value to that of i=4020...
c
c
c THINGS TO THINK ABOUT
c
c 1) blaze normalization; look at the "rat" output in 
c    the mod file (column 8) ; there are some places where 
c    all the spectra follow along... strange; there is 
c    still some residual strange behavior at the i>4000 end
c 2) the meteor is starting to be an issue... and it could
c    well be limiting our ultimate precision
c 3) we should pay some attention to which pixels get used
c    in the "bulk" velocity fit...
c
c-----------------------------------------------------------

c----------------------------------------------------------
c
c need to think about improving the continuum normalization
c

c-------------------------------------------------------------------------------
c
c default flags
c
#define _DOFUDGE_      .true.
#define _DOVORB_       .false.
#define _OUTPUT_MOD_    0
#define _OUTPUT_LSF_    1
#define _OUTPUT_VEL_    2
#define _OUTPUT_TEM_    1
#define _OUTPUT_NFO_    0
#define _OUTPUT_NRM_    0
#define _OUTPUT_LIN_    0
#define _NITo_          0
#define _NITf_         10


c------------------------------------------------------------------------------
c
c may be needed if any raw disk files are read in ; I can worry about
c big/little endian issues
c
#define _LINUX_ .true.
 
 
c------------------------------------------------------------------------------
c
c speed of light in m/s
c
#define _CMPS_ 299 792 458
 
c------------------------------------------------------------------------------
c
c These variables regulate the pixel arrays; only the
c pixels between _IBEG_ and _IEND_ will be used ; for the
c Keck spectra, the 4021st pixel contains no real data
c
#define _IDIM_ 4021
#define _IBEG_ 0001
#define _IEND_ 4020
 
 
c------------------------------------------------------------------------------
c
c These compiler variables define the points used for
c fitting the continuum.  There are _RATN_ free parameters.
c The fit assumes a simple spline between point 2 and NP-1.
c The last points are free, and the slope is just linear
c between points (1,2) and points (NP-1,NP).  The outer two
c points by construction have *no* bearing on what happens
c between (2,NP-1).  The routine find_rat() both finds the
c ratio between observations and model and rejects pixels
c that are more than _FFREJ_SIGCLIP_ sigmas off.  If _RATN_ c
c is 103 and P2=21 and M2=4000, then there are 4000-21 = 3979
c pixels in the main part, and there are 101 gridpoints...
c these are spaced by 39.79 pixels.
c
#define _RATA_ 0001
#define _RATB_ 0021
c#define _RATN_  103 old: 2009/11/11
c#define _RATN_   20 old: 2009/11/15
c#define _RATN_   50 old: 2010/07/26
#define _RATN_  100
#define _RATY_ 4000
#define _RATZ_ 4020
 
 
c------------------------------------------------------------------------------
c
c number of images to be used ; control variable "n"
c
#define _NDIM_ 999
 
 
c------------------------------------------------------------------------------
c
c number of supersampled pixels (i=1 ---> o=101) ; this
c is the fundamental way the spectrum is upersampled, then
c binned into pixels; control variable "o" ; I think in
c most things I hard-coded "101" as _O_I1_ (the value of
c "o" when "i = 1", since the LSF is +/- 10 pixels and
c we want to be able to go at least one LSF width to the
c left of pixel i=1 and to the right of pixel i=_IDIM_
c
#define _O_I1_    101
#define _ODIM_  40401
 
 
 
 
 
c------------------------------------------------------------------------------
c
c the size of the supersampled LSF
c
#define _JCEN_    101
#define _JDIM_    201
c
c------------------------------------------------------------------------------
c
c number of pivot points in the underlying LSF formulation; this
c can't really be changed on the fly, since much of the influence-function
c stuff is solved for with this number in mine
c
#define _QDIM_     13
c
c------------------------------------------------------------------------------
c
c the max number of independent samplings of something that will get
c fed into a least-squares process
c
#define _LDIM_   99999
c
c------------------------------------------------------------------------------
 
#define _FREJ_SIGCLIP_ 4.0
 
#define _MDIM_ 16
 
 
c------------------------------------------------------------------------------
c
c  _SDIM_ : the total number of template nodes ---> 45000
c _LWSTEP_: the node spacing in velocity/c
c   _WWO_ : the velocity at s = 1
c
#define _SDIM_   30000
#define _LWSTEP_ 6.0d-6
#define _WW0_    5000.0d0
 
 
c------------------------------------------------------------------------------
c
c _SFIT_ is the maximum number of local samples to extract the
c template ; this has to be hardcoded since we will invert
c a matrix of this size in constructing the template ; we'll
c put together the template one little chunk at a time; the edges
c of this are not as accurate.
c
#define _SFIT_     51
 
c------------------------------------------------------------------------------
c
c this should just be WAVN... why do I have two variables for
c the same thing?  I should probably change WAVN, since it really
c doesn't control the WLS, except to the extent that it controls
c the LSF samplings; this controls the LSF sampling; KBEG and KEND
c will be the same as WAVB and WAVZ
c
#define _KDIM_ 09
#define _I4KA_ 0021
#define _I4KZ_ 4001
 
#define _WAVA_ 0001
#define _WAVB_ 0021
#define _WAVN_   09
#define _WAVY_ 4000
#define _WAVZ_ 4020
 
 
c------------------------------------------------------------------------------
c
c max number of template points for each order ; this is hardcoded
c because we may have to invert such a matrix; at present, I do it
c in chunks, but there's no real reason not to do to the whole matrix.
c
#define _NTEMMX_ 1001
 
 
 
c------------------------------------------------------------------------------
c
c wavelength ranger over which the quality will be assessed
c
#define _WMINQ_ 5000.0
#define _WMAXQ_ 6000.0
 
 
 
      program find_tem_COMPILER
      implicit none
 
 
 
c---------------------------------------------------------
c
c below are the arrays that contain the information
c relevant to each pixel:  its data value, the
c model value, its wavelength, the template value,
c the normalization that matches model to data,
c the average of its nearby neighbors (to allow
c looking at fractional residuals), and whether it
c is good or not.
c
      real*4,    dimension(:,:,:), allocatable :: pix_imn
      real*4,    dimension(:,:,:), allocatable :: mod_imn
      real*8,    dimension(:,:,:), allocatable :: wav_imn
      real*4,    dimension(:,:,:), allocatable :: tem_imn
      real*4,    dimension(:,:,:), allocatable :: rat_imn
      real*4,    dimension(:,:,:), allocatable :: bar_imn
      integer*1, dimension(:,:,:), allocatable :: rej_imn
      integer*1, dimension(:,:,:), allocatable :: tel_imn
      integer*1, dimension(:,:,:), allocatable :: uuu_imn
      integer*1, dimension(:,:,:), allocatable :: met_imn

c
c this is the parametrized LSF; q=1-->13 for nodes
c
      real*4, dimension(:,:,:,:), allocatable ::   lsfqkmn
 
c----------------------------------------------
c
c only need to hold for the current image; this contains
c the gsn map for each pixel
c
      real   gsn_im(_IDIM_,_MDIM_)
 
      logical unim(_NDIM_) ! whether this image will be used
      logical all_ok
 
      real    lsfj(_JDIM_),lsfju           ! local fine-sam LSF
      real    lsfjim(_JDIM_,_IDIM_,_MDIM_) ! fine-sam LSF for this exposure
 
 
c----------------------------------------------
c
c here are all the loopy type variables
c
      integer i, ii         ! contols pixel index
      real    ri, fi
      integer n, nimo, Ns   ! controls image number
      integer m             ! controls order number1
      integer k             ! controls sub-order sector
      integer q             ! controls lsf node
      integer j, jj         ! controls fine-sam lsf index
      integer o, oo         ! controls fine-sam model index
      integer s,  is, su    ! controls template index
      integer sa, sb        ! s-range to fit spectrum
      integer s1, s2        ! local s-range
 
 
c-----------------------------------------------
c
c
c
      real    tem0(_SDIM_)
      real    tems(_SDIM_)   ! current template
      real    temz(_SDIM_)   ! current, scaled/shifted to solar
      integer temu(_SDIM_)
      real    temw           ! function to eval template
 
      integer L, Ls
      integer sal(_LDIM_)
      integer sbl(_LDIM_)
 
      real*8  w, ww
      real*8  r8s_ww
      real*8  ww_r8s
      real*8  ww_s
      real    sunw, ftsw
      real    ftso(_ODIM_)
      real    fxto(_ODIM_)
      real    bari(_IDIM_)
      real    modi(_IDIM_)
      real    mod0(_IDIM_)
      real    mod1(_IDIM_)
      real    rati(_IDIM_)
      real    sigi(_IDIM_)
 
      real*8  AA,BB,CC,DD,EE
 
      real gsig(25), gsigi(4021)
      common /common_gsig/gsig
      data gsig/ 0.40, 0.50, 0.60, 0.70, 0.80,
     .           0.90, 1.00, 1.10, 1.20, 1.30,
     .           1.40, 1.50, 1.60, 1.70, 1.80,
     .           1.90, 2.00, 2.10, 2.20, 2.30,
     .           2.40, 2.50, 2.65, 2.80, 3.00/
 
      integer gopt, dopt, sopt
 
      integer NNIT, NTEM
      integer  NIT
      logical  done
 
 
c--------------------------------------------
c
c just a summary of fits that were made
c
      real    coeffw(5), scoeffw
      real    coeffwx(5)
      real    coeffg(_RATN_)
 
      integer p
      real*8  psum, pbar, ssum, sbar, xsum
      integer find_locrej, nrej
 
 
c--------------------------------------------
c
c variables to store the velocity information
c
      real    valu
      real    velu, find_vel4sun  ! this is actually z...
      real*8  zn(_NDIM_)
      real*8  z0(_NDIM_)

      real*8,  dimension(:,:), allocatable :: vln
      real*8,  dimension(:,:), allocatable :: tln
      integer, dimension(:,:), allocatable :: iln
      integer, dimension(:,:), allocatable :: mln
  
      real*8  zbarn(_NDIM_), zsign(_NDIM_) ! averages from above, by image
      real*8  zbarl(_LDIM_), zsigl(_LDIM_) ! averages from above, by line
 
      integer iargc
      character*80 STRING
 
      character*80 RJSOURCE

      real qqq
      real   ftsi(1 400 001)
      common /ftsw_f__/ftsi

      logical SPEC_FIX
      data SPEC_FIX/.false./

      integer ss
      real    fudge(120001,16)
      real    cfdge
      real    rfdge
      integer ifdge

      integer MORD_TO_DO
      integer MMIN, MMAX
      data MMIN/0/
      data MMAX/0/

      character*59  LABEL59
      character*80  FILENAME
      character*255 FILENAME255
      character*255 REFDIR
      common /REFDIR_/REFDIR

      logical FIRST_FUDGE
      data FIRST_FUDGE/.true./

      logical DOFUDGE
      common /DOFUDGE_/DOFUDGE
      data DOFUDGE/_DOFUDGE_/

      logical DOVORB
      common /DOVORB_/DOVORB
      data DOVORB/_DOVORB_/

      integer OUTPUT_MOD
      common /OUTPUT_MOD_/OUTPUT_MOD
      data OUTPUT_MOD/_OUTPUT_MOD_/
     
      integer OUTPUT_LSF
      common /OUTPUT_LSF_/OUTPUT_LSF
      data OUTPUT_LSF/_OUTPUT_LSF_/
     
      integer OUTPUT_VEL
      common /OUTPUT_VEL_/OUTPUT_VEL
      data OUTPUT_VEL/_OUTPUT_VEL_/
     
      integer OUTPUT_TEM
      common /OUTPUT_TEM_/OUTPUT_TEM
      data OUTPUT_TEM/_OUTPUT_TEM_/

      integer OUTPUT_NFO
      common /OUTPUT_NFO_/OUTPUT_NFO
      data OUTPUT_NFO/_OUTPUT_NFO_/
     
      integer OUTPUT_NRM
      common /OUTPUT_NRM_/OUTPUT_NRM
      data OUTPUT_NRM/_OUTPUT_NRM_/
     
      integer OUTPUT_LIN
      common /OUTPUT_LIN_/OUTPUT_LIN
      data OUTPUT_LIN/_OUTPUT_LIN_/
   
      integer NITo 
      data NITo/_NITo_/

      integer NITf
      data NITf/_NITf_/

      integer NARG, NARGs
      character*80 LOGFILE



c------------------------------------------------------------------------------
c
      call getenv("GRAND_REFDIR", REFDIR)

      RJSOURCE = 'NULL'
      LOGFILE  = '<STDOUT>'
      LABEL59  = '<NOLABEL>'

      NARGs = iargc()

      NIT = 0

      if (iargc().lt.3) then
 6001    format(a)
 6002    format(t14,a)
 6003    format(2x,a,t25,a)
 6004    format(2x,a,t25,a,i1,a)
 6005    format(2x,a,t25,a,i2,a)
 6006    format(2x,a,t25,a,L1,a)
 6007    format(t6,a)
         write(6,6001) 'usage: grand obslist label m'
         write(6,6002) '[nito=] [nitf=] [out=] [fudge{-,+}] [vorb+]'
         write(6,6002) '[{0,1,2}{0,1,2}{0,1,2}{0,1,2}{0,1,2}{0,1,2}]'
         write(6,6002) '[vel{0,1,2}] [tem{0,1,2}] [mod{0,1,2}]'
         write(6,6002) '[lsf{0,1,2}] [nrm{0,1,2}] [lin{0,1,2}]'
         write(6,6001) ''
         write(6,6001) 'Simultaneously fit a set of observed' //
     .     ' spectra obtained with an iodine cell.'
         write(6,6001) 'Return radial velocities, deconvolved' //
     .     ' stellar template, wavelengths,'
         write(6,6001) 'model fit, line spread functions, and' //
     .     ' normalization vectors.'
         write(6,6001) ''
         write(6,6001) 'positional arguments:'
         write(6,6003) 'obslist',
     .     'file containing list of observed spectra'
         write(6,6003) 'label',
     .     'root to use when generating output files'
         write(6,6003) 'm', 'order number (starting with 1, not 0)'
         write(6,6001) ''
         write(6,6001) 'optional arguments:'
         write(6,6003) '{0,1,2} x 6',
     .     'output flags for {vel}{tem}{mod}{lsf}{nrm}{lin}'
         write(6,6004) 'nito=',
     .     'starting iteration number            (default: ',
     .     _NITo_, ')'
         write(6,6005) 'nitf=',
     .     'final iteration number               (default: ',
     .     _NITf_, ')'
         write(6,6003) 'out=',
     .     'file to contain diagnostics          (default: stdout)'
         write(6,6006) 'vorb{-,+}',
     .     'read vorb from column 5 of obslist   (default: ',
     .     _DOVORB_, ')'
         write(6,6006) 'fudge{-,+}',
     .     'apply fudge correction to model      (default: ',
     .     _DOFUDGE_, ')'
         write(6,6004) 'vel{0,1,2}',
     .     'output flag for velocities           (default: ',
     .     _OUTPUT_VEL_, ')'
         write(6,6004) 'tem{0,1,2}',
     .     'output flag for template spectrum    (default: ',
     .     _OUTPUT_TEM_, ')'
         write(6,6004) 'mod{0,1,2}',
     .     'output flag for model fit            (default: ',
     .     _OUTPUT_MOD_, ')'
         write(6,6004) 'lsf{0,1,2}',
     .     'output flag for line spread function (default: ',
     .     _OUTPUT_LSF_, ')'
         write(6,6004) 'nrm{0,1,2}',
     .     'output flag for normalization vector (default: ',
     .     _OUTPUT_NRM_, ')'
         write(6,6004) 'lin{0,1,2}',
     .     'output flag for line regions         (default: ',
     .     _OUTPUT_LIN_, ')'
         write(6,6001) ''
         write(6,6001) 'notes:'
         write(6,6001) '  Output flags values:'
         write(6,6007) '0: do not write output files'
         write(6,6007) '1: write output files for current iteration'
         write(6,6007) '2: write output files for every iteration'
         write(6,6001) '  Maximum length of label is 59 characters'
         write(6,6001) '  Interpretation of relative order number' //
     .     ' m for HIRES:'
         write(6,6007) ' 1: order 71 = 4977-5065 Angstroms'
         write(6,6007) '16: order 56 = 6310-6368 Angstroms'
         write(6,6001) '  Output flags can be specified two' //
     .     ' different ways:'
         write(6,6007) 'collectively as a six-digit integer' //
     .     ' (default: 210100)'
         write(6,6007) 'individually as a three-letter tag and' //
     .     ' a one-digit integer (e.g., mod1)'
         stop
         endif

      call getarg(1,RJSOURCE)

      print*,'---> sense_raw...' 
      call sense_raw(RJSOURCE,Ns)
      print*,'---> Ns: ',Ns

      print*,'ALLOCATE MAIN VARIABLES...'
      allocate(pix_imn(_IDIM_,_MDIM_,Ns))
      allocate(mod_imn(_IDIM_,_MDIM_,Ns))
      allocate(wav_imn(_IDIM_,_MDIM_,Ns))
      allocate(tem_imn(_IDIM_,_MDIM_,Ns))
      allocate(rat_imn(_IDIM_,_MDIM_,Ns))
      allocate(bar_imn(_IDIM_,_MDIM_,Ns))
      allocate(rej_imn(_IDIM_,_MDIM_,Ns))
      allocate(tel_imn(_IDIM_,_MDIM_,Ns))
      allocate(uuu_imn(_IDIM_,_MDIM_,Ns))
      allocate(met_imn(_IDIM_,_MDIM_,Ns))
      allocate(lsfqkmn(_QDIM_,_KDIM_,_MDIM_,Ns))

      allocate(vln(_LDIM_,Ns)) ! line by line velocities
      allocate(tln(_LDIM_,Ns)) ! line by line depth
      allocate(iln(_LDIM_,Ns)) ! line by line pixel index
      allocate(mln(_LDIM_,Ns)) ! line by line order number


      call getarg(2,LABEL59)
      call getarg(3,STRING)

      read(STRING,*) MMIN
      MMAX = MMIN

      write(FILENAME,'(2a,i2.2,a)') trim(LABEL59),'.',mmin,'.par'
      open(60,file=FILENAME,status='unknown')
      write(60,'(a)') '# REQUIRED ARGUMENTS'
      write(60,'(2a)') 'obslist_file: ', trim(RJSOURCE)
      write(60,'(2a)') 'label: ', trim(LABEL59)
      write(60,'(a,i0)') 'order: ', mmin

      i = len_trim(LABEL59)
      write(LABEL59(i+1:i+3),'(a,i2.2)') '.', mmin

      do NARG = 4, NARGs
         call getarg(NARG,STRING)
         if (STRING(1:5).eq.'nito='.or.
     .       STRING(1:5).eq.'NITo=') then
            read(STRING(6:80),*) NITo
            NIT = NITo
            goto 77
            endif
         if (STRING(1:5).eq.'nitf='.or.
     .       STRING(1:5).eq.'NITf=') then
            read(STRING(6:80),*) NITf
            goto 77
            endif
         if (STRING(1:4).eq.'out='.or.
     .       STRING(1:4).eq.'OUT='.or.
     .       STRING(1:8).eq.'LOGFILE=') then
             read(STRING(index(STRING,'=')+1:80),*) LOGFILE
             open(6,file=LOGFILE,status='unknown')
             goto 77
             endif
         if (STRING(1:6).eq.'FUDGE+'.or.
     .       STRING(1:6).eq.'fudge+') then
             DOFUDGE = .true.
             goto 77
             endif
         if (STRING(1:6).eq.'FUDGE-'.or.
     .       STRING(1:6).eq.'fudge-') then
             DOFUDGE = .false.
             goto 77
             endif
         if (STRING(1:3).eq.'vel'.or.
     .       STRING(1:3).eq.'VEL') then
             read(STRING(4:4),*) OUTPUT_VEL
             goto 77
             endif
         if (STRING(1:3).eq.'tem'.or.
     .       STRING(1:3).eq.'TEM') then
             read(STRING(4:4),*) OUTPUT_TEM
             goto 77
             endif
         if (STRING(1:3).eq.'mod'.or.
     .       STRING(1:3).eq.'MOD') then
             read(STRING(4:4),*) OUTPUT_MOD
             goto 77
             endif
         if (STRING(1:3).eq.'lsf'.or.
     .       STRING(1:3).eq.'LSF') then
             read(STRING(4:4),*) OUTPUT_LSF
             goto 77
             endif
         if (STRING(1:3).eq.'nrm'.or.
     .       STRING(1:3).eq.'NRM') then
             read(STRING(4:4),*) OUTPUT_NRM
             goto 77
             endif
         if (STRING(1:3).eq.'lin'.or.
     .       STRING(1:3).eq.'LIN') then
             read(STRING(4:4),*) OUTPUT_LIN
             goto 77
             endif
         if (STRING(1:1).eq.'0'.or.
     .       STRING(1:1).eq.'1'.or.
     .       STRING(1:1).eq.'2') then
             if (STRING(1:1).lt.'0'.or.string(1:1).gt.'2'.or.
     .           STRING(2:2).lt.'0'.or.string(2:2).gt.'2'.or.
     .           STRING(3:3).lt.'0'.or.string(3:3).gt.'2'.or.
     .           STRING(4:4).lt.'0'.or.string(4:4).gt.'2'.or.
     .           STRING(5:5).lt.'0'.or.string(5:5).gt.'2'.or.
     .           STRING(6:6).lt.'0'.or.string(6:6).gt.'2') then
                 write(6,*) 'INVALID NUMBER FIELD...'
                 write(6,*) '   ARG#: ',NARG
                 write(6,*) '       : ',STRING(1:6)
                 stop
                 endif
             read(STRING(1:1),*) OUTPUT_VEL
             read(STRING(2:2),*) OUTPUT_TEM
             read(STRING(3:3),*) OUTPUT_MOD
             read(STRING(4:4),*) OUTPUT_LSF
             read(STRING(5:5),*) OUTPUT_NRM
             read(STRING(6:6),*) OUTPUT_LIN
             goto 77
             endif
         if (STRING(1:5).eq.'VORB-'.or.
     .       STRING(1:5).eq.'vorb-') then
             DOVORB = .false.
             goto 77
             endif
         if (STRING(1:5).eq.'VORB+'.or.
     .       STRING(1:5).eq.'vorb+') then
             DOVORB = .true.
             goto 77
             endif

         write(*,'(''COMMAND-LINE ARGUMENT UNPARSABLE: '',80a)') STRING
  77     continue
         enddo

      write(60,'(/,a)') '# OPTIONAL ARGUMENTS'
      write(60,'(2a)') 'logfile: ', trim(LOGFILE)
      write(60,'(a,i0)') 'nito: ', NITo
      write(60,'(a,i0)') 'nitf: ', NITf
      write(60,'(a,L1)') 'fudge: ', DOFUDGE
      write(60,'(a,L1)') 'vorb: ', DOVORB
      write(60,'(a,i1)') 'output_vel: ', OUTPUT_VEL
      write(60,'(a,i1)') 'output_tem: ', OUTPUT_TEM
      write(60,'(a,i1)') 'output_mod: ', OUTPUT_MOD
      write(60,'(a,i1)') 'output_lsf: ', OUTPUT_LSF
      write(60,'(a,i1)') 'output_nrm: ', OUTPUT_NRM
      write(60,'(a,i1)') 'output_lin: ', OUTPUT_LIN
      call flush(60)

      write(6,'(''#  '')') 
      write(6,'(''# grand.F '')') 
      write(6,'(''#  '')') 
      do NARG = 0, NARGs
         call getarg(NARG,STRING)
         write(6,'(''#   ARG'',i2.2,'': '',80a)') NARG,STRING
         enddo
      write(6,'(''#                '')') 
      write(6,'(''#   COMPULSORY PARAMETERS ----'')')
      write(6,'(''#      OBSLIST:  '',80a)') RJSOURCE
      write(6,'(''#        LABEL:  '',9a)') trim(LABEL59)
      write(6,'(''#      ORDER M:  '',i2)') MMIN
      write(6,'(''#                '')') 
      write(6,'(''#   OPTIONAL PARAMETERS ----'')')
      write(6,'(''#      LOGFILE: '',a80)') LOGFILE
      write(6,31) 'NITo',_NITo_,NITo
      write(6,31) 'NITf',_NITf_,NITf
      write(6,32) _DOFUDGE_,DOFUDGE
      write(6,34) _DOVORB_,DOVORB
      write(6,'(''#      OUTPUT_   '')') 
      write(6,33) 'VEL', _OUTPUT_VEL_, OUTPUT_VEL
      write(6,33) 'TEM', _OUTPUT_TEM_, OUTPUT_TEM
      write(6,33) 'MOD', _OUTPUT_MOD_, OUTPUT_MOD
      write(6,33) 'LSF', _OUTPUT_LSF_, OUTPUT_LSF
      write(6,33) 'NRM', _OUTPUT_NRM_, OUTPUT_NRM
      write(6,33) 'LIN', _OUTPUT_LIN_, OUTPUT_LIN
   31 format('#         ',a4,  ': [',i2.2,'] ',i2.2)
   32 format(    '#      DOFUDGE:  [',l1,'] ',l1)
   33 format('#          ',a3, ':  [',i1,'] ',i1)
   34 format(    '#      DOVORB:  [',l1,'] ',l1)
      write(6,'(''#  '')') 
      write(6,'(''#  '')') 

 
c-------------------------------------------------------
c
c initial set-up
c
      do n = 1, _NDIM_
         unim(n) = .false.
         zn(n) = 0.0
         zbarn(n) = 0.000d0
         zsign(n) = 0.999d0
         do m = 1, _MDIM_
            do i = 1, _IDIM_
               tel_imn(i,m,n) = 0
               enddo
            enddo
         enddo
 
c-------------------------------------------------------
c
c
c     SPEC_FIX = .false.
c     if (NIT.lt.0) then
c        NITf = abs(NIT)
c        NIT = 0
c        SPEC_FIX = .true.
c        call read_tem(tems,99)
c        endif
 
      write(6,*) '      '
      write(6,*) ' NIT  : ',NIT
      write(6,*) ' LABEL: ',trim(LABEL59)
      write(6,*) '      '


 
 
c-------------------------------------------------------
c
c if we're not starting from scratch, then
c read in the data from the previous iteration
c
      if (NIT.ne.0) then
         write(6,*) 'Not from scratch...'
         call read_mod(pix_imn,mod_imn,wav_imn,
     .                 rat_imn,bar_imn,rej_imn,
     .                 unim,NIT-1,Ns,LABEL59)
         call read_vel(zn,NIT-1,Ns,z0,LABEL59)
         call read_tem(tems,NIT-1,LABEL59)
         call read_lsf(lsfqkmn,NIT-1,Ns,LABEL59)
         call find_lines(tems,Ls,sal,sbl,temu,LABEL59)
         call tell_rej(wav_imn,tel_imn,Ns)

         call smoo_tem(tems,NIT)

 
         write(6,*) 'IMAGES NOT USED: '
         s = 0
         do n = 1, Ns
            if (.not.unim(n)) then
               s = s + 1
               write(6,*) '----> ',n,s
               endif
            enddo
 
         write(6,*) 'remodel...  Ns: ' ,Ns
         call fill_mod_imn(mod_imn,wav_imn,lsfqkmn,tems,z0,unim,
     .                     Ns,MMIN,MMAX)
 
         write(6,*) 'find_uuu...'
         call find_uuuimn(pix_imn,mod_imn,rat_imn,
     .                    bar_imn,rej_imn,uuu_imn,met_imn,
     .                    wav_imn,z0,unim,Ns,MMIN,MMAX)

 
         goto 777 ! go and start mid-stream
         endif
 
c-------------------------------------------------------
c
c NIT=0... start from scratch...
c
      write(6,*) ' '

 
      call read_raw(RJSOURCE,
     .              pix_imn,rat_imn,tem_imn,
     .              rej_imn,met_imn,
     .              wav_imn,unim,z0,Ns,MMIN,MMAX)

      do n = 1, Ns
         zn(n) = z0(n)
         enddo

c-------------------------------------------------------
c
c write more information to parameter file
c
      write(60,'(/,a)') '# ENVIRONMENT VARIABLES'
      write(60,'(2a)') 'GRAND_REFDIR: ', trim(REFDIR)
      write(60,'(/,a)') '# REFERENCE FILES'
      call flush(60)
 
c-------------------------------------------------------
c
c initial iteration, find the rough model template + LSF
c that describes the data
c
      do n = 1, Ns
         if (.not.unim(n)) goto 445
         do m = MMIN, MMAX
 
         write(6,*) 'call find_initmod...'
 
         call find_initmod(pix_imn(1,m,n),
     .                     wav_imn(1,m,n),gopt,dopt,sopt,m)
 
         write(6,*)
         write(6,*) 'INITMOD: ',n,m,gopt,dopt,sopt
         write(6,*)
 
         do i = 0001, 4021
            gsn_im(i,m)    = gopt
            if (.not.SPEC_FIX) tem_imn(i,m,n) = 1.0
            enddo
 
c------------------------------------------------------------
c
c take the gaussians derived above and make a fine-scaled LSF
c
         call gsni2lsfji(gsn_im(1,m),lsfjim(1,1,m))
 
         call find_ftso(ftso,wav_imn(1,m,n))
         call find_modi(mod_imn(1,m,n),ftso,lsfjim(1,1,m)) ! find the rough model of the obsns
         call find_rat0(pix_imn(1,m,n),
     .                  mod_imn(1,m,n),
     .                  rat_imn(1,m,n))
 
 
 
c------------------------------------------------------------
c
c go through an iterative process wherein we successively
c improve the WLS, the LSF, and the template in a (hopefully)
c virtuous cycle.  In the end, we'll stick with this WLS,
c but will end up improving the LSF and template quite a bit
c via detailed least-squares; we reject data points that are
c too discordant.  This is done each time, so that if a
c point is removed, it can be redeemed.
c
         done = .false.
         do NNIT = 1, 15

            call tweak_wls(pix_imn(1,m,n),rej_imn(1,m,n),
     .                     wav_imn(1,m,n),rat_imn(1,m,n),
     .                     tem_imn(1,m,n),
     .                     lsfjim(1,1,m),coeffw)
            call polyfit_w(wav_imn(1,m,n),
     .                     wav_imn(1,m,n),AA,BB,CC,DD,EE)
            call find_ftso(ftso,wav_imn(1,m,n))
 
            if (NNIT.le.22) then
               call tweak_gsn(pix_imn(1,m,n),rej_imn(1,m,n),
     .                        wav_imn(1,m,n),rat_imn(1,m,n),
     .                        gsn_im(1,m),lsfjim(1,1,m),coeffg)
               call lsfji2lsfqk(lsfjim(1,1,m),lsfqkmn(1,1,m,n))
               endif
 
            call find_modi(mod_imn(1,m,n),ftso,lsfjim(1,1,m)) ! find the rough model of the obsns
 
            call find_tem_snglQ(pix_imn(1,m,n),rej_imn(1,m,n),
     .                          mod_imn(1,m,n),rat_imn(1,m,n),
     .                          tem_imn(1,m,n))
 
            nrej = find_locrej(pix_imn(1,m,n),mod_imn(1,m,n),
     .                         rej_imn(1,m,n),sigi)
 
            write(6,444) 0,n,m,NNIT,nrej,(coeffw(p),p=1,5),
     .                                   (lsfqkmn(7,k,m,n),k=1,_KDIM_)
 
  444       format('NIT',i2.2,':',i3.3,'.',i2.2,'_',i2.2,3x,i4,1x,
     .                         3x,05(f6.3,1x),
     .                         5x,99(f7.5,1x))
            scoeffw = 0.
            do p = 1, 5
               scoeffw = scoeffw + abs(coeffw(p))
               enddo
            if (scoeffw.lt.0.025.and.done) goto 443
            if (NNIT.ge.3..and.scoeffw.lt.0.025) done = .true.
            enddo!NNIT
 443        continue
 
         call find_bari(pix_imn(1,m,n),bar_imn(1,m,n),     ! average pixel value, for norm'n
     .                                 rej_imn(1,m,n))
 
         enddo!m
 
 
c
c this routine will compare the current spectrm to the sun,
c line-center by line-center, and will shift report the
c apparent velocity. This is probably only good to 1000 m/s,
c but it will be good enough to allow us to start the iterative
c procedure.
c
 
         write(6,*)
         write(6,'(''---> z0(N): '',f12.9)') z0(N)
         write(6,*)


         if (DOFUDGE) then
            if (FIRST_FUDGE) then
               do s = 000001, 120001
                  do m = 01, 16
                     fudge(s,m) = 0.000
                     enddo
                  enddo
               FILENAME255 = trim(REFDIR) // 
     .                       '/BSTAR_CORRECTION.2012.10.23'
               open(19,file=FILENAME255,status='old')
               write(60,'(2a)'), 'fudge_file: ',trim(FILENAME255)
               do s = 000001, 120001
                  read(19,*) ss,rfdge,(fudge(s,m),m=1,16)
                  enddo
               w = 5450.00
               s = 0001 + (w-5000)*100
               write(6,*) '  ---> w: ',w
               write(6,*) '  ---> s: ',s
               write(6,*) '  ---> 6: ',fudge(s,6)
               write(6,*) '  ---> 7: ',fudge(s,7)
               write(6,*) '  ---> 8: ',fudge(s,8)
               FIRST_FUDGE = .false.
               close(19)
               endif
            do m = MMIN, MMAX
               do i = 0001, _IDIM_
                  rfdge = 01 + (wav_imn(i,m,n)-5000.0)*100.0
                  ifdge = int(rfdge)
                  cfdge = 0.00
                  if (ifdge.ge.01.and.ifdge.le.120001) then
                     cfdge = fudge(ifdge,m) + (rfdge-ifdge)
     .                                       *(fudge(ifdge+1,m)
     .                                        -fudge(ifdge  ,m))
                     endif
                  if (i.eq.2000) 
     .            write(6,*)'-> i=2000  PRE-FUDG: ',pix_imn(i,m,n),cfdge
                  pix_imn(i,m,n) = pix_imn(i,m,n) - bar_imn(i,m,n)*cfdge
                  if (i.eq.2000) 
     .            write(6,*)'-> i=2000 POST-FUDG: ',pix_imn(i,m,n),cfdge
                  enddo
               enddo
            endif
 
  445    continue! not unim
         enddo!n

c--------------------------------------------------------
c
c combine all the velocity-shifted template observations
c to get our initial estimate of the "average" template.
c
      if (.not.SPEC_FIX) then
         call make_init_tem(tems,tem_imn,wav_imn,z0,Ns,MMIN,MMAX)
c        call cnorm_tem(tems,tems,NIT,LABEL59)
         call cnorm_tem_bar(tems,tems,NIT,LABEL59,wav_imn(0501,MMIN,1),
     .                                            wav_imn(3501,MMIN,1))
         endif

      
      call tell_rej(wav_imn,tel_imn,Ns) ! now that we have a wls, flag telluric lines
      call fill_mod_imn(mod_imn,wav_imn,lsfqkmn,tems,z0,unim,
     .                  Ns,MMIN,MMAX)
 
 
 
 
c-----------------------------------------------------------
c
c output the current data to a save-file, so that it
c can either be inepseced or so that the program can be
c restarted at an intermediate stage
c
 
      if (NIT.eq.0) then
         call save_mod(pix_imn,mod_imn,wav_imn,
     .                 rat_imn,bar_imn,
     .                 rej_imn,tel_imn,met_imn,
     .                 z0,unim,NIT,Ns,LABEL59,MMIN,MMAX)
         call save_vel(zn,NIT,Ns,z0,zbarn,zsign,unim,LABEL59)
         call save_tem(tems,NIT,tems,LABEL59)
         call save_lsf(lsfqkmn,NIT,Ns,unim,LABEL59,MMIN,MMAX)
         call info_nim(pix_imn,mod_imn,wav_imn,
     .                 rat_imn,bar_imn,rej_imn,
     .                 unim,NIT,Ns,LABEL59,MMIN,MMAX)
         endif
 
      call find_lines(tems,Ls,sal,sbl,temu,LABEL59)
 
 
c
c ok... now that we've completed the first iteration, let's
c move on and iterate again, this time using least-squares
c to determine (a) a more flexible LSF and (b) a better
c sampled template..
c
 
      NIT = 1
 
  777 continue
 
      write(6,'('' '')')
      write(6,'(''*********************************************'')')
      write(6,'(''*****                    '')')
      write(6,'(''*****   START ITERATION: '',i4,3x,i4)') NIT
      write(6,'(''*****                    '')')
      write(6,'(''*********************************************'')')
      write(6,'('' '')')
 
      do s = 0001, _SDIM_
         tem0(s) = tems(s)
         enddo
 
      do n = 1, Ns
         do m = MMIN, MMAX
            if (unim(n)) then
               write(6,'(''                                  '')')
               write(6,'(''**********************************'')')
               write(6,'(''*****  '')')
               write(6,'(''*****  '')')
               write(6,'(''*****  NIT = '',i3.3,1x,i3.3)')  NIT
               write(6,'(''*****    n = '',i3.3,2x,f12.8)') n,z0(n)
               write(6,'(''*****    m = '',i3.3)')          m
               write(6,'(''*****  '')')
               write(6,'(''*****  '')')
               write(6,'(''**********************************'')')
               write(6,'(''                                  '')')
 
               call find_fxto(fxto,wav_imn(1,m,n),z0(n),tems)      ! fxto is the super sam x10 resolution
               call lsfqk2lsfji(lsfqkmn(1,1,m,n),lsfjim(1,1,m))    ! make a fine-scale LSF
               call find_modi(mod_imn(1,m,n),fxto,lsfjim(1,1,m))   ! find the rough model of the obsns
               call find_rati_bar(pix_imn(1,m,n),mod_imn(1,m,n),   ! find the scale that takes
     .                            bar_imn(1,m,n),rat_imn(1,m,n),   !    model*rat ---> pix
     .                            rej_imn(1,m,n),met_imn(1,m,n),NIT)                  
               call find_lsfq_sp(pix_imn(1,m,n),                   ! using the WLS and normalization,
     .                           rej_imn(1,m,n),                   !    find the spatially var LSF
     .                           rat_imn(1,m,n),fxto,              !    that best fits the fxto model
     .                           lsfqkmn(1,1,m,n),lsfjim(1,1,m))   !    to the data
               unim(n) = unim(n).and.
     .                   all_ok(lsfqkmn(1,1,m,n),(_QDIM_)*(_KDIM_))
               endif
            enddo!m
         write(6,*) '----------> IMAGE OK? ',unim(n)
         enddo!n
 
      close(90)
 
c   7  continue
 
 
 
 
 
c------------------------------------------------------------
c
c this is the range over which we'll solve for the
c template, between template points sa and sb
c
 
 
      call find_uuuimn(pix_imn,mod_imn,rat_imn,
     .                 bar_imn,rej_imn,uuu_imn,met_imn,
     .                 wav_imn,z0,unim,Ns,MMIN,MMAX)
 
      sa = r8s_ww(wav_imn(_WAVB_,MMIN,1))
      sb = r8s_ww(wav_imn(_WAVY_,MMAX,1))

c     sa = r8s_ww(wav_imn(0001,MMIN,1))
c     sb = r8s_ww(wav_imn(4020,MMAX,1))

      if (.true.) then
         call fill_mod_imn(mod_imn,wav_imn,lsfqkmn,tems,z0,unim,
     .                     Ns,MMIN,MMAX)
         call find_twktem(tems,
     .                    pix_imn,mod_imn,
     .                    rat_imn,bar_imn,uuu_imn,
     .                    wav_imn,z0,unim,Ns,MMIN)
         call find_lsqtem(tems,sa,sb,
     .                    pix_imn,rat_imn,uuu_imn,
     .                    wav_imn,z0,lsfqkmn,Ns,MMIN,MMAX)
         call smoo_tem(tems,NIT)
         endif 
 
      if (NIT.le.5.and.(.not.SPEC_FIX)) then
c        call cnorm_tem(tems,tems,NIT,LABEL59)
         call cnorm_tem_bar(tems,tems,NIT,LABEL59,wav_imn(0501,MMIN,1),
     .                                            wav_imn(3501,MMIN,1))
         endif

 
    7 continue
 
c-------------------------------------------------------------
c
c compute a new version of the current model
c
      call fill_mod_imn(mod_imn,wav_imn,lsfqkmn,tems,z0,unim,
     .                  Ns,MMIN,MMAX)
      call fill_rat_imn(pix_imn,mod_imn,
     .                  bar_imn,rat_imn,
     .                  rej_imn,met_imn,unim,Ns,
     .                  NIT,MMIN,MMAX)
 
 
 
 
c-------------------------------------------------------------
c
c re-compute the best velocity; this is no longer used in
c the modeling; it isn't even a nuisance parameter.  It's
c just a diagnostic, at this point.
c
      do n = 1, Ns

         if (unim(n)) then 
            write(6,*) '            '
            write(6,*) '-------> n: ',n,z0(n)
            write(6,*) '            '
 
           do m = MMIN, MMAX
              call lsfqk2lsfji(lsfqkmn(1,1,m,n),lsfjim(1,1,m))
              enddo
 
            zn(n) = z0(n)

            write(6,*) '-------> TWEAK_VEL: ' 
            call tweak_vel(pix_imn(1,1,n),mod_imn(1,1,n),        ! this does a one-param search for
     .                     wav_imn(1,1,n),                       !    velocity, using the entire
     .                     rat_imn(1,1,n),rej_imn(1,1,n),        !    spectrum ; there should be a better
     .                     met_imn(1,1,n),
     .                     tems,zn(n),lsfjim,MMIN,MMAX)          !    way that focuses on best lines
 
            write(6,*) '-------> VLINE_FIT: ' 
            call vlinefit_LSQ(pix_imn(1,1,n),wav_imn(1,1,n),
     .                        rat_imn(1,1,n),rej_imn(1,1,n),
     .                        tems,z0(n),lsfjim,
     .                        Ls,sal,sbl,n,
     .                        vln(1,n),tln(1,n),iln(1,n),mln(1,n),
     .                        MMIN,MMAX)
            endif
 
         enddo

 
      write(6,*) '  call anal_vlines: '
      write(6,*) '                Ls: ',Ls
      write(6,*) '                Ns: ',Ns
      call anal_vlines(vln,tln,iln,mln,  ! careful analysis of the quality of each line;
     .                 unim,Ns,Ls,NIT,   ! compute the average for each exposure using the best
     .                 zbarn,zsign,      ! lines in the best way (still needs to be optimized)
     .                 zbarl,zsigl,LABEL59)
 
c-------------------------------------------------------------
c
c save the current model ; eventually, I'll prolly do this at
c the end of the iteration, but for now, I'd like the chance
c of saving things mid-step, so I can mebbe start mid-step
c
      call save_mod(pix_imn,mod_imn,wav_imn,
     .              rat_imn,bar_imn,
     .              rej_imn,tel_imn,met_imn,
     .              z0,unim,NIT,Ns,LABEL59,MMIN,MMAX)
      call save_vel(zn,NIT,Ns,z0,zbarn,zsign,unim,LABEL59)
      call save_tem(tems,NIT,tem0,LABEL59)
      call save_lsf(lsfqkmn,NIT,Ns,unim,LABEL59,MMIN,MMAX)
      call info_nim(pix_imn,mod_imn,wav_imn,
     .              rat_imn,bar_imn,rej_imn,
     .              unim,NIT,Ns,LABEL59,MMIN,MMAX)
 
      NIT = NIT + 1
      if (NIT.le.NITf) goto 777

c-------------------------------------------------------------
c
c write internal parameters to parameter file
c
 6011    format(2x,a,i1,t16,': ',a)
         write(6,6006) 'mod',OUTPUT_MOD,'model output',
     .                 '1=most recent interation'
      write(60,'(/,a)') '# OBSERVATION INTERNAL PARAMETERS'
      write(60,'(a,i0)') '_IDIM_: ', _IDIM_
      write(60,'(a,i0)') '_IBEG_: ', _IBEG_
      write(60,'(a,i0)') '_IEND_: ', _IEND_
      write(60,'(/,a)') '# ECHELLE WAVELENGTH INTERNAL PARAMETERS'
      write(60,'(a,i0)') '_KDIM_: ', _KDIM_
      write(60,'(a,i0)') '_WAVA_: ', _WAVA_
      write(60,'(a,i0)') '_WAVB_: ', _WAVB_
      write(60,'(a,i0)') '_WAVY_: ', _WAVY_
      write(60,'(a,i0)') '_WAVZ_: ', _WAVZ_
      write(60,'(a,i0)') '_WAVN_: ', _WAVN_
      write(60,'(/,a)') '# TEMPLATE INTERNAL PARAMETERS'
      write(60,'(a,i0)') '_SDIM_: ', _SDIM_
      write(60,'(a,i0)') '_SFIT_: ', _SFIT_
      write(60,'(a,i0)') '_NTEMMX_: ', _NTEMMX_
      write(60,'(a,f0.2)') '_WW0_: ', _WW0_
      write(60,'(a,f0.8)') '_LWSTEP_: ', _LWSTEP_
      write(60,'(/,a)') '# LSF INTERNAL PARAMETERS'
      write(60,'(a,i0)') '_JCEN_: ', _JCEN_
      write(60,'(a,i0)') '_JDIM_: ', _JDIM_
      write(60,'(a,i0)') '_QDIM_: ', _QDIM_
      write(60,'(/,a)') '# CONTINUUM NORMALIZATION INTERNAL PARAMETERS'
      write(60,'(a,i0)') '_RATA_: ', _RATA_
      write(60,'(a,i0)') '_RATB_: ', _RATB_
      write(60,'(a,i0)') '_RATN_: ', _RATN_
      write(60,'(a,i0)') '_RATY_: ', _RATY_
      write(60,'(a,i0)') '_RATZ_: ', _RATZ_
      write(60,'(/,a)') '# FITTING INTERNAL PARAMETERS'
      write(60,'(a,f0.2)') '_FREJ_SIGCLIP_: ', _FREJ_SIGCLIP_
      write(60,'(a,f0.1)') '_WMINQ_: ', _WMINQ_
      write(60,'(a,f0.1)') '_WMAXQ_: ', _WMAXQ_
      write(60,'(/,a)') '# MISC INTERNAL PARAMETERS'
      write(60,'(a,L1)') '_LINUX_: ', _LINUX_
      close(60)
 
      stop 'stuck the landing!'
      end
 
 
c------------------------------------------------------------
c
c include files
c
 
c
c *****************************************************
c **                                                 **
c ** dnvertNxN -- subroutine                         **
c **                                                 **
c ** perform simple row operations on a to convert   **
c ** it into the identity matrix; perform the exact  **
c ** same operations on b                            **
c **                                                 **
c *****************************************************
c
 
      subroutine dinvertNxN(a,b,N,NDIM)
      implicit none
 
      integer NDIM        ! the dimensions of them both
      real*8 a(NDIM,NDIM) ! the input matrix; this will NOT be destroyed
      real*8 b(NDIM,NDIM) ! the output matrix;
      integer N           ! the dimensions of them both
 
      integer NMAX         ! change this before compiling to increase
      parameter (NMAX=117) ! the max no of dimensions
      real*8 c(NMAX,NMAX)   ! the copy of a that IS destroyed
 
      integer i,j
      real*8 h
      integer ii, jj
 
 
      if (N.gt.NMAX) then
         write(6,*) ' '
         write(6,*) 'routine: invertNxN '
         write(6,*) ' '
         write(6,*) '   must bomb because N was input:',N
         write(6,*) '   which  is  greater  than NMAX:',NMAX
         write(6,*) ' '
         write(6,*) ' '
         stop
         endif
 
 
c
c set C = A,
c and B = I (the identity matrix)
c
c
      do i = 1, N
      do j = 1, N
         c(i,j) = a(i,j)
         b(i,j) = 0.0
         enddo
         enddo
      do i = 1, N
         b(i,i) = 1.0
         enddo
 
 
c
c stage  I: reduce the lower left half of C to
c           zeros, and the diagonal to ones
c
      do j = 1, N
         h = c(j,j)
         do i = 1, N
            c(i,j) = c(i,j)/h
            b(i,j) = b(i,j)/h
            enddo
         do jj = j+1, N
            h = c(j,jj)
            do i = 1, N
               c(i,jj) = c(i,jj) - h*c(i,j)
               b(i,jj) = b(i,jj) - h*b(i,j)
               enddo
            enddo
         enddo
 
 
c
c stage II: reduce the upper right half of C to
c           zeros
c
      do j = 1, N-1
         do jj = j+1,N
            h = c(jj,j)
            do ii = 1, N
               c(ii,j) = c(ii,j) - h*c(ii,jj)
               b(ii,j) = b(ii,j) - h*b(ii,jj)
               enddo
            enddo
         enddo
 
c
c now C is the identity matrix...
c and B is the inverse of A...
c and A is still A...
c
      return
      end
 
 
c
c *****************************************************
c **                                                 **
c ** dnvertNxN -- subroutine                         **
c **                                                 **
c ** perform simple row operations on a to convert   **
c ** it into the identity matrix; perform the exact  **
c ** same operations on b                            **
c **                                                 **
c *****************************************************
c
 
      subroutine dinvertNxN_0OK(a,b,N,NDIM)
      implicit none
 
      integer NDIM        ! the dimensions of them both
      real*8 a(NDIM,NDIM) ! the input matrix; this will NOT be destroyed
      real*8 b(NDIM,NDIM) ! the output matrix;
      integer N           ! the dimensions of them both
 
      integer NMAX         ! change this before compiling to increase
      parameter (NMAX=117) ! the max no of dimensions
      real*8 c(NMAX,NMAX)   ! the copy of a that IS destroyed
 
      integer i,j
      real*8 h
      integer ii, jj
 
      integer z(NMAX)
      real*8 zsum
 
      if (N.gt.NMAX) then
         write(6,*) ' '
         write(6,*) 'routine: invertNxN '
         write(6,*) ' '
         write(6,*) '   must bomb because N was input:',N
         write(6,*) '   which  is  greater  than NMAX:',NMAX
         write(6,*) ' '
         write(6,*) ' '
         stop
         endif
 
 
c
c set C = A,
c and B = I (the identity matrix)
c
c
      do i = 1, N
         z(i) = 1
         zsum = 0.
         do j = 1, N
            zsum = zsum + abs(a(i,j)) + abs(a(j,i))
            enddo
         if (zsum.eq.0) then
            z(i) = 0
            a(i,i) = 1.0
            endif
         enddo
 
      do i = 1, N
      do j = 1, N
         c(i,j) = a(i,j)
         b(i,j) = 0.0
         enddo
         enddo
      do i = 1, N
         b(i,i) = 1.0
         enddo
 
 
c
c stage  I: reduce the lower left half of C to
c           zeros, and the diagonal to ones
c
      do j = 1, N
         h = c(j,j)
         do i = 1, N
            c(i,j) = c(i,j)/h
            b(i,j) = b(i,j)/h
            enddo
         do jj = j+1, N
            h = c(j,jj)
            do i = 1, N
               c(i,jj) = c(i,jj) - h*c(i,j)
               b(i,jj) = b(i,jj) - h*b(i,j)
               enddo
            enddo
         enddo
 
 
c
c stage II: reduce the upper right half of C to
c           zeros
c
      do j = 1, N-1
         do jj = j+1,N
            h = c(jj,j)
            do ii = 1, N
               c(ii,j) = c(ii,j) - h*c(ii,jj)
               b(ii,j) = b(ii,j) - h*b(ii,jj)
               enddo
            enddo
         enddo
 
c
c now C is the identity matrix...
c and B is the inverse of A...
c and A is still A...
c
 
      do i = 1, N
      do j = 1, N
         if (z(i).eq.0) then
            b(i,j) = 0.
            a(i,j) = 0.
            endif
         enddo
         enddo
 
      return
      end
 
 
 
c-----------------------------------------------
c
c this will take the initial template and will find
c something akin to a pixel-by-pixel normalization
c for it; not perfect, but not bad ; the normalized
c template will then be temi/nrmi
c
c-----------------------------------------------
 
      subroutine find_rat0(pixi,modi,nrmi)
      implicit none
 
      real   pixi(_IDIM_)
      real   modi(_IDIM_)
      real   nrmi(_IDIM_)
 
      real    temi(_IDIM_)
      real*8  normal(5,5)
      real*8  lamron(5,5)
      real*8  coeffs(5)
      real*8  vector(5)
      real    projn_i(5,_IDIM_)
      logical use(_IDIM_)
 
      integer reju, rejn
      integer i, n
      integer ni,nj
      integer n1,n2
      real*8  rn
      real    ssum, ssig
      integer NIT
      integer nsum, nuse
 
      do i = 0001, _IDIM_
         temi(i) = pixi(i)/modi(i)
         enddo
 
c
c will evaluate this at 5 points, with
c a spline in between
c
      do i = 0001, _IDIM_
         rn = 1 + 4.*(i-_RATB_)/(_RATY_-_RATB_)
         call rspline_fp(rn,5,projn_i(1,i),n1,n2)
         use(i) = .true.
         enddo
 
c
c iterate 9x to remove things inconsistent
c with the "modal" value
c
      do NIT = 1, 9
 
 
c
c construct the normal equations using
c only the points that are near the
c continuum
c
         do ni = 1, 5
         do nj = 1, 5
            normal(ni,nj) = 0.
            do i = _IBEG_, _IEND_
               if (use(i))
     .         normal(ni,nj) = normal(ni,nj) + projn_i(ni,i)
     .                                        *projn_i(nj,i)
               enddo
            enddo
            enddo
 
         call dinvertNxN(normal,lamron,5,5)
 
         do n = 1, 5
            vector(n) = 0
            do i = _IBEG_, _IEND_
               if (use(i))
     .         vector(n) =
     .         vector(n) + temi(i)*projn_i(n,i)
               enddo
            enddo
 
c
c use the inverted matrix to get the value of
c each pivot-point node
c
         do n = 1, 5
            coeffs(n) = 0
            do nj = 1, 5
               coeffs(n) = coeffs(n) + lamron(n,nj)*vector(nj)
               enddo
            enddo
 
c
c use the new solution to determine the continuum
c value at each pixel
c
         do i = 0001, _IDIM_
            nrmi(i) = 0.
            do n = 1, 5
               nrmi(i) = nrmi(i) + coeffs(n)*projn_i(n,i)
               enddo
            enddo
 
c
c use abs value to define a sigma that isn't biased
c by the extreme outliers
c
         nsum = 0
         ssum = 0
         do i = 0001, _IDIM_
            if (use(i)) then
               ssum = ssum + abs((temi(i)-nrmi(i))/nrmi(i))
               nsum = nsum + 1
               endif
            enddo
         ssig = ssum/nsum
 
 
c
c deselect things that don't look like continuum ; these
c flags won't be used beyond this routine
c
         nuse = 0
         do i = 0001, _IDIM_
            if (abs((temi(i)-nrmi(i))/nrmi(i)).gt.2.5*ssig)
     .          use(i) = .false.
            if (use(i)) nuse = nuse + 1
            enddo
 
         enddo
 
      return ! we now have nrmi(i), which represents the continuum
      end
 
 
 
c-----------------------------------------------
c
c this will take the initial template and will find
c something akin to a pixel-by-pixel normalization
c for it; not perfect, but not bad ; the normalized
c template will then be temi/nrmi
c
c-----------------------------------------------
 
      subroutine find_bar0(temi,nrmi)
      implicit none
 
      real   temi(_IDIM_)
      real   nrmi(_IDIM_)
 
      real*8  normal(5,5)
      real*8  lamron(5,5)
      real*8  coeffs(5)
      real*8  vector(5)
      real    projn_i(5,_IDIM_)
      logical use(_IDIM_)
 
      integer reju, rejn
      integer i, n
      integer ni,nj
      integer n1,n2
      real*8  rn
      real    ssum, ssig
      integer NIT
      integer nsum, nuse
 
 
c
c will evaluate this at 5 points, with
c a spline in between
c
      do i = 0001, _IDIM_
         rn = 1 + 4.*(i-_RATB_)/(_RATY_-_RATB_)
         call rspline_fp(rn,5,projn_i(1,i),n1,n2)
         use(i) = .true.
         enddo
 
c
c iterate 9x to remove things inconsistent
c with the "modal" value
c
      do NIT = 1, 9
 
 
c
c construct the normal equations using
c only the points that are near the
c continuum
c
         do ni = 1, 5
         do nj = 1, 5
            normal(ni,nj) = 0.
            do i = _IBEG_, _IEND_
               if (use(i))
     .         normal(ni,nj) = normal(ni,nj) + projn_i(ni,i)
     .                                        *projn_i(nj,i)
               enddo
            enddo
            enddo
 
         call dinvertNxN(normal,lamron,5,5)
 
         do n = 1, 5
            vector(n) = 0
            do i = _IBEG_, _IEND_
               if (use(i))
     .         vector(n) =
     .         vector(n) + temi(i)*projn_i(n,i)
               enddo
            enddo
 
c
c use the inverted matrix to get the value of
c each pivot-point node
c
         do n = 1, 5
            coeffs(n) = 0
            do nj = 1, 5
               coeffs(n) = coeffs(n) + lamron(n,nj)*vector(nj)
               enddo
            enddo
 
c
c use the new solution to determine the continuum
c value at each pixel
c
         do i = 0001, _IDIM_
            nrmi(i) = 0.
            do n = 1, 5
               nrmi(i) = nrmi(i) + coeffs(n)*projn_i(n,i)
               enddo
            enddo
 
c
c use abs value to define a sigma that isn't biased
c by the extreme outliers
c
         nsum = 0
         ssum = 0
         do i = 0001, _IDIM_
            if (use(i)) then
               ssum = ssum + abs((temi(i)-nrmi(i))/nrmi(i))
               nsum = nsum + 1
               endif
            enddo
         ssig = ssum/nsum
 
 
c
c deselect things that don't look like continuum ; these
c flags won't be used beyond this routine
c
         nuse = 0
         do i = 0001, _IDIM_
            if (abs((temi(i)-nrmi(i))/nrmi(i)).gt.2.5*ssig)
     .          use(i) = .false.
            if (use(i)) nuse = nuse + 1
            enddo
 
         enddo
 
      return ! we now have nrmi(i), which represents the continuum
      end
 
 
 
c=======================================================================
c
c This routine will take in the pix and rej array and will compute a
c local average of the pixel values.  This is useful for constructing
c scaled residuals that aren't too biased by whether the particular
c pixel is high or low.  The compiler parameter _RATN_ tells you how
c many nodes will be computed.  Around the edges, it just extrapolates.
c
      subroutine find_bari(pixi,bari,reji)
      implicit none
 
      real   pixi(_IDIM_)   ! INPUT  pixel array
      real   bari(_IDIM_)   ! OUTPUT average array
      integer*1 reji(_IDIM_)   ! INPUT  array of pixels to use
 
      integer i, i1, i2
      real    psum
      integer nsum
      real barp(_RATN_)
      integer p
      real rp, fp
      integer  ip
      integer WIDTH
 
 
c-----------------------------------------------------
c
c compute the average pixel value about several nodes ; use
c overlapping bins (x2) ; don't worry about double-counting
c at the margins, this is only a diagnostic
c
      do p = 001, _RATN_
         if (p.ge.0002.and.p.le._RATN_-1) then
            i1 = max(_RATA_,_RATB_+(p-1)*(_RATY_-_RATB_)/(_RATN_-2))
            i2 = min(_RATZ_,_RATB_+(p+1)*(_RATY_-_RATB_)/(_RATN_-2))
c
c 2013-Jan-10, Valenti, The following code creates intervals that
c start at _RATB_ (21-101) and end at _RATY_ (3919-4000). The code
c above makes a first bin (61-142) and last bin (4000-4081) that
c start after _RATB_ and extend beyond _RATY_.
c
c           i1 = _RATB_+(p-2)*(_RATY_-_RATB_)/(_RATN_-1)
c           i2 = _RATB_+p*(_RATY_-_RATB_)/(_RATN_-1)
c
            goto 1 ! quicker, since most pix not near edge
            endif
         if (p.eq. 0001) then
            i1 = _RATA_
            i2 = _RATB_
            endif
         if (p.eq._RATN_) then
            i1 = _RATY_
            i2 = _RATZ_
            endif
    1    continue
         psum = 0.
         nsum = 0
         do i = i1, i2
            if (reji(i).ne.1) then    ! only reject if bad; ok if it's on edge
               psum = psum + pixi(i)
               nsum = nsum + 1
               endif
            enddo
         barp(p) = 99999.0
         if (nsum.gt.0) barp(p) = psum/nsum  ! error handling
         enddo
 
 
c-----------------------------------------------------
c
c now fill each pixel value with this average
c
      do i = 0001 , _IDIM_
         if (i.le._RATB_) then
            bari(i) = barp(1) + (barp(2)-barp(1))*(  i   -_RATA_)/
     .                                            (_RATB_-_RATA_)
            endif
         if (i.ge._RATY_) then
            bari(i) = barp(_RATN_-2) +
     .               (barp(_RATN_-1)-barp(_RATN_-2))*(  i   -_RATY_)/
     .                                               (_RATZ_-_RATY_)
            endif
         if (i.gt._RATB_.and.i.lt._RATY_) then
            rp = 2 + (_RATN_-3.)*(i-_RATB_)/(_RATY_-_RATB_)  ! the "3." is
            ip = int(rp)                                     ! important...
            fp = rp-ip
            bari(i) = barp(ip) + fp*(barp(ip+1)-barp(ip))
            endif
         enddo
 
      return
      end
 
 
c=======================================================================
c
c This routine uses least squares to determine the ratio between the
c observations pixi() and the model modi() along the order from pixels
c _RATZ_ to _RATA_.  It returns this in the array rati(), which scales
c things such that:   pix(i) ~ rati(i)*modi(i).
c
c It also initializes reji() and it selectively rejects the bad pixels;
c one at a time, but within regions, so it should be relatively fast.
c
c The variable "p" goes from 1 ---> NP to go from i = _RATA_ ----> _RATZ_,
c the full range of useful pixel values.  For each pixel, p1 and p2 will
c tell you which nodes are relevant, ie whether proj() has any power there.
c Proj() will tell you how each node point contributes to rat() at that
c location.
c
      subroutine find_rati_bar(pixi,modi,bari,rati,reji,meti,NIT)
      implicit none
 
      real       pixi(_IDIM_)        ! INPUT : the array of pixel values
      real       modi(_IDIM_)        ! INPUT : the unscaled model values
      real       bari(_IDIM_)
      real       rati(_IDIM_)        ! OUTPUT: the scaling for each pixel
      integer*1  reji(_IDIM_)        ! OUTPUT: the determined quality
      integer*1  meti(_IDIM_)        ! OUTPUT: the determined quality
      integer    NIT                 !         of each pixel
 
      integer i                      ! pixel number
      integer p
      real*8  barp(201), sigp(201)
      real*8  xl(999), bar, sig
      integer Ls, Lu 
      integer ip
      real*8  r
      real    rp, fp
      real    ptot, mtot, rtot, rat
      integer ntot
      real    rmax
      integer imax
      integer nrej, irej

      irej = 1

      do i = 0001, 4021
         reji(i) = 0
         if (meti(i).ge.2) reji(i) = 1
         enddo

      do p = 002, 200
         ip = 1 + (p-1)*20
    1    continue
         ptot = 0.
         mtot = 0.
c
c 2013-Jan-10, Valenti, The pixel intervals for each segment are
c centered at [1], 21, 41, ..., 3961, 3981, [4001] and span pixel
c interval 1-41, 21-61, ..., 3941-3981, 3961-4001. Note the overlap
c and the failure to reach pixel 4020, which is the value of _RATZ_.
c
         do i = max(0001,ip-20), min(4001,ip+20)
            if (reji(i).eq.0.and.meti(i).le.1) then
               ptot = ptot + pixi(i)
               mtot = mtot + modi(i)
               endif 
            enddo   
         rat = ptot/mtot
         ntot = 0
         rtot = 0.
         rmax = 0.
         imax = 0.
         do i = max(0001,ip-20), min(4001,ip+20)
            if (reji(i).eq.0.and.meti(i).le.1) then
               r = abs(pixi(i)-rat*modi(i))/bari(i)
               ntot = ntot + 1
               rtot = rtot + r
               if (r.gt.rmax) then
                  rmax = r
                  imax = i
                  endif
               endif 
            enddo
         sig = rtot/ntot
         if ((NIT.le.1.and.rmax.gt.6.0*sig).or.
     .       (NIT.eq.2.and.rmax.gt.5.0*sig).or.
     .       (NIT.ge.3.and.rmax.gt.4.0*sig)) then
            reji(imax) = 1
            goto 1
            endif
         barp(p) = rat
         barp(p) = rat
         sigp(p) = sig
         enddo

      barp(001) = barp(002) - (barp(003)-barp(002))
      barp(201) = barp(200) - (barp(199)-barp(200))
      sigp(001) = sigp(002)
      sigp(201) = sigp(200)

      nrej = 0
      do i = 0001, 4021
         rp = 1 + (i-1)/20.00
         ip = int(rp)
         if (ip.gt.200) ip = 200
         fp = rp-ip
c
c 2013-Jan-10, Valenti, Pixel 1 maps to rp=1.0, which is correct.
c Pixel 21 maps to rp=2.0, which is correct. Pixel 4001 maps to
c rp=201.0, which is correct. Because the last measured node is
c p=200 at pixel 3981, all pixels beyond 3981 (i.e. 3982-4020)
c are handled by linear extrapolation.
c
         bar     = barp(ip) + fp*(barp(ip+1)-barp(ip)) 
         sig     = sigp(ip) + fp*(sigp(ip+1)-sigp(ip)) 
         rati(i) = bar
         r = (pixi(i)-rati(i)*modi(i))/bari(i)
         if (reji(i).eq.1) then
            nrej = nrej + 1
            irej = i
            endif
         enddo

      write(*,'(''find_rati_bar:  '',i4,1x,i4.4,1x,i1)') 
     .            nrej,irej,reji(irej)
 
      return
      end
 
 
 
 
c==============================================================================
c
c This routine will go through the image pixels in the same bins we
c used to find_rat, and will look for outliers and will remove them.
c
      integer function find_reji(pixi,modr,reji)
      implicit none
 
      real   pixi(_IDIM_) ! INPUT  pixel values
      real   modr(_IDIM_) ! INPUT  scaled model
      integer*1 reji(_IDIM_) ! OUTPUT pixel quality array 0 = good
 
      real*8  ssum
      integer nsum, ntot
      real    a, amax, sig
      integer imax
      integer p, i
      integer i1, i2
      integer j1, j2
 
      find_reji = 0
 
 
 
c--------------------------------------------------
c
c go through the pixels between each pair of nodes
c to see if there are any that need to be rejected
c
      do p = 01, _RATN_-1
 
c
c define the i-range for each node pair ; the minus
c one on i2 ensures that we won't double count pixels
c
         i1 = _RATB_ + (p-1)*(_RATY_-_RATB_)/
     .                       (_RATN_-2)
         i2 = _RATB_ + (p  )*(_RATY_-_RATB_)/
     .                       (_RATN_-2)-1
c
c 2013-Jan-10, Valenti, The code above makes a first bin (61-101) and
c last bin (4000-4039) that start after _RATB_ and extend beyond _RATY_.
c The code below creates intervals that start at _RATB_ (21-60) and
c end at _RATY_ (3960-4000).
c
c           i1 = _RATB_+(p-2)*(_RATY_-_RATB_+1)/(_RATN_-2)
c           i2 = _RATB_+(p-1)*(_RATY_-_RATB_+1)/(_RATN_-2)-1
c
         if (p.eq.001) then
            i1 = _RATA_
            i2 = _RATB_-1
            endif
         if (p.eq._RATN_-1) then
            i1 = _RATY_
c
c 2013-Jan-10, Valenti, include both _RATB_ and _RATY_ or neither
c
c           i1 = _RATY_+1
c
            i2 = _RATZ_
            endif
 
c
c now, go through these pixels; compute the rms
c between the observations and the best-scaled
c model
c
         nsum = 0
         ssum = 0.0
         amax = 0.
         imax = 0
         do i = i1, i2
            if (reji(i).eq.0) then
               a = abs(pixi(i)-modr(i))
               nsum = nsum + 1
               ssum = ssum + a
               if (a.gt.amax) then
                  amax = a
                  imax = i
                  endif
               endif
            enddo
         sig = ssum/nsum
 
 
c
c see if any of the pixels between (i1:i2) are more
c than 3.5 sigma off, relative to their neighbors
c
         if (amax.gt.(_FREJ_SIGCLIP_)*sig) then
            reji(imax) = 1
            find_reji = find_reji + 1
            endif
 
         enddo  ! loop through all pairs of node points
 
 
      return ! report back how many were newly flagged bad
 
      end
 
 
 
c=======================================================================
c
c This function will return the spline value at location rx in the array
c    yp(1:Ns); the routine uses rspline_fp() which gives the projection
c    of the spline nodes on point at location rx
c
      real function rspline(rx,yp,Ns)
      implicit none
 
      integer Ns                 ! INPUT: the number of elements in y()
      real*8 rx                  ! INPUT: the location within the y array
      real   yp(Ns)              ! INPUT: the y array
 
      integer ic                 ! the smallest possible array location
      real fp6(6)                ! the maximum number of contributing elements
      integer ip                 !
      integer ip1, ip2           ! the range of relevant elements
 
      ic = int(rx)-2
      if (ic.lt.   1) ic = 1
      if (ic.gt.Ns-5) ic = Ns-5
 
c
c This routine will determine the how each of the nodes in the spline
c will contribute to the value of the spline at location rx.  At most,
c 6 can contribute, in that we fit a 5-point quartic to each pixel for
c its value/slope and we use these two quantities between two pixels.
c All the subtraction and array bounds below are designed to allow me
c to just use local variables here and not to have to generate a whole
c new array with dimension Ns.
c
      call rspline_fp(1+rx-ic,6,fp6(1),ip1,ip2)
 
c
c Now, just add up the contributions to the value of
c the spline at location rx.
c
      rspline = 0.
      do ip = ip1, ip2
         rspline = rspline + fp6(ip)*yp(ip+ic-1)
         enddo
 
      return
      end
 
 
c==============================================================================
c
c Subroutine rspline_fp() will return the projection of each spline
c node on the splined value at point rx, for a domain that goes from
c x = 1 to x = NP.
c
c The value of the spline of the function y(1:N) at location rx is just
c going to be:
c
c             ip2
c     y(rx) = SUM  y(ip)*fp(ip)
c           ip=ip1
c
c To evaluate the spline we must determine:
c (1) which two nodes we want to spine between (x1,x2)
c (2) what is the value of the function at the endpoints
c (3) what is the value of the function's derivative at the endpoints
c
c This provides 4 constraints on f(x) between x=x1 and x=x2, and
c we can therefore solve for a cubic function of the form:
c
c     f(x) = A + B*x + C*x**2 + D*x**3
c
      subroutine rspline_fp(rx,Ns,fp,ip1,ip2)
      implicit none
 
      real*8  rx          ! INPUT:  location at which spline is to be eval
      integer Ns          ! INPUT:  number of nodes in spline
      real    fp(Ns)      ! OUTPUT: contribution of each node to spline val
      integer ip1, ip2    ! OUTPUT: the range over which nodes contribute
 
      integer ix
      real    fx, fx1, fx2, fx3, fx4
      real A05QUAR(5) ! these are the lsq solutions
      real B05QUAR(5) !    for a quartic fit to five
      real C05QUAR(5) !    points:
      real D05QUAR(5) !
      real E05QUAR(5) ! f(x) = A + B*x + C*x**2 + D*x**3 + E*x**4
      integer i,ip
 
c----------------------------------------------------
c
c Once-for all determined values for the lsq fit of a
c quartic to 5 points.  If point 3 is in the middle,
c then:
c
c             A = SUM A05QUAR(i)*pix(i)
c             B = SUM A05QUAR(i)*pix(i) ....
c
c then:  f(x) = A + B*(x-3) + C*(x-3)**2 + ...
c
c We will use this two ways.  We will also use it to
c evaluate the spline within the inner/outer three
c gridpoints.
c
      data A05QUAR/ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000/
      data B05QUAR/ 0.083333,-0.666667, 0.000000, 0.666667,-0.083333/
      data C05QUAR/-0.041667, 0.666667,-1.250000, 0.666667,-0.041667/
      data D05QUAR/-0.083333, 0.166667, 0.000000,-0.166667, 0.083333/
      data E05QUAR/ 0.041667,-0.166667, 0.250000,-0.166667, 0.041667/
 
 
c---------------------------------------------------------------------
c
c The MAP array is the compact least-squares solution to what form
c the spline function should have between two points.  It comes from
c specifying the value of the function and of its derivatives at each
c node.  When we want to evaluate the function between two node points,
c we find the cubic polynomial that goes through the two points with
c the predetermined slopes.  The following "MAP" will accomplish this.
c If we want the value between 3 and 4, using Y(1:6), then summing
c i = 1 to 6, we get:
c
c              A = SUM  MAP(i,1)*Y(i)
c              B = SUM  MAP(i,2)*Y(i)
c              C = SUM  MAP(i,3)*Y(i)
c              D = SUM  MAP(i,4)*Y(i)
c
c Then, between 3 and 4, we'll see that:
c
c     y(x) = A + B*(x-3) + C*(x-3)**2 + D*(x-3)**3
c
c It really is that easy.  The great thing is this makes finding the
c best-fitting spline a PURELY linear process.
c
      real MAP(6,4)
      data MAP/ 0.00000, 0.00000, 1.00000, 0.00000, 0.00000, 0.00000,
     .          0.08333,-0.66667, 0.00000, 0.66667,-0.08333, 0.00000,
     .         -0.16666, 1.25001,-2.33333, 1.66666,-0.50001, 0.08333,
     .          0.08333,-0.58334, 1.33333,-1.33333, 0.58334,-0.08333/
 
 
c-------------------------------------------------------------
c
c This is the meaty part; the point is well within the bounds
c of the function, and we can use our full spline power
c
      if (rx.ge.03.and.rx.lt.Ns-2) then
         ix = int(rx)
         fx = rx-ix
         ip1 = ix-2
         ip2 = ix+3
         fx1 = fx
         fx2 = fx*fx
         fx3 = fx*fx*fx
         fx4 = fx*fx*fx*fx
         do ip = ip1, ip2                   ! as described above, the MAP
            fp(ip) = MAP(1+ip-ip1,1)*(1.0)  ! matrix tells how the values
     .             + MAP(1+ip-ip1,2)*(fx1)  ! of the function in the 6 pts
     .             + MAP(1+ip-ip1,3)*(fx2)  ! surrounding the location of
     .             + MAP(1+ip-ip1,4)*(fx3)  ! interest contrib to the spline
            enddo
         return
         endif
 
c-------------------------------------------------------------
c
c Below, the point is somewhere within the first three.  To
c determine the spline, we just find the best-fitting quartic
c for the first 5 points, and return the weights that will
c give us the value
c
      if (rx.le.01+2) then
         ix = 3
         fx = rx-ix
         fx1 = fx
         fx2 = fx*fx
         fx3 = fx*fx*fx
         fx4 = fx*fx*fx*fx
         ip1 = 1
         ip2 = 5
         do ip = ip1, ip2
            fp(ip) = (1.0)*A05QUAR(ip-ip1+1)
     .             + (fx1)*B05QUAR(ip-ip1+1)
     .             + (fx2)*C05QUAR(ip-ip1+1)
     .             + (fx3)*D05QUAR(ip-ip1+1)
     .             + (fx4)*E05QUAR(ip-ip1+1)
            enddo
         return
         endif
 
c-------------------------------------------------------------
c
c Same as above, but for the outer 3 points.
c
      if (rx.ge.Ns-2) then
         ix = Ns-2
         fx = rx-ix
         fx1 = fx
         fx2 = fx*fx
         fx3 = fx*fx*fx
         fx4 = fx*fx*fx*fx
         ip1 = Ns-4
         ip2 = Ns
         do ip = ip1, ip2
            fp(ip) = (1.0)*A05QUAR(ip-ip1+1)
     .             + (fx1)*B05QUAR(ip-ip1+1)
     .             + (fx2)*C05QUAR(ip-ip1+1)
     .             + (fx3)*D05QUAR(ip-ip1+1)
     .             + (fx4)*E05QUAR(ip-ip1+1)
            enddo
         return
         endif
 
      end
 
 
 
c--------------------------------------------------------------
c
c this is the workhorse routine that will take raw observations
c and will construct a template....
c
c
      subroutine find_loctem(loctem,s1,s2,
     .                       pix_imn,rat_imn,rej_imn,
     .                       wav_imn,zn,lsfqkmn,Ns,MMIN,MMAX)
      implicit none
 
      real loctem(_SFIT_)                            ! OUTPUT: temlate to be filled
      integer   s1, s2                               ! INPUT:  range to fit
      real      pix_imn(_IDIM_,_MDIM_,Ns)            ! INPUT:  observations
      real      rat_imn(_IDIM_,_MDIM_,Ns)            ! INPUT:  ratio between obs/mod
      integer*1 rej_imn(_IDIM_,_MDIM_,Ns)            ! INPUT:  reject this pixel?
      real*8    wav_imn(_IDIM_,_MDIM_,Ns)            ! INPUT:  iodine wavelength
      real*8    zn(Ns)
      real      lsfqkmn(_QDIM_,_KDIM_,_MDIM_,Ns)     ! INPUT:  LSF
      integer   Ns                                   ! INPUT:  number of images
      integer   MMIN, MMAX 
 
      real*8  ww1, ww2, ww_s, ww_rs
      real*8  rs_ww
      integer Ss
      real*8  r8s_ww
 
      integer sa, sb
      integer L, Ls
      real    pixl(_LDIM_)
      real*8  wavl(_LDIM_)
      real*8  wwvl(_LDIM_)
      real    ratl(_LDIM_)
      real    teml(_LDIM_)
      real    barl(_LDIM_)
      real    modl(_LDIM_)
      real    issl(_LDIM_)
 
      real    proj(_SFIT_)
      real    modls(_SFIT_,_LDIM_)
 
      real    ftsw, sunw
      real    ftso(_ODIM_)
      real    temo(_ODIM_)
      real    ftsv(_ODIM_)
      real    modo(_ODIM_)
 
      real*8  normal(_SFIT_,_SFIT_)
      real*8  lamron(_SFIT_,_SFIT_)
      real*8  vector(_SFIT_)
      real*8  csum, trum
      real    coeffs(_SFIT_)
 
      real    lsfj(_JDIM_)
 
      integer m
      integer n
 
      integer o, oo
      real    ri, fi
      real*8  w, ww
      integer ii
      integer i
      integer s, is
      integer j
 
      real*8  wavo(_ODIM_)
      real*8  wwvo(_ODIM_)
      real*8  rs, r8s
      real*8  mtot
      integer i1, i2, i_w
      integer o1, o2, oi
      real    rspline
 
      ww1 = ww_s(s1)
      ww2 = ww_s(s2)
 
      Ss = s2-s1+1
 
      if (Ss.gt._SFIT_) then
         write(6,*) '         '
         write(6,*) 'ONLY HARD-CODED TO FIT PART OF'
         write(6,*) 'THE TEMPLATE AT ONCE:       '
         write(6,*) '     s1: ',s1
         write(6,*) '     s2: ',s2
         write(6,*) '     Ss: ',Ss
         write(6,*) ' _SFIT_: ',_SFIT_
         write(6,*) '         '
         stop
         endif
 
 
 
      Ls = 0
      do n = 01, Ns
         do m = MMIN, MMAX
            if (wav_imn(_RATB_,m,n)/(1+zn(n)).gt.ww1) goto 888 ! will need to do a better job
            if (wav_imn(_RATY_,m,n)/(1+zn(n)).lt.ww2) goto 888 ! deremining which pix should contrib
            i1 = i_w(ww1*(1+zn(n)),wav_imn(1,m,n))
            i2 = i_w(ww2*(1+zn(n)),wav_imn(1,m,n))
c
c we now know that the pixels
c between i1 and i1 for (n,m)
c
 
 
c
c custom-make lsf...
c
c old       call lsfqk2lsfj_i((i1+i2)/2,lsfqkmn(1,1,2,1),lsfj)
            call lsfqk2lsfj_i((i1+i2)/2,lsfqkmn(1,1,m,n),lsfj)
 
 
c
c set up the arrays needed to multiply and
c integrate over the iodine
c
            o1 = 101 + (i1-1)*10
            o2 = 101 + (i2-1)*10
            do o = max(00001,o1-150), min(_ODIM_,o2+150)
               ri = 1 + (o-101)*0.1
               ii = int(ri)
               fi = ri-ii
               w  = wav_imn(ii,m,n)+fi*(wav_imn(ii+1,m,n)
     .                                 -wav_imn(ii  ,m,n))
               ww = w/(1+zn(n))
               wavo(o) = w
               wwvo(o) = ww
               ftso(o) = ftsw(w)
               enddo
c
c Go through each pixel; find the
c projection of each template node
c on the pixel.  Essentially,
c construct the pixel value as a
c linear combination of template-node
c values.  This allows a straightforward
c least-squares solution.
c
            do i = i1-1, i2+1
               if (rej_imn(i,m,n).ne.1.and.
     .             rej_imn(i,m,n).ne.9) then
                   Ls = Ls + 1
                   if (Ls.gt._LDIM_) stop 'Ls.gt._LDIM_ '
                   pixl(Ls) = pix_imn(i,m,n)
                   wavl(Ls) = wav_imn(i,m,n)
                   wwvl(Ls) = wavl(Ls)/(1+zn(n))
                   ratl(Ls) = rat_imn(i,m,n)
                   do s = 01, Ss
                      modls(s,Ls) = 0.
                      proj(s) = 0.
                      enddo
                   oi = 101 + (i-1)*10
                   do j = -100, 100
                      ww = wwvo(oi+j)                      ! teml waw at this locn
                      r8s = r8s_ww(ww)                     ! abs locn w/in template
                      rs  = 1 + (r8s-s1)                   ! loc locn w/in template
                      if (rs.lt.01) rs = 1.00
                      if (rs.gt.Ss) rs = Ss
                      call rspline_fp(rs,Ss,proj,sa,sb)    ! projn of each node on pix
                      do s = sa, sb
                         modls(s,Ls)=                      ! mod ~ SUM modls(s)*tem(s)
     .                   modls(s,Ls)+ratl(Ls)*lsfj(101+j)
     .                               *ftso(oi+j)*proj(s)
                         enddo
                      enddo
                   endif
               if (Ls.eq.Ls/0500*0500.and.(.false.)) then
                   mtot = 0.
                   do s = 01, Ss
                      mtot = mtot + modls(s,Ls)
                      enddo
                   write(6,198) n,m,i,oi,Ls,wav_imn(i,m,n)/(1+zn(n)),
     .                          int(pixl(Ls)+0.5),
     .                          1.000,
c    .                          temi_imn(i,m,n),
     .                          ratl(Ls),int(mtot+0.5),
     .                          (int((ratl(Ls)*modls(is,Ls))
     .                                 /1e4+0.5),is=01,Ss)
  198              format(i3.3,1x,i2.2,1x,i4.4,1x,i5.5,1x,i5,1x,
     .                    f9.4,1x,i6,1x,
     .                    f7.5,1x,f7.5,3x,
     .                    i6,5x,51i1)
                   endif
               enddo
  888       continue
            enddo!m
         enddo!n
 
 
c------------------------------------------------------------------
c
c
c
c     write(6,*) 'make normal matrix...   Ss: ',Ss
      do i = 01, Ss
      do j = 01, Ss
         normal(i,j) = 0.0d0
         do L = 1, Ls
            normal(i,j) = normal(i,j) + modls(i,L)*modls(j,L)
            enddo
         enddo
         enddo
 
 
c------------------------------------------------------------------
c
c
c
c     write(6,*) 'call dinvert...'
      call dinvertNxN(normal,lamron,Ss,51)
c     write(6,*) 'exit dinvert...'
 
c------------------------------------------------------------------
c
c
c
      do i = 01, Ss
         vector(i) = 0.0d0
         do L = 1, Ls
            vector(i) = vector(i) + pixl(L)*modls(i,L)
            enddo
         enddo
 
cc------------------------------------------------------------------
cc
cc
cc
      do is = 01, Ss
         csum = 0.
         do j = 01, Ss
            csum = csum + vector(j)*lamron(is,j)
            enddo
         coeffs(is) = csum
         ww = ww_s(s1+is-1)
         loctem(is) = coeffs(is)
         enddo
 
 
c------------------------------------------------------------------
c
c
c
      if (.false.) then
         do L = 1, Ls
            modl(L) = 0.
            do s = 01, Ss
               modl(L) = modl(L) + coeffs(s)*modls(s,L)
               enddo
            ww = wavl(L)
            r8s = r8s_ww(ww)                     ! abs locn w/in template
            rs  = 1 + (r8s-s1)                   ! loc locn w/in template
            if (rs.lt.01) rs = 1.00
            if (rs.gt.Ss) rs = Ss
            issl(L) = rspline(rs,coeffs,Ss)
            write(96,196) L,wavl(L),pixl(L),modl(L),barl(L),
     .                    teml(L),issl(L)
  196       format(i4,1x,f9.4,1x,f8.1,1x,f8.1,1x,f8.1,1x,
     .             f8.6,1x,f8.6)
            enddo
            close(96)
         endif
 
      return
      end
 
 
c------------------------------------------------------
c
c this
c
c
      subroutine find_lsqtem(tems,sa,sb,
     .                       pix_imn,nrm_imn,rej_imn,
     .                       wav_imn,zn,lsfqkmn,
     .                       NIMs,MMIN,MMAX)
      implicit none
 
      real    tems(_SDIM_)
      integer sa, sb
      real    pix_imn(_IDIM_,_MDIM_,NIMs)
      real    nrm_imn(_IDIM_,_MDIM_,NIMs)
      integer*1  rej_imn(_IDIM_,_MDIM_,NIMs)
      real*8  wav_imn(_IDIM_,_MDIM_,NIMs)
      real*8  zn(NIMs)
      real    lsfqkmn(_QDIM_,_KDIM_,_MDIM_,NIMs)
      integer NIMs
      integer MMIN, MMAX
 
      integer s
      integer s1, s2, su
      real    teml(_SDIM_)
      real    temh(_SDIM_)
      real*8  ww_s
      integer u
      integer sau, sbu
      real*8  wmin, wmax
      integer i, m, n
      real*8  r8s_ww
      integer ss

      wmin = 9e9
      wmax = 0e0
      do i = 0001, 4020
      do m = MMIN, MMAX
      do n = 0001, NIMs
         if (wav_imn(i,m,n).gt.wmax) wmax = wav_imn(i,m,n)
         if (wav_imn(i,m,n).lt.wmin) wmin = wav_imn(i,m,n)
         enddo
         enddo
         enddo

      sau = int(r8s_ww(wmin))+1
      sbu = int(r8s_ww(wmax))+1

      sau = sa+1
      sbu = sb-1

      write(6,*) 'FIND LSQTEM: '
      write(6,119) sa,sau,wmin,ww_s(sau)
      write(6,119) sb,sbu,wmax,ww_s(sbu)
  119 format(10x,i6.6,1x,i6.6,3x,f12.6,3x,f12.6)

      do s = 0001, _SDIM_
         temh(s) = tems(s)
         enddo
 
c
c go through the template nodes and solve
c for a new one every 25 (the domain actually
c goes out to 50)
c
      do s = sau, sbu+100,30
         ss = s
         if (ss.gt.sbu-50) ss = sbu-50
c
c this finds a least-squares solution to the
c LSF using a local region, at most 51 nodes
c wide
c
         s1 = ss
         s2 = ss+50
         call find_loctem(teml(s1),s1,s2,
     .                    pix_imn,nrm_imn,rej_imn,
     .                    wav_imn,zn,lsfqkmn,NIMs,MMIN,MMAX)

c
c only the middle part of the solved-for
c domain is good, so only fill the "good"
c array with the meat
c
         do su = s1, s2
            u = 0
            if ((su.ge.s1+10.or.s1.eq.sau).and. 
     .          (su.le.s2-10.or.s2.eq.sbu)) u = 1
  144       format(1x,i6.6,1x,f12.6,1x,f10.6,1x,f10.6,1x,i1)
            if (u.eq.1) tems(su) = teml(su)
            write(44,144) su,ww_s(su),teml(su),temh(su),u
            enddo
         write(6,190) s,s1,s2,ww_s(s1),ww_s(s2),
     .                        teml((s1+s2)/2),
     .                        tems((s1+s2)/2),
     .                        temh((s1+s2)/2),
     .                        teml((s1+s2)/2)-
     .                        tems((s1+s2)/2)
  190    format(i6.6,1x,i6.6,1x,i6.6,1x,f9.4,1x,f9.4,1x,6f10.6)
         enddo
 
      do s = 00001, _SDIM_
         if (.not.(tems(s).gt.-0.1)) tems(s) = 1.000
         if (.not.(tems(s).lt. 9.0)) tems(s) = 1.000
         enddo

      close(44)
      close(43)
 
      return
      end
 
 
 
c=======================================================================
c
c This routine uses least squares to determine the ratio between the
c observations pixi() and the model modi() along the order from pixels
c _WAVZ_ to _WAVA_.  It returns this in the array rati(), which scales
c things such that:   pix(i) ~ rati(i)*modi(i).
c
c It also initializes reji() and iteratively rejects the bad pixels;
c one at a time, but within regions, so it should be relatively fast.
c
c The variable "p" goes from 1 ---> NP to go from i = _WAVA_ ----> _WAVZ_,
c the full range of useful pixel values.  For each pixel, p1 and p2 will
c tell you which nodes are relevant, ie whether proj() has any power there.
c Proj() will tell you how each node point contributes to rat() at that
c location.
c
      subroutine find_wavi(pixi,modi,wavi,wwvi,lsfj,rati,reji)
      implicit none
 
      real   pixi(_IDIM_)           ! INPUT : the array of pixel values
      real   modi(_IDIM_)           ! INPUT : the unscaled model values
      real*8 wavi(_IDIM_)
      real*8 wwvi(_IDIM_)
      real   lsfj(_JDIM_)
      real   rati(_IDIM_)           ! OUTPUT: the scaling for each pixel
      integer*1 reji(_IDIM_)           ! OUTPUT: the determined quality
                                     !         of each pixel
 
      real*8 wavx(_IDIM_)
      real   modr(_IDIM_)           ! the scaled model
      real   modx(_IDIM_)
      real   modd(_IDIM_)
      real   modu(_IDIM_)
      integer i                      ! pixel number
      real*8  rx                     ! the location within the rat array,
                                     !    from 1 ---> _WAVN_
      real     proj(_WAVN_,_IDIM_)   ! this contains the proj of each
                                     !    node on each pixel
      real*8 normal(_WAVN_,_WAVN_)   ! the normal matrix
      real*8 lamron(_WAVN_,_WAVN_)   ! the inverse of the normal matrix
      real*8 identy(_WAVN_,_WAVN_)   ! the inverse of the normal matrix
      real*8 vector(_WAVN_)          ! <label?>
      real*8 coeffs(_WAVN_)          ! the output coefficients!
      integer pi, pj                 ! loopy variables
      integer p, p1, p2              ! the current limits on p
      integer NIT                    ! counter of iterations
      integer NREJ, NTOT             ! current new rejections; total
      integer find_reji              ! return value of new rejections
 
      integer i1p(_WAVN_)
      integer i2p(_WAVN_)
 
      integer o, j
      real    ri
      integer ii
 
      real    temo(_ODIM_)
      real    ftso(_ODIM_)
      real    ftsx(_ODIM_)
      real    wavo(_ODIM_)
      real*8  w, ww
      real    ftsw
      real    sunw
 
 
      do p = 1, _WAVN_
         i1p(p) = _IEND_
         i2p(p) = _IBEG_
         enddo
 
c
c Go through the pixels of the current order and initialize them to
c all be "good", and for each to have initially no projection of node
c points (later, to be quick, we will only add projection to points
c that have any without reinitializing the array)
c
c
      do i = 0001, _IDIM_
         do p = 1, _WAVN_
            proj(p,i) = 0.               ! initialize proj array to zeros
            enddo
         enddo
 
 
 
      do i = _IBEG_, _IEND_
 
         p1 = 0
         p2 = 0
 
c
c Now, deal with the middle portion of the image---the meat.
c Use the inner _WAVN_-2 pionts in the spline, hence the (2,1)
c and the rx-1.  We're only using a spline between points
c (2,_WAVN_-1).  Outside of this, we're just using something
c simple and linear.  It's continuous in value, but not in
c slope or anything else.  I wanted the end points to do something
c reasonable, but to have *no* impact on the rest of the fit.
c
         if (i.ge._WAVB_.and.
     .       i.le._WAVY_) then
            rx = 1 + (_WAVN_-1.)*(  i   -_WAVB_)/
     .                           (_WAVY_-_WAVB_)
            call rspline_fp(rx,_WAVN_,proj(1,i),pi,pj)
            p1 = pi
            p2 = pj
            endif
 
 
c
c save time; if we don't have to go through all the
c pixels for each and every node, this will save gobs
c of time
c
         do p = p1, p2
            if (i.lt.i1p(p)) i1p(p) = i
            if (i.gt.i2p(p)) i2p(p) = i
            enddo
 
c         if (i.le.0050.or.i.ge.3950.or.(.false.))
c     .   write(6,222) i,p1,p2,rx,
c     .                (proj(p,i),p=0001,_WAVN_)
c  222    format(i4.4,1x,i3,1x,i3,1x,f7.3,5x,
c     .          21(1x,f7.5))
 
         enddo
 
 
 
      do o = 00001, 40401
         ri = 1 + (o-101)*0.1
         ii = int(ri)
         if (ii.le.0001) ii = 0001
         if (ii.ge.4019) ii = 4019
         w  = wavi(ii) + (ri-ii)*(wavi(ii+1)-wavi(ii))
         ww = wwvi(ii) + (ri-ii)*(wwvi(ii+1)-wwvi(ii))
         ftso(o) = ftsw(w)
         temo(o) = sunw(ww)
         wavo(o) = w
         enddo
      ftsx(00001) = 0
      ftsx(40401) = 0
      do o = 00002, 40400
         ftsx(o) = (ftso(o+1)-ftso(o-1))/2.00
         enddo
      do i = _WAVB_, _WAVY_
         modu(i) = 0.
         modx(i) = 0.
         o = 101 + (i-1)*10
         do j = -100, 100
            if (o+j.ge.00001.and.o+j.le.40401) then
               modu(i) = modu(i)+lsfj(101+j)*ftso(o+j)*temo(o+j)*rati(i)
               modx(i) = modx(i)+lsfj(101+j)*ftsx(o+j)*temo(o+j)*rati(i)
               endif
            enddo
         enddo
 
 
c----------------------------------------------------------------
c
c compute the normal matrix; the set of equations just look
c like:
c
c         NORMAL(N*N) x COEFFS(N) = VECTOR(N)
c
c Where VECTOR(N) is the right side of the equation, and
c COEFFS(N) is the best-fitting set of parameters
c
      do pi = 1, _WAVN_
      do pj = 1, _WAVN_
         normal(pi,pj) = 0.
         do i = i1p(pi), i2p(pi)            ! only thru relevant pix
            if (reji(i).eq.0)               ! construct the normal
     .      normal(pi,pj) =                 !    matrix
     .      normal(pi,pj) + modx(i)**2*proj(pi,i)*proj(pj,i)
            enddo
         enddo
         enddo
 
c
c compute the inverse of the normal equations; essentially
c the covariance matrix;
c
      call dinvertNxN(normal,lamron,_WAVN_,_WAVN_)
 
 
      do pi = 1, _WAVN_
      do pj = 1, _WAVN_
         identy(pi,pj) = 0.
         do i = 1, _WAVN_
            identy(pi,pj) = identy(pi,pj) + normal(pi,i)*lamron(i,pj)
            enddo
          enddo
          enddo
 
 
 
c
c compute the right-hand side of the normal equations
c
      do pi = 1, _WAVN_
         vector(pi) = 0.0
         do i = _WAVB_, _WAVY_
            if (reji(i).eq.0)
     .      vector(pi) = vector(pi) + (pixi(i)-modu(i))*
     .                                 modx(i)*proj(pi,i)
            enddo
         enddo
 
c
c solve for the optimal coefficients for the nodes
c
      do pi = 1, _WAVN_
         coeffs(pi) = 0.0
         do pj = 1, _WAVN_
            coeffs(pi) = coeffs(pi) + vector(pj)*lamron(pj,pi)
            enddo
         enddo
 
 
      do i = _WAVB_, _WAVY_
         modd(i) = 0.
         do pi = 1, _WAVN_
            modd(i) = modd(i) + coeffs(pi)*proj(pi,i)
            enddo
         enddo
 
      return
      end
 
 
 
 
c---------------------------------------
c
c
      real function rspline_tem(rk,gk)
      implicit none
 
      real*8 rk
      real   gk(140001)
 
      real*8 p(140001)
      integer k, kk
      real*8 x0, x1, x2, x3
 
      kk = int(rk)
      if (kk.lt.000001+03) kk = 000001 + 03
      if (kk.gt.140001-03) kk = 140001 - 03
 
      x0 = 1.
      x1 = rk-kk
      x2 = x1**2
      x3 = x1**3
 
      p(kk-1) = 0.00*x0 - 0.5*x1 + 1.0*x2 - 0.5*x3
      p(kk  ) = 1.00*x0 + 0.0*x1 - 2.5*x2 + 1.5*x3
      p(kk+1) = 0.00*x0 + 0.5*x1 + 2.0*x2 - 1.5*x3
      p(kk+2) = 0.00*x0 + 0.0*x1 - 0.5*x2 + 0.5*x3
 
      rspline_tem = p(kk-1)*gk(kk-1)
     .            + p(kk  )*gk(kk  )
     .            + p(kk+1)*gk(kk+1)
     .            + p(kk+2)*gk(kk+2)
 
      return
      end
 
 
 
 
c---------------------------------------
c
c this returns the location in the template
c domain that corresponds to the input wavelength
c
      real*8 function r8s_ww(ww)
      implicit none
 
      real*8 ww
 
      r8s_ww = 1 + dlog10(ww/(_WW0_))/(_LWSTEP_)/(0.4342942647d0)
 
      return
      end
 
 
c---------------------------------------
c
c this returns the wavelength corresponding
c to the r8s location
c
      real*8 function ww_r8s(r8s)
      implicit none
 
      real*8 r8s
      real*8 r8su
 
      r8su = r8s
      if (r8s.lt.000001) r8su = 1
      if (r8s.gt._SDIM_) r8su = _SDIM_
 
      ww_r8s = (_WW0_)*10**((_LWSTEP_)*(r8su-1)*0.4342942647d0)
 
      return
      end
 
 
c----------------------------------------
c
c    this new version will read in a regularly
c    sampled value...
c
c    for a given wavelength w; it reads it in from
c    a file that I store in a special place
c
      real function sunw(w)
      implicit none
 
      real*8 w
 
      real   suns(1 400 001), ff
      common /sunw_f__/suns
      real*8  rs, fs, ww
      integer ss, s
      real   rspline_tem
      character*255 FILENAME255
      character*255 REFDIR
      common /REFDIR_/REFDIR
 
      logical startsun
      common /startsun__/startsun
      data startsun/.true./
 
      if (startsun) then
         FILENAME255 = trim(REFDIR) // 
     .                '/nsosun.tem'
         write(60,'(2a)'), 'solar_file: ',trim(FILENAME255)
         call flush(60)
         open(19,file=FILENAME255,status='old')
         read(19,*)
         read(19,*)
         do ss = 000 001, 140 001
            read(19,*) s, ww, ff
            suns(s) = ff
            enddo
         close(19)
         startsun = .false.
         endif
 
      if (w.lt.5000.00) then
         sunw = 1.0
         return
         endif
      if (w.gt.6400.00) then
         sunw = 1.0
         return
         endif
 
      rs = 1 + (w-5000.00)*100
      ss = int(rs)
      if (ss.lt.000001) ss = 000001
      if (ss.gt.140000) ss = 140001
      fs = rs-ss
 
      sunw = suns(ss) + fs*(suns(ss+1)-suns(ss))
      sunw = rspline_tem(rs,suns)
 
      return
      end
 
 
c---------------------------------------
c
c this returns the wavelength corresponding
c to the template node location
c
      real*8 function ww_s(s)
      implicit none
 
      integer s
      integer su
      real*8 q
 
      su = s
      if (s.lt.000001) su = 1
      if (s.gt._SDIM_) su = _SDIM_
 
      ww_s = (_WW0_)*10**((_LWSTEP_)*(su-1)*0.4342942647d0)
 
      return
      end
 
 
 
c----------------------------------------
c
c    this new version will read in a regularly
c    sampled value...
c
c    for a given wavelength w; it reads it in from
c    a file that I store in a special place
c
      real function ftsw(w)
      implicit none
 
      real*8 w
 
      real   ftsi(1 400 001)
      common /ftsw_f__/ftsi
      real*8  rs, ff, ww
      integer ii
 
      character*255 FILENAME255
      character*255 REFDIR
      common /REFDIR_/REFDIR
 
      logical start
      common /start__/start
      data start/.true./
 
      if (start) then
         FILENAME255 = trim(REFDIR) // 
     .                '/kecknso52p5.sam'
         write(6,*) '     REFDIR: ',REFDIR
         write(6,*) '       TRIM: ',TRIM(REFDIR)
         write(6,*) 'FILENAME255: ',FILENAME255
         write(60,'(2a)'), 'iodine_file: ',trim(FILENAME255)
         call flush(60)
         open(19,file=FILENAME255,status='old')
         do ii = 0000001, 1 400 001
            read(19,*) ww, ff
            ftsi(ii) = ff
            enddo
         close(19)
         start = .false.
         endif
 
      rs = 1 + (w-4950.00d0)*1000
      ii = int(rs)
      ff = rs-ii
 
      if (ii.lt.0 000 001) then
         ftsw = 1e5
         return
         endif
      if (ii.gt.1 400 000) then
         ftsw = 1e5
         return
         endif
 
      ftsw = ftsi(ii) + ff*(ftsi(ii+1)-ftsi(ii))
 
 
      return
      end
 
 
 
c------------------------------------------------------------
c
c
      integer function i_w(ww,wwi)
      implicit none
      real*8 ww
      real*8 wwi(_IDIM_)
 
      integer ia, ib, ic
 
      ia = _IBEG_
      ic = _IEND_
    1 ib = (ia+ic)/2
      if (wwi(ib).le.ww) ia = ib
      if (wwi(ib).ge.ww) ic = ib
      if (ic.gt.ia+1) goto 1
 
      i_w = ib
 
      return
      end
 
 
c#include "51.ROUTINES/find_lsfq_sp.f"
 
c-----------------------------------------------
c
c this routine finds the lsfq that matches
c    the pix array to the ftso array
c
      subroutine find_lsfq_sp(pixi,reji,rati,fxto,lsfqk,lsfji)
      implicit none
 
      real   pixi(_IDIM_)
      integer*1 reji(_IDIM_)
      real   rati(_IDIM_)
      real   fxto(_ODIM_)  ! contains the iodine and temlate
      real   lsfqk(_QDIM_,_KDIM_)
      real   lsfji(_JDIM_,_IDIM_)
 
      integer o, j, i
      integer q, qq
      integer k, kk
      real*8  rk
      real*8  a_iqk(_IDIM_,_QDIM_,_KDIM_)
      real    proji_k(_KDIM_), projk(_KDIM_)
      integer k1_i(_IDIM_), k1
      integer k2_i(_IDIM_), k2
 
      real*8 normal(_QDIM_,_KDIM_,_QDIM_,_KDIM_)
      real*8 lamron(_QDIM_,_KDIM_,_QDIM_,_KDIM_)
      real*8 identy(_QDIM_,_KDIM_,_QDIM_,_KDIM_)
      real*8 vector(_QDIM_,_KDIM_)
      real*8 coeffs(_QDIM_,_KDIM_), csum, csumq(_QDIM_)
 
      real   lsfq(_QDIM_)
      real   lsfjk(_JDIM_,_KDIM_)
 
      integer*2  le_jq_x10000(13,201)
      common /le_jq_x10000_/le_jq_x10000
 
      logical qk_use(_QDIM_,_KDIM_), unuse
 
 
c--------------------------------------
c
c go through and compute how the model is
c consctructed by linear superpositions
c of (q,k) functions
c
      do i = _WAVB_, _WAVY_
         if (reji(i).eq.0) then
            o = 101 + (i-1)*10
            rk = 1 + (_WAVN_-1.0)*
     .               (  i   -_WAVB_)/
     .               (_WAVY_-_WAVB_)
            do k = 1, _KDIM_
               proji_k(k) = 0.0
               enddo
            call rspline_fp(rk,_WAVN_,proji_k,k1_i(i),k2_i(i))
            do q = 0001, _QDIM_
            do k = 0001, _KDIM_
               a_iqk(i,q,k) = 0.
               do j = -100, 100
                  a_iqk(i,q,k) =
     .            a_iqk(i,q,k) + rati(i)*fxto(o+j)*
     .                           le_jq_x10000(q,101+j)/1.0e5*
     .                           proji_k(k)
                  enddo!j
               enddo!k
               enddo!q
            endif
         enddo!i
 
 
c
c come back here if we reject new (qk) nodes
c
c
      do q = 0001, _QDIM_
      do k = 0001, _KDIM_
         qk_use(q,k) = .true.
         enddo
         enddo
 
    1 continue
 
 
c-------------------------------
c
c compute the normal equations
c
      do  q = 0001, _QDIM_
      do  k = 0001, _KDIM_
         do  qq =  0001, _QDIM_
         do  kk =  0001, _KDIM_
             normal(q,k,qq,kk) = 0.0d0
             do i = _WAVB_, _WAVY_
                if (reji(i).eq.0)
     .          normal(q,k,qq,kk) =
     .          normal(q,k,qq,kk) + a_iqk(i,q,k)*a_iqk(i,qq,kk)
                enddo!i
            enddo!kk
            enddo!qq
         enddo!k
         enddo!q
 
      call dinvertNxN_0OK(normal,lamron,(_QDIM_)*(_KDIM_),
     .                                  (_QDIM_)*(_KDIM_))
 
 
c------------------------------
c
c find the vector of condition
c
      do q = 01, _QDIM_
      do k = 01, _KDIM_
         vector(q,k) = 0.0d0
         do i = _WAVB_, _WAVY_
            if (reji(i).eq.0) vector(q,k) =
     .                        vector(q,k) + pixi(i)*a_iqk(i,q,k)
            enddo!i
         enddo!k
         enddo!q
 
 
c------------------------------
c
c find the solution
c
      do k = 01, _KDIM_
      do q = 01, _QDIM_
         coeffs(q,k) = 0.0d0
         do qq = 01, _QDIM_
         do kk = 01, _KDIM_
            coeffs(q,k) = coeffs(q,k) +
     .                    lamron(q,k,qq,kk)*vector(qq,kk)
            enddo
            enddo
         lsfqk(q,k) = coeffs(q,k)
         enddo
         write(6,115) k,(qk_use(q,k),q=1,_QDIM_),
     .                  ( lsfqk(q,k),q=1,_QDIM_)
  115    format(5x,'LSFQ_K',i1,2x,13l1,2x,
     .                    6(1x,f7.5),3x,f7.5,3x,
     .                    6(f7.5,1x))
         enddo
 
      unuse = .false.
      do k = 01, _KDIM_
         if ((     qk_use(04,k)).and.
     .       (.not.qk_use(03,k)).and.
     .       (.not.qk_use(02,k)).and.
     .       (.not.qk_use(01,k))) then
             if (lsfqk(04,k).lt.0.or.
     .           lsfqk(04,k).gt.lsfqk(05,k)) then
                 unuse = .true.
                 qk_use(04,k) = .false.
                 do i = 1, _IDIM_
                    a_iqk(i,04,k) = 0.
                    enddo
                 endif
             endif
         if ((     qk_use(10,k)).and.
     .       (.not.qk_use(11,k)).and.
     .       (.not.qk_use(12,k)).and.
     .       (.not.qk_use(13,k))) then
             if (lsfqk(10,k).lt.0.or.
     .           lsfqk(10,k).gt.lsfqk(09,k)) then
                 unuse = .true.
                 qk_use(10,k) = .false.
                 do i = 1, _IDIM_
                    a_iqk(i,10,k) = 0.
                    enddo
                 endif
             endif
         if ((     qk_use(03,k)).and.
     .       (.not.qk_use(02,k)).and.
     .       (.not.qk_use(01,k))) then
             if (lsfqk(03,k).lt.0.or.
     .           lsfqk(03,k).gt.lsfqk(04,k)) then
                 unuse = .true.
                 qk_use(03,k) = .false.
                 do i = 1, _IDIM_
                    a_iqk(i,03,k) = 0.
                    enddo
                 endif
             endif
         if ((     qk_use(11,k)).and.
     .       (.not.qk_use(12,k)).and.
     .       (.not.qk_use(13,k))) then
             if (lsfqk(11,k).lt.0.or.
     .           lsfqk(11,k).gt.lsfqk(10,k)) then
                 unuse = .true.
                 qk_use(11,k) = .false.
                 do i = 1, _IDIM_
                    a_iqk(i,11,k) = 0.
                    enddo
                 endif
             endif
         if ((     qk_use(02,k)).and.
     .       (.not.qk_use(01,k))) then
             if (lsfqk(02,k).lt.0.or.
     .           lsfqk(02,k).gt.lsfqk(03,k)) then
                 unuse = .true.
                 qk_use(02,k) = .false.
                 do i = 1, _IDIM_
                    a_iqk(i,02,k) = 0.
                    enddo
                 endif
             endif
         if ((     qk_use(12,k)).and.
     .       (.not.qk_use(13,k))) then
             if (lsfqk(12,k).lt.0.or.
     .           lsfqk(12,k).gt.lsfqk(11,k)) then
                 unuse = .true.
                 qk_use(12,k) = .false.
                 do i = 1, _IDIM_
                    a_iqk(i,12,k) = 0.
                    enddo
                 endif
             endif
         if (qk_use(01,k)) then
             if (lsfqk(01,k).lt.0.or.
     .           lsfqk(01,k).gt.lsfqk(02,k)) then
                 unuse = .true.
                 qk_use(01,k) = .false.
                 do i = 1, _IDIM_
                    a_iqk(i,01,k) = 0.
                    enddo
                 endif
             endif
         if (qk_use(13,k)) then
             if (lsfqk(13,k).lt.0.or.
     .           lsfqk(13,k).gt.lsfqk(12,k)) then
                 unuse = .true.
                 qk_use(13,k) = .false.
                 do i = 1, _IDIM_
                    a_iqk(i,13,k) = 0.
                    enddo
                 endif
             endif
c        write(6,'(i4,3x,13l1)') k,(qk_use(q,k),q=1,13)
         enddo
      write(6,*) ' '
 
c
c don't have isolated "good" points...
c
      do q = 01, _QDIM_
         if (     qk_use(q, 0001   ).and.
     .      (.not.qk_use(q, 0002   ))) then
                  qk_use(q, 0001) = .false.
                  do i = 1, _IDIM_
                     a_iqk(i,q, 0001 ) = 0.
                     enddo
                  endif
         if (     qk_use(q,_KDIM_  ).and.
     .      (.not.qk_use(q,_KDIM_-1))) then
                  qk_use(q,_KDIM_)=.false.
                  do i = 1, _IDIM_
                     a_iqk(i,q,_KDIM_) = 0.
                     enddo
                  endif
         do k = 2, _KDIM_-1
            if (qk_use(q,k)         .and.
     .          (.not.qk_use(q,k-1)).and.
     .          (.not.qk_use(q,k+1))) then
                  qk_use(q,k) = .false.
                  do i = 1, _IDIM_
                     a_iqk(i,q,k) = 0.
                     enddo
                  endif
            enddo
         enddo
 
c      if (unuse) goto 1
 
      do k = 01, _KDIM_
         call lsfq2lsfj(lsfqk(1,k),lsfjk(1,k))
         enddo
 
      do i = 0001, _IDIM_
         rk = 1 + (_WAVN_-1.0)*(  i   -_WAVB_)/
     .                         (_WAVY_-_WAVB_)
         if (rk.lt. 0001 ) rk = 1.
         if (rk.gt._WAVN_) rk = _WAVN_
         call rspline_fp(rk,_WAVN_,projk,k1,k2)
         do q = 1, _QDIM_
            lsfq(q) = 0.
            do k = k1, k2
               lsfq(q) = lsfq(q) + projk(k)*lsfqk(q,k)
               enddo
            enddo
         call lsfq2lsfj(lsfq,lsfji(1,i))
         enddo
 
      return
      end
 
 
c-----------------------------------------------
c
c this will find what I used to call the ftsv
c function; it is the template multipied by
c the iodine, sampled with the regular pixel
c scale, but with x10 resolution, and going
c from -10 to +10 of the i range, so that
c we can integrate over the LSF
c
      subroutine find_fxto(fxto,wavi,velu,tems)
      implicit none
 
      real   fxto(_ODIM_)
      real*8 wavi(_IDIM_)
      real*8 velu
      real   tems(_SDIM_)
 
      real*8  w, ww
      integer o
      real ftsw
      real temw
      real r, f
      integer i
 
      do o = 0001, _ODIM_
         r = 1 + (o-101)*0.1
         i = int(r)
         f = r-i
         if (i.lt. 0001   ) i =  0001
         if (i.gt._IDIM_-1) i = _IDIM_-1
         w  = wavi(i) + f*(wavi(i+1)-wavi(i))
         ww = w/(1+velu)
         fxto(o) = ftsw(w)*temw(ww,tems)
         enddo
 
      return
      end
 
 
c-----------------------------------------------
c
c this will find what things will look like
c x10 supersampled with no template multiplication;
c just the iodine at the observed wavelength
c
      subroutine find_ftso(fxto,wavi)
      implicit none
 
      real   fxto(_ODIM_)
      real*8 wavi(_IDIM_)
 
      real*8  w
      integer o
      real ftsw
      real r, f
      integer i
 
      do o = 0001, _ODIM_
         r = 1 + (o-101)*0.1
         i = int(r)
         f = r-i
         if (i.lt. 0001   ) i =  0001
         if (i.gt._IDIM_-1) i = _IDIM_-1
         w  = wavi(i) + f*(wavi(i+1)-wavi(i))
         fxto(o) = ftsw(w)
         enddo
 
      return
      end
 
c-------------------------------------------
c
c this will convert the parametrized LSF into
c      the fully sampled version
c
      subroutine lsfq2lsfj(lsfq,lsfj)
      implicit none
 
      real  lsfq(013)
      real  lsfj(201)
      real*8 summ, sumj
 
      integer j, q, k
 
      integer*2  le_jq_x10000(13,201)
      common /le_jq_x10000_/le_jq_x10000
 
      summ = 0.
      do j = 001, 201
         sumj = 0.
         do q = 01, 13
            sumj = sumj + lsfq(q)*le_jq_x10000(q,j)/1.0e5
            enddo
         lsfj(j) = sumj
         summ = summ + sumj
         enddo
      do j = 1, 201
         lsfj(j) = lsfj(j)/summ
         enddo
      do q = 1, 13
         lsfq(q) = lsfq(q)/summ
         enddo
 
      return
      end
 
 
c------------------------------------------
c
c
c
      subroutine find_modi(modi,fxto,lsfji)
      implicit none
      real modi(_IDIM_)
      real fxto(_ODIM_)
      real lsfji(_JDIM_,_IDIM_)
 
      integer i, o, j
 
      do i = 0001, _IDIM_
         o = 101 + (i-1)*10
         modi(i) = 0.
         do j = 001, _JDIM_
            modi(i) = modi(i) + fxto(o+j-_JCEN_)*lsfji(j,i)
            enddo
         enddo
 
      return
      end
 
c------------------------------------------
c
c
c
c      subroutine find_modi_(modi,fxto,lsfj)
c      implicit none
c      real modi(_IDIM_)
c      real fxto(_ODIM_)
c      real lsfj(_JDIM_)
c
c      integer i, o, j
c
c      do i = 0001, _IDIM_
c         o = 101 + (i-1)*10
c         modi(i) = 0.
c         do j = -100, 100
c            modi(i) = modi(i) + fxto(o+j)*lsfj(101+j)
c            enddo
c         enddo
c
c      return
c      end
c
 
c------------------------------------------------------
c
c this routine pins down the wavelength solution and gets
c us within the ballpark for the LSF shape, searching 
c for (g,d,s),
c
c    "g" is the gaussian that best describes the lsf,
c    "d" is the offset of the spectrum from nominal in wavelength,
c    "s" is the off-nominal slope for the spectrum
c
c------------------------------------------------------
 
 
      subroutine find_initmod(pixi,wavi,gopt,dopt,sopt,m)
      implicit none
 
      real   pixi(4021)
      real*8 wavi(4021)
      integer gopt
      integer dopt
      integer sopt
      integer m
 
      integer k, i
      real*8  sump
      integer sumi
      real*8  wavo(40401)
      real    ftso(40401), ftsw
      real    modo(40401)
      real    e
      integer g, d, s
      integer o, oo
      real    ri
      real*8  ro
      integer ii
      real    ratu
      real    lsfj_gg
      real    psum, msum
      real    emin
      integer gmin, dmin, smin
      real    modr
      real*8  wav_rim_frompoly
 
 
c
c Set up the input spectum which will be broadened
c by the LSF and possibly shifted or stretched
c
c
      do o = 00001, 40401
         ri = 1 + (o-101)*0.1
         ii = int(ri)
         if (ii.lt.0001) ii = 0001
         if (ii.gt.4020) ii = 4020
         wavo(o) = wavi(ii) + (ri-ii)*(wavi(ii+1)-wavi(ii))
         wavo(o) = wav_rim_frompoly(ri,m)
         ftso(o) = ftsw(wavo(o))
         enddo
 
      emin = 9e9
 
c---------------------------------------------------------
c
c loop through possible values of g, d, and s ; I have found
c it easiest to do a crude search in (d,s) then a refined
c search.   We don't need infinite precision; we really only
c need the offset/slope to 0.3 pixel.... for now.  This will
c be improved in the next stage.  We just need to get close
c enough to find the minimum with least squares.
c
      do g =   03,  23, 05
         do o = 00001, 40401
            modo(o) = 0.
            do oo = max(00001,o-100), min(40401,o+100)
               modo(o) = modo(o) + ftso(oo)*lsfj_gg(101+oo-o,g,g)
               enddo
            enddo
 
         do d = -100, 100, 20
         do s = -100, 100, 20
            e = 0.
            do i = 0001+100, 4021-100, 7
               o  = 101 + (i-1)*10 + d + (o-20201)*s/20201
               psum = 0.
               msum = 0.
               do ii = -3, +3
                  psum = psum + pixi(i+ii)
                  msum = msum + modo(o+ii*10)
                  enddo
               ratu = psum/msum
               do ii = -3, 3
                  e = e + abs(pixi(i+ii)-ratu*modo(o+ii*10))
                  enddo
               enddo!i
            if (e.lt.emin) then
               emin = e
               smin = s
               dmin = d
               gmin = g
               endif
            enddo!d
            enddo!s
 
         dopt = dmin
         sopt = smin
 
         do d = dopt-15, dopt+15, 03
         do s = sopt-15, sopt+15, 03
            e = 0.
            do i = 0001+100, 4021-100, 7
               o = 101 + (i-1)*10
               o = o + d + (o-20201)*s/20201
               psum = 0.
               msum = 0.
               do ii = -3, +3
                  psum = psum + pixi(i+ii)
                  msum = msum + modo(o+ii*10)
                  enddo
               ratu = psum/msum
               do ii = -3, 3
                  e = e + abs(pixi(i+ii)-ratu*modo(o+ii*10))
                  enddo
               enddo!i
            if (e.lt.emin) then
               emin = e
               smin = s
               dmin = d
               gmin = g
               endif
            enddo!d
            enddo!s
 
         enddo!g
 
      gopt = gmin
      dopt = dmin
      sopt = smin
 
c------------------------------------------
c
c compute the WLS for the current solution
c
      do i = 0001, 4021
         o  = 101 + (i-1)*10
         ro = o + dopt + (o-20201.0)*sopt/20201.0
         oo = int(ro)
         if (oo.lt.00001) oo = 00001
         if (oo.gt.40400) oo = 40400
         wavi(i) = wavo(oo) + (ro-oo)*(wavo(oo+1)-wavo(oo))
c
c redo... new way
c
         ri = 1 + (ro-101)*0.1
         wavi(i) = wav_rim_frompoly(ri,m)
         enddo
 
      return
      end
 
 
      subroutine polyfitk(xi,xo,
     .                    APOLY,BPOLY,CPOLY,DPOLY,EPOLY)
      implicit none
 
      real*8 xi(41)
      real*8 xo(41)
      real*8 APOLY, BPOLY,CPOLY, DPOLY, EPOLY

      real*8 xl(41)
 
      real*8 abcde(5,41)
      data abcde/
     .  0.0323313d0, 0.1467147d0,-0.4554908d0,-0.3437667d0, 0.6875333d0,
     .  0.0085082d0, 0.0852787d0,-0.1710239d0,-0.2406367d0, 0.3437667d0,
     . -0.0082901d0, 0.0321658d0, 0.0434336d0,-0.1507285d0, 0.0793308d0,
     . -0.0189799d0,-0.0130620d0, 0.1966324d0,-0.0733462d0,-0.1155167d0,
     . -0.0244280d0,-0.0508429d0, 0.2968503d0,-0.0077939d0,-0.2499916d0,
     . -0.0254516d0,-0.0816148d0, 0.3518920d0, 0.0466242d0,-0.3327830d0,
     . -0.0228183d0,-0.1058158d0, 0.3690892d0, 0.0906041d0,-0.3720534d0,
     . -0.0172464d0,-0.1238839d0, 0.3553006d0, 0.1248416d0,-0.3754388d0,
     . -0.0094044d0,-0.1362574d0, 0.3169117d0, 0.1500326d0,-0.3500484d0,
     .  0.0000884d0,-0.1433741d0, 0.2598354d0, 0.1668730d0,-0.3024650d0,
     .  0.0106627d0,-0.1456721d0, 0.1895113d0, 0.1760586d0,-0.2387446d0,
     .  0.0217983d0,-0.1435896d0, 0.1109059d0, 0.1782855d0,-0.1644167d0,
     .  0.0330246d0,-0.1375645d0, 0.0285130d0, 0.1742493d0,-0.0844841d0,
     .  0.0439208d0,-0.1280350d0,-0.0536469d0, 0.1646461d0,-0.0034230d0,
     .  0.0541154d0,-0.1154390d0,-0.1320262d0, 0.1501717d0, 0.0748169d0,
     .  0.0632863d0,-0.1002147d0,-0.2035503d0, 0.1315221d0, 0.1468127d0,
     .  0.0711613d0,-0.0828001d0,-0.2656176d0, 0.1093930d0, 0.2096679d0,
     .  0.0775174d0,-0.0636333d0,-0.3160995d0, 0.0844803d0, 0.2610129d0,
     .  0.0821813d0,-0.0431523d0,-0.3533406d0, 0.0574800d0, 0.2990044d0,
     .  0.0850292d0,-0.0217952d0,-0.3761583d0, 0.0290879d0, 0.3223259d0,
c modify the following to get the average of A to come out to be 1.0000
c    .  0.0859867d0, 0.0000000d0,-0.3838430d0, 0.0000000d0, 0.3301875d0,
     .  0.0859871d0, 0.0000000d0,-0.3838430d0, 0.0000000d0, 0.3301875d0,
     .  0.0850292d0, 0.0217952d0,-0.3761583d0,-0.0290879d0, 0.3223259d0,
     .  0.0821813d0, 0.0431523d0,-0.3533406d0,-0.0574800d0, 0.2990044d0,
     .  0.0775174d0, 0.0636333d0,-0.3160995d0,-0.0844803d0, 0.2610129d0,
     .  0.0711613d0, 0.0828001d0,-0.2656176d0,-0.1093930d0, 0.2096679d0,
     .  0.0632863d0, 0.1002147d0,-0.2035503d0,-0.1315221d0, 0.1468127d0,
     .  0.0541154d0, 0.1154390d0,-0.1320262d0,-0.1501717d0, 0.0748169d0,
     .  0.0439208d0, 0.1280350d0,-0.0536469d0,-0.1646461d0,-0.0034230d0,
     .  0.0330246d0, 0.1375645d0, 0.0285130d0,-0.1742493d0,-0.0844841d0,
     .  0.0217983d0, 0.1435896d0, 0.1109059d0,-0.1782855d0,-0.1644167d0,
     .  0.0106627d0, 0.1456721d0, 0.1895113d0,-0.1760586d0,-0.2387446d0,
     .  0.0000884d0, 0.1433741d0, 0.2598354d0,-0.1668730d0,-0.3024650d0,
     . -0.0094044d0, 0.1362574d0, 0.3169117d0,-0.1500326d0,-0.3500484d0,
     . -0.0172464d0, 0.1238839d0, 0.3553006d0,-0.1248416d0,-0.3754388d0,
     . -0.0228183d0, 0.1058158d0, 0.3690892d0,-0.0906041d0,-0.3720534d0,
     . -0.0254516d0, 0.0816148d0, 0.3518920d0,-0.0466242d0,-0.3327830d0,
     . -0.0244280d0, 0.0508429d0, 0.2968503d0, 0.0077939d0,-0.2499916d0,
     . -0.0189799d0, 0.0130620d0, 0.1966324d0, 0.0733462d0,-0.1155167d0,
     . -0.0082901d0,-0.0321658d0, 0.0434336d0, 0.1507285d0, 0.0793308d0,
     .  0.0085082d0,-0.0852787d0,-0.1710239d0, 0.2406367d0, 0.3437667d0,
     .  0.0323313d0,-0.1467147d0,-0.4554908d0, 0.3437667d0, 0.6875333d0/
      real*8 a,b,c,d,e
      integer k
      real*8 fk
      real*8 asum,abar
 
      asum = 0.0d0
      do k = 1, 41
         asum = asum + xi(k)
         enddo
      abar = asum/41.0
 
      a = 0.
      b = 0.
      c = 0.
      d = 0.
      e = 0.
      do k = 1, 41
         a = a + abcde(1,k)*(xi(k)-abar)
         b = b + abcde(2,k)*(xi(k)-abar)
         c = c + abcde(3,k)*(xi(k)-abar)
         d = d + abcde(4,k)*(xi(k)-abar)
         e = e + abcde(5,k)*(xi(k)-abar)
         enddo
 
      do k = 1, 41
         fk = (k-21)/20.
         xo(k) = abar +
     .           A + B*fk +
     .               C*fk*fk +
     .               D*fk*fk*fk +
     .               E*fk*fk*fk*fk
         enddo
 
      APOLY = A + ABAR
      BPOLY = B
      CPOLY = C
      DPOLY = D
      EPOLY = E
 
      return
      end
 
 
 
c---------------------------------------
c
c will take an array with 41 elements and
c    will perform a spline to it, returning
c    the value of the rk_th element ; this is
c    the frequency with with I sampled the
c    original WLS
c
      real*8 function dsplinek(rk,gk)
      implicit none
 
      real*8 gk(41)
      real*8 rk
 
      real*8 p(41)
      integer k, kk
      real*8 x0, x1, x2, x3
 
      kk = int(rk)
      if (kk.lt.03) kk = 03
      if (kk.gt.39) kk = 39
 
      x0 = 1.
      x1 = rk-kk
      x2 = x1**2
      x3 = x1**3
 
      p(kk-1) = 0.00*x0 - 0.5*x1 + 1.0*x2 - 0.5*x3
      p(kk  ) = 1.00*x0 + 0.0*x1 - 2.5*x2 + 1.5*x3
      p(kk+1) = 0.00*x0 + 0.5*x1 + 2.0*x2 - 1.5*x3
      p(kk+2) = 0.00*x0 + 0.0*x1 - 0.5*x2 + 0.5*x3
 
      dsplinek = p(kk-1)*gk(kk-1)
     .         + p(kk  )*gk(kk  )
     .         + p(kk+1)*gk(kk+1)
     .         + p(kk+2)*gk(kk+2)
 
      return
      end
 
 
      subroutine interp_rk2ri(vk,vi)
      implicit none
 
      real vk(41)
      real vi(4021)
 
      integer ik, i
      real    rk, fk
 
      do i = 1, 4021
         rk = 1 + (i-1)*0.01
         if (rk.le. 1.00) rk =  1.00
         if (rk.gt.40.95) rk = 40.95
         ik = int(rk)
         fk = rk-ik
         vi(i) = vk(ik) + fk*(vk(ik+1)-vk(ik))
         enddo
 
      return
      end
 
 
      subroutine smoo_ik(vik,vsk)
      implicit none
 
      integer vik(41)
      real    vsk(41)
      integer k
 
      vsk(01) = (vik(01))/1.
      vsk(02) = (vik(01)+vik(02)+vik(03))/3.
      vsk(03) = (vik(01)+vik(02)+vik(03)+vik(04)+vik(05))/5.
      vsk(41) = (vik(41))/1.
      vsk(40) = (vik(41)+vik(40)+vik(39))/3.
      vsk(39) = (vik(41)+vik(40)+vik(39)+vik(38)+vik(37))/5.
 
      do k = 04, 38
         vsk(k) = -0.095*vik(k-3)
     .          +  0.143*vik(k-2)
     .          +  0.286*vik(k-1)
     .          +  0.333*vik(k  )
     .          +  0.286*vik(k+1)
     .          +  0.143*vik(k+2)
     .          -  0.095*vik(k+3)
         enddo
 
      return
      end
 
 
 
      subroutine smoo_rk(vik,vsk)
      implicit none
 
      real vik(41)
      real vsk(41)
      integer k
 
      vsk(01) = (vik(01))/1.
      vsk(02) = (vik(01)+vik(02)+vik(03))/3.
      vsk(03) = (vik(01)+vik(02)+vik(03)+vik(04)+vik(05))/5.
      vsk(41) = (vik(41))/1.
      vsk(40) = (vik(41)+vik(40)+vik(39))/3.
      vsk(39) = (vik(41)+vik(40)+vik(39)+vik(38)+vik(37))/5.
 
      do k = 04, 38
         vsk(k) = -0.095*vik(k-3)
     .          +  0.143*vik(k-2)
     .          +  0.286*vik(k-1)
     .          +  0.333*vik(k  )
     .          +  0.286*vik(k+1)
     .          +  0.143*vik(k+2)
     .          -  0.095*vik(k+3)
         enddo
 
      return
      end




c--------------------------------------
c
c
      subroutine smoo_quad5(xon,xsn,Ns)
      implicit none
    
      integer Ns
      real*8 xon(Ns)
      real*8 xsn(Ns)

      real*8 hhh(999)
      real*8 A, B, C
      integer N, NN

      if (Ns.gt.999) stop 'Ns.gt.999 smoo_quad5'

      do N = 1, Ns
         hhh(N) = xon(N)
         enddo
 
      do N = 1, Ns
         NN = N
         if (NN.le.01+2) NN = 01 + 2
         if (NN.ge.Ns-2) NN = Ns - 2
         A = -0.086*hhh(NN-2)
     .     +  0.343*hhh(NN-1)
     .     +  0.486*hhh(NN  )
     .     +  0.343*hhh(NN+1)
     .     -  0.086*hhh(NN+2)
         B = -0.200*hhh(NN-2)
     .     -  0.100*hhh(NN-1)
     .     +  0.000*hhh(NN  )
     .     +  0.100*hhh(NN+1)
     .     +  0.200*hhh(NN+2)
         C =  0.143*hhh(NN-2)
     .     -  0.071*hhh(NN-1)
     .     -  0.144*hhh(NN  )
     .     -  0.071*hhh(NN+1)
     .     +  0.143*hhh(NN+2)
         xsn(N) = A + B*(N-NN) + C*(N-NN)**2
         enddo

      return
      end



c--------------------------------------
c
c
      subroutine smoo_quad7(xon,xsn,Ns)
      implicit none
    
      integer Ns
      real*8 xon(Ns)
      real*8 xsn(Ns)

      real*8 hhh(999)
      real*8 A, B, C
      integer N, NN

      if (Ns.gt.999) stop 'Ns.gt.999 smoo_quad7'

      do N = 1, Ns
         hhh(N) = xon(N)
         enddo
 
      do N = 1, Ns
         NN = N
         if (NN.le.01+3) NN = 01 + 3 
         if (NN.ge.Ns-3) NN = Ns - 3
         A = -0.095*hhh(NN-3)
     .     +  0.143*hhh(NN-2)
     .     +  0.286*hhh(NN-1)
     .     +  0.333*hhh(NN  )
     .     +  0.286*hhh(NN+1)
     .     +  0.143*hhh(NN+2)
     .     -  0.095*hhh(NN+3)
         B = -0.107*hhh(NN-3)
     .     -  0.071*hhh(NN-2)
     .     -  0.036*hhh(NN-1)
     .     +  0.000*hhh(NN  )
     .     +  0.036*hhh(NN+1)
     .     +  0.071*hhh(NN+2)
     .     +  0.107*hhh(NN+3)
         C =  0.060*hhh(NN-3)
     .     +  0.000*hhh(NN-2)
     .     -  0.036*hhh(NN-1)
     .     -  0.048*hhh(NN  )
     .     -  0.036*hhh(NN+1)
     .     +  0.000*hhh(NN+2)
     .     +  0.060*hhh(NN+3)
         xsn(N) = A + B*(N-NN) + C*(N-NN)**2
         enddo

      return
      end

 
c--------------------------------------------------------
c
c
c
      real function lsfj_gg(ju,gu,ggu)
      implicit none
 
      integer ju, gu, ggu
 
      real gsig(25)
      common /common_gsig/gsig
 
 
      logical first_call
      common /first_call/first_call
      data first_call/.true./
 
      real lsfjgg(201,25,25)
      integer j, jj
      integer g1,g2
      real suma
      real sumx
      real dj, rj, fj, ff
 
      real psum, msum, esum
      integer guse, ouse
      real    r, ruse
      real lsfa(201)
      real lsfb(201)
 
      if (first_call) then
         do g1 = 1, 25
         do g2 = 1, 25
            do j = 001, 201
               r = (j-101)*0.1
               if (j.le.101) lsfa(j) = exp(-0.5*(r/gsig(g1))**2)
               if (j.ge.101) lsfa(j) = exp(-0.5*(r/gsig(g2))**2)
               enddo
            suma = 0.
            sumx = 0.
            do j = 001, 201
               lsfb(j) = 0.
               do jj = max(001,j-5), min(201,j+5)
                  ff = 1.0
                  if (abs(j-jj).eq.5) ff = 0.5
                  lsfb(j) = lsfb(j) + ff*lsfa(j)
                  enddo
               suma = suma + lsfb(j)
               sumx = sumx + lsfb(j)*(j-101)
               enddo
            dj = sumx/suma
            do j = 001, 201
               rj = j + dj
               if (rj.lt.001) rj = 001.00
               if (rj.gt.201) rj = 201.00
               jj = int(rj)
               if (jj.lt.001) jj = 001
               if (jj.gt.200) jj = 200
               fj = rj-jj
               lsfjgg(j,g1,g2) = (lsfb(jj)+fj*(lsfb(jj+1)-lsfb(jj)))/
     .                            suma
               enddo
            enddo!g2
            enddo!g1
         first_call = .false.
         endif
 
      lsfj_gg = lsfjgg(ju,gu,ggu)
 
      return
      end
 
 
      real function lsfj_ga(ju,rgu,rau)
      implicit none
 
      integer ju
      real    rgu
      real    rau
 
      real    rg1, rg2
      integer ig1, ig2
      real    fg1, fg2
      real    lsfj_gg
 
      rg1 = rgu-0.5*rau
 
      ig1 = int(rg1)
      if (ig1.lt.01) ig1 = 01
      if (ig1.gt.24) ig1 = 24
      fg1 = rg1 - ig1
 
      rg2 = rgu+0.5*rau
      ig2 = int(rg2)
      if (ig2.lt.01) ig2 = 01
      if (ig2.gt.24) ig2 = 24
      fg2 = rg2 - ig2
 
      lsfj_ga = (1-fg1)*(1-fg2)*lsfj_gg(ju,ig1  ,ig2  )
     .        + (1-fg1)*( fg2 )*lsfj_gg(ju,ig1  ,ig2+1)
     .        + ( fg1 )*(1-fg2)*lsfj_gg(ju,ig1+1,ig2  )
     .        + ( fg1 )*( fg2 )*lsfj_gg(ju,ig1+1,ig2+1)
 
      return
      end
 
 
 
 
      subroutine init_wavim(wav_im,MMIN,MMAX)
      implicit none
      real*8 wav_im(_IDIM_,_MDIM_)
      integer MMIN, MMAX 

      integer m
      real*8  fts0km(41,16)
      common /fts0_/fts0km
      integer i, k
      real*8 rk
      real*8 wavk(41), wavks(41)
      real*8 dsplinek
      real*8 APOLY,BPOLY,CPOLY,DPOLY,EPOLY
      real   xk
      real*8 wav_rim_frompoly


      if (.false.) then ! old way 
      do m = MMIN, MMAX
         do k = 1, 41
            wavk(k) = fts0km(k,m)
            enddo
         call polyfitk(wavk,wavks,
     .                 APOLY,BPOLY,CPOLY,DPOLY,EPOLY)
         do i = 0001, 4021, 1
            rk = 1 + (i-1)/100.
            wav_im(i,m) = dsplinek(rk,wavks)
            xk = (i-2001.0)/2000.0
            write(99,199) i,rk,wav_im(i,m),
     .                      xk,APOLY + 
     .                         BPOLY*xk +
     .                         CPOLY*xk*xk +
     .                         DPOLY*xk*xk*xk +
     .                         EPOLY*xk*xk*xk*xk,
     .                         wav_rim_frompoly(i*1.00,m)
  199       format(1x,i4.4,1x,f10.6,1x,f16.5,'d0',
     .                     1x,f10.7,1x,f16.5,'d0',1x,f16.5,'d0')
            enddo
         enddo
         endif

      if (.true.) then ! old way 
         do m = MMIN, MMAX
            do i = 0001, 4021
               wav_im(i,m) = wav_rim_frompoly(i*1.00,m)
               enddo
            enddo
         endif
 
      return
      end
 
 
c----------------------------------------------------
c
c this routine will take the observations and
c the model and will determine which data points
c are bad, based on the local set of pixels
c
 
      integer function find_locrej(pixi,modi,reji,sigi)
      implicit none
 
      real   pixi(_IDIM_)
      real   modi(_IDIM_)
      integer*1 reji(_IDIM_)
      real   sigi(_IDIM_)
 
      integer i, ic, ii
      real psum, ssum
      real pbar, sbar
      integer nsum
 
 
c
c go through pixel by pixel ; for each pixel
c look at the local +5 to -5 pixels (ignoring
c itself) and see if it's consistent with them
c this is vulnerable to missing bad pixels next
c to other bad pixels
c
      find_locrej = 0
      do i = 0001, _IDIM_
         ic = i
         if (ic.lt. 0001 +6) ic =  0001  + 5
         if (ic.gt._IDIM_-6) ic = _IDIM_ - 5
         psum = 0.
         ssum = 0.
         nsum = 0
         do ii = ic-5, ic+5
            if (ii.ne.i) then
               psum = psum + pixi(ii)
               ssum = ssum + abs(pixi(ii)-modi(ii))
               nsum = nsum + 1
               endif
            enddo
         pbar = psum/nsum
         sbar = ssum/nsum
         sigi(i) = (pixi(i)-modi(i))/sbar
         if (reji(i).ne.9) reji(i) = 0
         if (i.lt._IBEG_) reji(i) = 4
         if (i.gt._IEND_) reji(i) = 4
         if (i.lt._RATB_) reji(i) = 2
         if (i.gt._RATY_) reji(i) = 3
         if (abs(pixi(i)-modi(i)).gt.4.5*sbar) then
            find_locrej = find_locrej + 1
            reji(i) = 1
            endif
         enddo
 
      return
      end
 
 
c
c
c
c
      subroutine gsni2lsfji(gsni,lsfji)
      implicit none
 
      real gsni(_IDIM_)
      real lsfji(_JDIM_,_IDIM_)
 
      real lsfj_ga
      integer i, j
 
      do i = 1, _IDIM_
         do j = 1, _JDIM_
            lsfji(j,i) = lsfj_ga(j,gsni(i),0.)
            enddo
         enddo
 
      return
      end
 
 
      real function temw(ww,tems)
      implicit none
 
      real*8 ww
      real   tems(_SDIM_)
 
      real*8 rs
      real rspline
      real*8 r8s_ww
 
      rs = r8s_ww(ww)
      temw = 1.00
      if (rs.lt. 0001 ) return
      if (rs.gt._SDIM_) return
 
      temw = rspline(rs,tems,_SDIM_)
 
      return
      end
 
 
      integer function iclip(ival,ilo,ihi)
      implicit none
 
      integer ival
      integer ilo
      integer ihi
 
      iclip = ival
      if (iclip.gt.ihi) iclip = ihi
      if (iclip.lt.ilo) iclip = ilo
      if (.not.(iclip.lt.ihi).and.
     .    .not.(iclip.gt.ilo)) iclip = ihi
 
      return
      end
 
 
 
c--------------------------------------------------------
c
c this rouine will do a find for velocity using
c least squares and all the pixels, it may be
c better to do this line-by-line
c
      subroutine tweak_vel(pixim,modim,wavim,
     .                     ratim,rejim,metim,
     .                     tems,velu,lsfjim,
     .                     MMIN, MMAX)
      implicit none
 
      real      pixim(_IDIM_,_MDIM_)
      real      modim(_IDIM_,_MDIM_)
      real*8    wavim(_IDIM_,_MDIM_)
      real      ratim(_IDIM_,_MDIM_)
      integer*1 rejim(_IDIM_,_MDIM_)
      integer*1 metim(_IDIM_,_MDIM_)
      real      tems(_SDIM_)
      real*8    velu
      real      lsfjim(_JDIM_,_IDIM_,_MDIM_)
      integer   MMIN, MMAX
 
      real fxto_m(_ODIM_,_MDIM_)
      real fxtp_m(_ODIM_,_MDIM_)
      real fxtm_m(_ODIM_,_MDIM_)
      real modu_im(_IDIM_,_MDIM_)
      real modx_im(_IDIM_,_MDIM_)
      real*8  w, ww, wwp, wwm
      integer m, i, j, o,oo
      real    ftsu, ftsw, temw
      integer iclip
      real*8  suma, sumb
      integer sumn
      real    ri
      integer ii
      integer NIT
      real*8  veluu
 
      veluu = velu
 
      write(6,*) ' '
      write(6,*) ' '
 
      do NIT = 1, 4
      do m = MMIN, MMAX
         do o = 00001, _ODIM_
            ri = 1.0 + (o-101)*0.1
            ii = iclip(int(ri),_IBEG_,_IEND_-1)
            w  =  wavim(ii,m) + (ri-ii)*(wavim(ii+1,m)-wavim(ii,m))
c           ww =   w/(1+veluu)
c           wwp=  ww/(1+0.1666d-8) ! +0.5 m/s
c           wwm=  ww/(1-0.1666d-8) ! -0.5 m/s
            ww =   w/(1+veluu)
            wwp=   w/(1+veluu+0.1666d-8) ! +0.5 m/s
            wwm=   w/(1+veluu-0.1666d-8) ! -0.5 m/s
            ftsu = ftsw(w)
            fxto_m(o,m) = temw(ww ,tems)*ftsu
            fxtp_m(o,m) = temw(wwp,tems)*ftsu
            fxtm_m(o,m) = temw(wwm,tems)*ftsu
            enddo
         do i = 0001, _IDIM_
            modu_im(i,m) = 0.
            modx_im(i,m) = 0.
            o = 101 + (i-1)*10
            do j = -100, +100
               modu_im(i,m) =
     .         modu_im(i,m) + lsfjim(101+j,i,m)*fxto_m(o+j,m)
               modx_im(i,m) =
     .         modx_im(i,m) + lsfjim(101+j,i,m)*fxtp_m(o+j,m)
     .                      - lsfjim(101+j,i,m)*fxtm_m(o+j,m)
               enddo
            enddo
         enddo
 
      suma = 0.
      sumb = 0.
      sumn = 0
      do m = MMIN, MMAX
c        do i = _RATB_, _RATY_
         do i = 0100, 3900
            if (rejim(i,m).eq.0.and.metim(i,m).eq.0) then
               suma = suma + ratim(i,m)*modx_im(i,m)
     .                      *(pixim(i,m)-ratim(i,m)*modu_im(i,m))
               sumb = sumb + ratim(i,m)*modx_im(i,m)
     .                      *ratim(i,m)*modx_im(i,m)
               sumn = sumn + 1
               endif
            enddo
         enddo
 
      write(6,'(15x,''VELNIT: '',i2,1x,5f15.10,4x,2f8.2,1x,i4.4)')
     .                           NIT,velu,
     .                            veluu-velu,
     .                            veluu,suma/sumb*0.3d-8,
     .                            veluu+suma/sumb*0.3d-8,
     .                            suma/sumb*0.3d-8*3e8,
     .                            (veluu-velu)*3e8, sumn
      veluu = veluu + suma/sumb*0.3d-8
 
      enddo!NIT
 
      velu = veluu
 
      do m = MMIN, MMAX
      do i =  0001 , _IDIM_
         modim(i,m) = modu_im(i,m)
         enddo
         enddo

      if (.false.) 
     .write(19,'(15x,''VELNIT: '',i2,1x,5f15.10,4x,2f8.2)')
     .                           NIT,velu,
     .                            veluu-velu,
     .                            veluu,suma/sumb*0.3d-8,
     .                            veluu+suma/sumb*0.3d-8,
     .                            suma/sumb*0.3d-8*3e8,
     .                            (veluu-velu)*3e8
 
      veluu = veluu + suma/sumb*0.3d-8
 
      return
      end
 
 
c----------------------------------------------------
c
c
c
c
c
 
      subroutine lsfqk2lsfj_i(i,lsfqk,lsfj)
      implicit none
 
      integer i
      real lsfqk(_QDIM_,_KDIM_)
      real lsfj(_JDIM_)
 
      real*8 rk
      integer q
      real lsfq(_QDIM_)
      integer k, k1, k2
      real projk(_KDIM_)
 
      rk = 1 + (_WAVN_-1.0)*
     .         (  i   -_WAVB_)/
     .         (_WAVY_-_WAVB_)
 
      do k = 1, _KDIM_
         projk(k) = 0.0
         enddo
 
      call rspline_fp(rk,_WAVN_,projk,k1,k2)
 
      do q = 1, _QDIM_
         lsfq(q) = 0.
         enddo
      do k = k1, k2
         do q = 1, _QDIM_
            lsfq(q) = lsfq(q) + projk(k)*lsfqk(q,k)
            enddo
         enddo
 
      call lsfq2lsfj(lsfq,lsfj)
 
c      write(6,111) i,rk,k1,k2,(projk(k),k=1,_KDIM_),(lsfq(q),q=1,_QDIM_),
c     .               lsfj(101)
  111 format(i4.4,1x,f4.2,2x,i1,1x,i1,3x,9(1x,f4.2),
     .                                3x,13(1x,f5.3),
     .                                3x,f8.6)
 
      return
      end
 
 
c----------------------------------------------------
c
c
c
c
c
 
      subroutine lsfqk2lsfji(lsfqk,lsfji)
      implicit none
 
      real lsfqk(_QDIM_,_KDIM_)
      real lsfji(_JDIM_,_IDIM_)
 
      real*8 rk
      integer q
      real lsfq(_QDIM_)
      integer k, k1, k2
      real projk(_KDIM_)
      integer i
 
 
      do i = 0001, _IDIM_
 
         rk = 1 + (_WAVN_-1.0)*
     .            (  i   -_WAVB_)/
     .            (_WAVY_-_WAVB_)
 
         do k = 1, _KDIM_
            projk(k) = 0.0
            enddo
 
         call rspline_fp(rk,_WAVN_,projk,k1,k2)
 
         do q = 1, _QDIM_
            lsfq(q) = 0.
            enddo
 
         do k = k1, k2
            do q = 1, _QDIM_
               lsfq(q) = lsfq(q) + projk(k)*lsfqk(q,k)
               enddo
            enddo
 
         call lsfq2lsfj(lsfq,lsfji(1,i))
 
         enddo!i
 
      return
      end
 
 
c--------------------------------------------------------
c
c this routine will read in an rdsk file into a pixel
c    array
c
      subroutine rdsk(RDSKFILE,pix)
      implicit none
 
      character*200 RDSKFILE
      real*8        pix(4021,16)
 
      byte buffb(4)
      integer ios
 
      byte      b4(4)
      character c4(4)
      real      r1
      integer   i1
      integer*2 i2(2)
      equivalence (c4,b4,r1,i1,i2)
      integer irec
 
      integer i,j,n
      character*80 FITSFILE
      integer jold
 
      integer NMOS
      character*99999 CHAROUT
      character*80    CHAR80
      character*19    CHARDATE
      integer NCHAR, NCHARs
      integer NLINE
      character*80 LOGFILE
      byte stream(316518)
      character*10 STRING

      real blaze_abcde_m(5,16)
      common /blaze_abcde_m_/blaze_abcde_m
      data blaze_abcde_m /
     .     0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .     0.404957, -0.005007, -0.271641,  0.000815,  0.073954,
     .     0.436549, -0.035502, -0.264960,  0.012835,  0.068061,
     .     0.447705, -0.037597, -0.270370,  0.007823,  0.075862,
     .     0.454244, -0.036874, -0.246967,  0.007651,  0.048908,
     .     0.464779, -0.034572, -0.265386, -0.000611,  0.069471,
     .     0.471788, -0.033998, -0.264400, -0.006306,  0.066897,
     .     0.471220, -0.038637, -0.256811, -0.005251,  0.064193,
     .     0.471379, -0.039233, -0.249445, -0.012743,  0.065507,
     .     0.460334, -0.031035, -0.241018, -0.018622,  0.068456,
     .     0.464803, -0.026182, -0.264137, -0.020164,  0.085511,
     .     0.450931, -0.025705, -0.240858, -0.013854,  0.070892,
     .     0.441154, -0.020534, -0.235190, -0.017902,  0.073820,
     .     0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .     0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .     0.000000,  0.000000,  0.000000,  0.000000,  0.000000/

      real*8  ff, fx
      real*8  pbar, rbar
      real*8  psum
      integer nsum
      integer m

c---------------------------------------------------------------
c
c
c 
      open(10,file=RDSKFILE,status='old',
     .     err=900,recl=316518,form='UNFORMATTED',
     .     access='DIRECT')
 
c
c this line will read in the entire fits file
c    into long array variable
c
      read(10,rec=1,iostat=ios) stream
 
      irec = irec/4*4 + 1
      jold = 0
      irec = 4*4021*16/4*4 - 8 + 1
      irec = 1
   1  continue
      n = (irec-25)/4 + 1
      j = (n-1)/4021 + 1
      i = n - (j-1)*4021
 
      b4(1) = stream((irec-1)+1)
      b4(2) = stream((irec-1)+2)
      b4(3) = stream((irec-1)+3)
      b4(4) = stream((irec-1)+4)
      if (_LINUX_) then
         b4(1) = stream((irec-1)+4)
         b4(2) = stream((irec-1)+3)
         b4(3) = stream((irec-1)+2)
         b4(4) = stream((irec-1)+1)
         endif
 
      irec = irec + 4
      if (j.le.00) goto 1
      if (i.le.00) goto 1
      if (j.gt.16) goto 2
      pix(i,j) = r1
      if (j.ne.jold) then
         jold = j
         endif
      goto 1
    2 continue


c
c do the avarage blaze normalization, as taken
c from b-stars
c
      do m = 02, 13
         do i = 0001, 4021
            fx = (i-2001.0)/2000.0
            ff = blaze_abcde_m(1,m)
     .         + blaze_abcde_m(2,m)*fx
     .         + blaze_abcde_m(3,m)*fx*fx
     .         + blaze_abcde_m(4,m)*fx*fx*fx
     .         + blaze_abcde_m(5,m)*fx*fx*fx*fx
            ff = ff/blaze_abcde_m(1,m)
            pix(i,m) = pix(i,m)/ff
            enddo
         enddo


c
c for simplicity (?) normalize each spectrum
c to have an average pixel value of 50000
c between i=501 and i=3501
c
      do m = 02, 13
         psum = 0.0d0
         nsum = 0
         do i = 0501, 3501
            psum = psum + pix(i,m)
            nsum = nsum + 1
            enddo
         pbar = psum/nsum
         rbar = pbar/50000
         do i = 0001, 4021
            pix(i,m) = pix(i,m)/rbar
            enddo
         pix(4021,m) = pix(4020,m)
         enddo
 
      return
 
  900 continue
      write(6,'(''               '')')
      write(6,'(''ERROR OPENING  '')')
      write(6,'(''     RDSKFILE: '',200a)') RDSKFILE
      write(6,'(''               '')')
      stop
      end
 
 
 
c--------------------------------------------------------
c
c this routine will read in an rdsk file into a pixel
c    array
c
      subroutine rdsk_DP(RDSKFILE,pix)
      implicit none
 
      character*80 RDSKFILE
      real*8       pix(4021,16)
 
      byte buffb(4)
      integer ios
 
      byte      b8(8)
      character c8(8)
      real*8    r8
      integer   i4(2)
      integer*2 i2(4)
      equivalence (c8,b8,r8,i4,i2)
      integer irec
 
      integer i,j,n
      character*80 FITSFILE
      integer jold
 
      integer NMOS
      integer NARG, NARGs
      character*99999 CHAROUT
      character*80    CHAR80
      character*19    CHARDATE
      integer NCHAR, NCHARs
      integer NLINE
      byte stream(514743)
      character*10 STRING
 
      open(10,file=RDSKFILE,status='old',
     .     err=900,recl=514743,form='UNFORMATTED',
     .     access='DIRECT')
 
 
c
c this line will read in the entire fits file
c    into one variable
c
      read(10,rec=1,iostat=ios) stream
 
c     write(6,'(80a)') (stream(i),i=1,80)
c     if (.true.) stop
 
      irec = irec/4*4 + 1
      jold = 0
      irec = 4*4021*16/4*4 - 8 + 1
      irec = 1
   1  continue
      n = (irec-25)/8 + 1
      j = (n-1)/4021 + 1
      i = n - (j-1)*4021
 
      b8(1) = stream((irec-1)+1)
      b8(2) = stream((irec-1)+2)
      b8(3) = stream((irec-1)+3)
      b8(4) = stream((irec-1)+4)
      b8(5) = stream((irec-1)+5)
      b8(6) = stream((irec-1)+6)
      b8(7) = stream((irec-1)+7)
      b8(8) = stream((irec-1)+8)
      if (_LINUX_) then
         b8(1) = stream((irec-1)+8)
         b8(2) = stream((irec-1)+7)
         b8(3) = stream((irec-1)+6)
         b8(4) = stream((irec-1)+5)
         b8(5) = stream((irec-1)+4)
         b8(6) = stream((irec-1)+3)
         b8(7) = stream((irec-1)+2)
         b8(8) = stream((irec-1)+1)
         endif
 
      irec = irec + 8
      if (j.le.00) goto 1
      if (i.le.00) goto 1
      if (j.gt.16) goto 2
      pix(i,j) = r8
      if (j.ne.jold) then
         jold = j
         endif
      goto 1
    2 continue
 
      return
 
  900 continue
      write(6,*) '---> ERROR: '
      end
 
 
c-------------------------------------------
c
c this will convert the parametrized LSF into
c      the fully sampled version
c
      subroutine lsfj2lsfq(lsfj,lsfq)
      implicit none
 
      real  lsfj(_JDIM_)
      real  lsfq(_QDIM_)
 
      real*8 sumq, sumj
 
      integer j, q, k
 
      integer*2  le_jq_x10000(13,201)
      common /le_jq_x10000_/le_jq_x10000
 
      real*8 normal(13,13)
      real*8 lamron(13,13)
      real*8 vector(13)
      real*8 coeffs(13)
 
      real   lsfje(_JDIM_)
      integer qa, qb
 
      real*8 proj_jq(13,201)
 
      do j = 001, _JDIM_
      do q = 001, _QDIM_
         proj_jq(q,j) = le_jq_x10000(q,j)/1.0e5
         enddo
         enddo
 
      do qa = 1, _QDIM_
      do qb = 1, _QDIM_
         normal(qa,qb) = 0.0d0
         do j = 1, _JDIM_
            normal(qa,qb) = normal(qa,qb) + proj_jq(qa,j)*proj_jq(qb,j)
            enddo
         enddo
         enddo
 
 
      call dinvertNxN(normal,lamron,_QDIM_,_QDIM_)
 
      do q = 1, _QDIM_
         vector(q) = 0.0d0
         do j = 1, _JDIM_
            vector(q) = vector(q) + proj_jq(q,j)*lsfj(j)
            enddo
         enddo
 
      do q = 1, _QDIM_
         coeffs(q) = 0.0d0
         do qa = 1, _QDIM_
            coeffs(q) = coeffs(q) + lamron(q,qa)*vector(qa)
            enddo
         enddo
 
      sumj = 0.0
      do j = 001, _JDIM_
         sumj = sumj + lsfj(j)
         enddo
      sumq = 0.0
      do q = 001, _QDIM_
         lsfq(q) = coeffs(q)
         sumq = sumq + lsfq(q)
         enddo
 
      return
      end
 
 
 
c------------------------------------------------
c
c
      subroutine make_init_tem(temu,nrm_imn,wav_imn,zn,NIMs,
     .                         MMIN,MMAX)
      implicit none
 
      real    temu(_SDIM_)
      real    nrm_imn(_IDIM_,_MDIM_,NIMs)
      real*8  wav_imn(_IDIM_,_MDIM_,NIMs)
      real*8  zn(NIMs)
      integer NIMs
      integer MMIN, MMAX
 
      integer s, i, m, n
      real    tems(_SDIM_)
      real    temt(_SDIM_)
      integer temn(_SDIM_)
      real*8  ww, ww_s, w
      real    wt
      integer i_w
 
      do s = 00001, _SDIM_
         tems(s) = 0.000
         temt(s) = 0.000
         temn(s) = 0
         ww = ww_s(s)
         do n = 1, NIMs
            do m = MMIN, MMAX
               w = ww*(1+zn(n))
               if (wav_imn(_IBEG_,m,n).gt.w) goto 15
               if (wav_imn(_IEND_,m,n).lt.w) goto 15
               i = i_w(w,wav_imn(1,m,n))
               wt = min(100.00,float(i-_IBEG_),float(_IEND_-i))/100.0
               tems(s) = tems(s) + wt*nrm_imn(i,m,n)
               temt(s) = temt(s) + wt
               temn(s) = temn(s) + 1
  15           continue
               enddo
            enddo
         temu(s) = 1.00
         if (temt(s).ge.0.1) temu(s) = tems(s)/temt(s)
         if (s.eq.s/1000*1000)
     .   write(6,155) s,ww_s(s),tems(s),temt(s),temn(s),temu(s)
 155     format('make_init_tem: ',i6,1x,f12.5,1x,
     .                            f11.6,1x,f9.3,1x,i4,1x,f9.6)
         enddo
 
      return
      end
 
 
 
c---------------------------------------------------
c
c
c
      subroutine save_mod(pix_imn,mod_imn,wav_imn,
     .                    nrm_imn,bar_imn,
     .                    rej_imn,tel_imn,met_imn,
     .                    zn,unim,NIT,NIMs,LABEL59,MMIN,MMAX)
      implicit none
 
      real       pix_imn(_IDIM_,_MDIM_,NIMs)
      real       mod_imn(_IDIM_,_MDIM_,NIMs)
      real*8     wav_imn(_IDIM_,_MDIM_,NIMs)
      integer*1  rej_imn(_IDIM_,_MDIM_,NIMs)
      integer*1  tel_imn(_IDIM_,_MDIM_,NIMs)
      integer*1  met_imn(_IDIM_,_MDIM_,NIMs)
      real       nrm_imn(_IDIM_,_MDIM_,NIMs)
      real       bar_imn(_IDIM_,_MDIM_,NIMs)
      real*8     zn(NIMs)
      logical    unim(NIMs)
 
      real*8  ww
      integer NIT, NIMs, NITu
      character*59 LABEL59
      integer MMIN, MMAX
 
      integer m, i, n
      character*80 FILENAME

      integer OUTPUT_MOD
      common /OUTPUT_MOD_/OUTPUT_MOD

      if (OUTPUT_MOD.eq.0) return

      NITu = NIT
      if (OUTPUT_MOD.eq.1) NITu = 99
      write(FILENAME,'(''.'',i2.2,''.mod'')') NITu
      FILENAME = TRIM(LABEL59) // FILENAME
 
      write(6,'(''SAVE_MOD: '',80a)') FILENAME
 
      open(10,file=FILENAME,status='unknown')
      do n = 1, NIMs
         if (unim(n)) then
            do m = MMIN, MMAX
               do i = 0001, _IDIM_
                  write(10,110) n,m,i,pix_imn(i,m,n),mod_imn(i,m,n),
     .                                wav_imn(i,m,n),
     .                                wav_imn(i,m,n)/(1+zn(n)),
     .                                nrm_imn(i,m,n),
     .                                bar_imn(i,m,n),
     .                                rej_imn(i,m,n),
     .                                tel_imn(i,m,n),
     .                                met_imn(i,m,n)
  110             format(1x,i3.3,1x,i2.2,1x,i4.4,1x,f8.1,1x,f8.1,1x,
     .                      f11.6,1x,f11.6,1x,f8.6,1x,f12.4,1x,
     .                      i1,1x,i1,1x,i1)
                  enddo
               enddo
            endif
         enddo
      close(10)
 
      return
      end
 
 
c-------------------------------------------------
c
c
c
      subroutine save_vel(zn,NIT,NIMs,z0,zbarn,zsign,unim,LABEL59)
      implicit none

      real*8 zn(NIMs)
      integer NIT
      integer NIMs
      real*8  z0(NIMs)
      real*8  zbarn(NIMs)
      real*8  zsign(NIMs)
      logical unim(NIMs)
      character*59 LABEL59
 
      integer n
      character*80 FILENAME
      real*8 zbar

      integer NITu
      integer OUTPUT_VEL
      common /OUTPUT_VEL_/OUTPUT_VEL

      if (OUTPUT_VEL.eq.0) return

      NITu = NIT
      if (OUTPUT_VEL.eq.1) NITu = 99

 
      zbar = 0.
      do n = 1, NIMs
         zbar = zbar + (zn(n)-z0(n))
         enddo
      zbar = zbar/NIMs
c     do n = 1, NIMs
c        zn(n)  = zn(n)  - zbar
c        zbarn(n) = zbarn(n) - zbar
c        enddo
 
      write(6,'(''          '')')
      write(6,'(''          '')')
      write(6,'(''   ZBAR:  '',f12.10)') zbar
      write(6,'(''          '')')
      write(6,'(''          '')')


      write(FILENAME,'(''.'',i2.2,''.vel'')') NITu
      FILENAME = TRIM(LABEL59) // FILENAME

      write(6,'(''SAVE_VEL: '',80a)') FILENAME
      write(6,'(''#   ZBAR: '',f12.10)') zbar
      open(10,file=FILENAME,status='unknown')
      do n = 1, NIMs
         if (unim(n)) then
            write(10,110) n,zn(n),z0(n),
     .                      zbarn(n),zsign(n)
  110       format(1x,i3.3,1x,f12.10,1x,f12.10,
     .                     1x,f12.10,1x,f12.10)
            endif
         enddo
      close(10)
 
      return
      end
 
 
 
c-----------------------------------------------------------
c
c
c
      subroutine save_tem(tems,NIT,temz,LABEL59)
      implicit none
 
      real tems(_SDIM_)
      integer NIT
      real temz(_SDIM_)
      character*59 LABEL59
 
      integer      s
      real*8       ww_s
      character*80 FILENAME
      real         sunw
      real         temp

      integer NITu
      integer OUTPUT_TEM
      common /OUTPUT_TEM_/OUTPUT_TEM

      if (OUTPUT_TEM.eq.0) return

      NITu = NIT
      if (OUTPUT_TEM.eq.1) NITu = 99
      write(FILENAME,'(''.'',i2.2,''.tem'')') NITu
      FILENAME = TRIM(LABEL59) // FILENAME
 
      temp = sunw(5500.0d0)
 


      write(6,'(''SAVE_TEM: '',80a)') FILENAME
      open(10,file=FILENAME,status='unknown')
      do s = 00001, _SDIM_
         write(10,110) s,ww_s(s),tems(s),temz(s),
     .                 sunw(ww_s(s))
  110    format(1x,i5.5,1x,f10.5,1x,f8.6,4x,
     .              f8.6,1x,f8.6,1x,f8.6,1x,f8.6)
         enddo
      close(10)
 
      return
      end
 
 
c-----------------------------------------------------------
c
c
c
      subroutine save_lsf(lsfqkmn,NIT,NIMs,unim,LABEL59,MMIN,MMAX)
      implicit none
 
      real lsfqkmn(_QDIM_,_KDIM_,_MDIM_,NIMs)
      integer NIT
      integer NIMs
      logical unim(NIMs)
      character*59 LABEL59
      integer MMIN, MMAX
 
      integer q,k,m,n
      character*80 FILENAME

      integer NITu
      integer OUTPUT_LSF
      common /OUTPUT_LSF_/OUTPUT_LSF
      NITu = NIT
      if (OUTPUT_LSF.eq.1) NITu = 99
 
      write(FILENAME,'(''.'',i2.2,''.lsf'')') NITu
      FILENAME = TRIM(LABEL59) // FILENAME

      write(*,'(''SAVE: '',80a)') FILENAME
      open(10,file=FILENAME,status='unknown')
      do n = 1, NIMs
         if (unim(n)) then
            do m = MMIN, MMAX
               do k = 001, _WAVN_
                  write(10,110) n,m,k,(lsfqkmn(q,k,m,n),q=1,13)
  110             format(1x,i3.3,1x,i2.2,1x,i1,1x,13(f8.5,1x))
                  enddo
               enddo
            endif
         enddo
      return
      end
 
 
 
c---------------------------------------------------
c
c
c
      subroutine read_mod(pix_imn,mod_imn,wav_imn,
     .                    nrm_imn,bar_imn,rej_imn,
     .                    unim,NIT,NIMs,LABEL59)
      implicit none
 
      real    pix_imn(_IDIM_,_MDIM_,NIMs)
      real    mod_imn(_IDIM_,_MDIM_,NIMs)
      real*8  wav_imn(_IDIM_,_MDIM_,NIMs)
      integer*1  rej_imn(_IDIM_,_MDIM_,NIMs)
      real    nrm_imn(_IDIM_,_MDIM_,NIMs)
      real    bar_imn(_IDIM_,_MDIM_,NIMs)
      logical unim(NIMs)
      integer NIT, NIMs
      character*59 LABEL59

      character*200 STRING
 
      real*8  ww
 
      integer m, i, n, no
      character*80 FILENAME
      real pixu, modu, wavu, wwavu, nrmu, baru
      integer reju
      integer q

      integer OUTPUT_MOD
      common /OUTPUT_MOD_/OUTPUT_MOD
      integer NITu
  
      NITu = NIT
      if (OUTPUT_MOD.eq.1) NITu = 99
      write(FILENAME,'(''.'',i2.2,''.mod'')') NITu
      FILENAME = TRIM(LABEL59) // FILENAME
       
      do n = 1, NIMs
         unim(n) = .false.
         enddo
 
      open(10,file=FILENAME,status='old')
      write(6,'(''READ_MOD: '',80a)') FILENAME
      NIMs = 0
      no = -1
    1 read(10,'(200a)',end=2) STRING
      do q = 1, 200 ! avoid data with NaNs
         if (STRING(q:q).eq.'n'.or.
     .       STRING(q:q).eq.'a'.or.
     .       STRING(q:q).eq.'N') goto 1
         enddo
      read(STRING,*,end=2) n,m,i,pixu,modu,wavu,wwavu,
     .                     nrmu,baru,reju
      pix_imn(i,m,n) = pixu
      mod_imn(i,m,n) = modu
      wav_imn(i,m,n) = wavu
      nrm_imn(i,m,n) = nrmu
      bar_imn(i,m,n) = baru
      rej_imn(i,m,n) = reju
      if (n.gt.NIMs) NIMs = n
      unim(n) = .true.
      if (n.ne.no.and.n.eq.n/10*10) write(6,*) '---> n: ',n
      no = n
      goto 1
    2 close(10)
 
      write(6,*) '           '
      write(6,*) 'IMAGES NOT TO USE...'
      write(6,*) '           '
      no = 0
      do n = 1, NIMs
         if (.not.unim(n)) then
            no = no + 1
            write(6,*) '---> NIM : ',n,no
            endif
         enddo
      write(6,*) '           '
      write(6,*) '---> NIMs: ',no
      write(6,*) '           '
 
      return
      end
 
 
 
c---------------------------------------
c
c
      subroutine read_vel(zn,NIT,NIMs,z0,LABEL59)
      implicit none
      real*8 zn(NIMs)
      integer NIT
      integer NIMs
      real*8 z0(NIMs)
      character*59 LABEL59
 
      integer n, nn
      character*80 FILENAME

      integer NITu
      integer OUTPUT_VEL
      common /OUTPUT_VEL_/OUTPUT_VEL

      NITu = NIT
      if (OUTPUT_VEL.eq.1) NITu = 99
 
      write(FILENAME,'(''.'',i2.2,''.vel'')') NITu
      FILENAME = TRIM(LABEL59) // FILENAME
      write(6,'(''READ_VEL: '',80a)') FILENAME
      open(10,file=FILENAME,status='old')
      do n = 1, NIMs
        read(10,*) nn,zn(n),z0(n)
        if (n.ne.nn) stop
        enddo
 
      return
      end
 
 
c-----------------------------------------------------------
c
c
c
      subroutine read_tem(tems,NIT,LABEL59)
      implicit none
 
      real tems(_SDIM_)
      integer NIT
      character*59 LABEL59
 
      integer      s,ss
      real*8       ww
      character*80 FILENAME

      integer NITu
      integer OUTPUT_TEM
      common /OUTPUT_TEM_/OUTPUT_TEM

      if (OUTPUT_TEM.eq.0) return

      NITu = NIT
      if (OUTPUT_TEM.eq.1) NITu = 99
      write(FILENAME,'(''.'',i2.2,''.tem'')') NITu
      FILENAME = TRIM(LABEL59) // FILENAME
 
      write(6,'(''READ_TEM: '',80a)') FILENAME
      open(10,file=FILENAME,status='old')
      do s = 00001, _SDIM_
         read(10,*) ss,ww,tems(s)
         if (s.ne.ss) stop 'SAVE_TEM:  s.ne.ss'
         enddo
      close(10)
 
      return
      end
 
 
 
c-----------------------------------------------------------
c
c
c
      subroutine read_lsf(lsfqkmn,NIT,NIMs,LABEL59)
      implicit none
 
      real lsfqkmn(_QDIM_,_KDIM_,_MDIM_,NIMs)
      integer NIT, NIMs
      character*59 LABEL59
 
      integer q,k,m,n
      character*80 FILENAME
      integer ku,mu,nu
      real    lsfq(_QDIM_)

      integer NITu
      integer OUTPUT_LSF
      common /OUTPUT_LSF_/OUTPUT_LSF
      NITu = NIT
      if (OUTPUT_LSF.eq.1) NITu = 99
 
      write(FILENAME,'(''.'',i2.2,''.lsf'')') NITu
      FILENAME = TRIM(LABEL59) // FILENAME

      write(6,'(''READ_LSF: '',80a)') FILENAME
      open(10,file=FILENAME,status='old')
      nu = 0
    1 continue
      read(10,*,end=2) n,m,k,(lsfq(q),q=1,_QDIM_)
      do q = 1, _QDIM_
         lsfqkmn(q,k,m,n) = lsfq(q)
         enddo
      if (nu.ne.n.and.n.eq.n/25*25) write(6,*) '  n: ',n
      nu = n
      goto 1
    2 continue
 
      return
      end
 
 
 
 
c-------------------------------------------------------------
c
c this will read in the file for the raw data and will return
c the number of images it will need to read in...
c
      subroutine sense_raw(RJSOURCE,NIMs)
      implicit none
    
      character*80 RJSOURCE
      integer      NIMs

      character*200 STRING
      integer NIM

      open(20,file=RJSOURCE,status='old')
      NIMs = 0
    1 continue
      read(20,'(200a)',end=4) STRING
      if (STRING(1:1).ge.'0'.and.STRING(1:1).le.'9') then
          read(STRING,*) NIM
          if (NIM.gt.NIMs) NIMs = NIM
          endif
      goto 1
    4 continue
      close(20)

      return
      end

 
c-------------------------------------------------------------
c
c this will read in the data, using the input file provided
c on the command line
c
      subroutine read_raw(RJSOURCE,
     .                    pix_imn,nrm_imn,
     .                    tem_imn,uuu_imn,met_imn,
     .                    wav_imn,unim,z0,NIMs,MMIN,MMAX)
      implicit none
 
      character*80 RJSOURCE
      real      pix_imn(_IDIM_,_MDIM_,NIMs)
      real      nrm_imn(_IDIM_,_MDIM_,NIMs)
      real      tem_imn(_IDIM_,_MDIM_,NIMs)
      integer*1 uuu_imn(_IDIM_,_MDIM_,NIMs)
      integer*1 met_imn(_IDIM_,_MDIM_,NIMs)
      real*8 wav_imn(_IDIM_,_MDIM_,NIMs)
      logical unim(NIMs)
      real*8  z0(NIMs)
      integer NIMs
      integer MMIN, MMAX
 
      real*8 pixarray(4021,16)
 
      character*200 STRING
      character*200 FILEI
      character*200 FILEU
      real*8 d, barycorr, vorb, minBC, maxBC
      real*8 velsys, velorb, velear, velsum         !velocities, m/s
      real*8 betsys, betorb, betear, betsum         !beta=vel/c
 
      character*200 RJDIR
      integer nrjdir, i1, i2
 
      integer i, m, n
      character*200 METEORFILE

      logical DOVORB
      common /DOVORB_/DOVORB
 
      do n = 1, NIMs
         unim(n) = .false.
         do m = 1, _MDIM_
         do i = 1, _IDIM_
            pix_imn(i,m,n) = 0.0
            nrm_imn(i,m,n) = 1.
            tem_imn(i,m,n) = 1.
            uuu_imn(i,m,n) = 9
            met_imn(i,m,n) = 0
            enddo
            enddo
         enddo
 
      RJDIR  = './'
      nrjdir = 2

      METEORFILE = 'NONE'
 
      minBC =  99999.0
      maxBC = -99999.0
      NIMs = 0
      write(6,'(''READ_RAW: '',200a)') RJSOURCE
      open(20,file=RJSOURCE,status='old')
      read(20,'(200a)',end=4) STRING
      if (STRING(1:8).ne.'VSYST = ') then
         write(6,*) '                                             '
         write(6,*) 'The first line in the RJ file                '
         write(6,*) '                                             '
         write(6,*) 'has to be VSYST = x , where x is             '
         write(6,*) 'the systemic velocity in m/s.                '
         write(6,*) '                                             '
         write(6,*) '   eg:                                       '
         write(6,*) '                                             '
         write(6,*) '  VSYST = -16400 m/s                         '
         write(6,*) '  00001 "SIGDRA/rj01.46" 185144 4199.221 ... '
         write(6,*) '  00002 "SIGDRA/rj01.74" 185144 4180.827 ... '
         write(6,*) '    .                                        '
         write(6,*) '    .                                        '
         write(6,*) '    .                                        '
         write(6,*) '                                             '
         stop
         endif
      read(STRING(8:18),*) velsys
      write(6,'(''---> SYSTEMIC VELOCITY: '',f10.2,'' m/s'')') velsys
   3  continue
      read(20,'(200a)',end=4) STRING
      if (STRING(1:1).eq.'#') goto 3
      if (STRING(1:6).eq.'METEOR') then
         read(STRING(8:200),*) METEORFILE
         goto 3
         endif
      if (STRING(1:5).eq.'RJDIR') then
         i1 = 0
         i2 = 0
         write(6,*) 'STRING: ',STRING
         do n = 1, 200
            if (STRING(n:n).eq.'"'.and.i1.ne.0
     .                            .and.i2.eq.0) i2 = n-1
            if (STRING(n:n).eq.'"'.and.i1.eq.0) i1 = n+1
           !write(6,*) n,STRING(n:n),i1,i2
            enddo
         if (i1.eq.0.or.i2.eq.0) then
            write(6,*) 'tried to assign RJDIR, but failed.'
            stop
            endif
         RJDIR  = STRING(i1:i2)
         nrjdir = i2-i1+1
         write(6,'(''--->             RJDIR: '',200a)') RJDIR(1:nrjdir)
         write(6,'(''--->            NRJDIR: '',  i2)') nrjdir
         goto 3
         endif
      if (DOVORB) then
         read(STRING,*,end=4) n, FILEI,d,barycorr,vorb
         else
         read(STRING,*,end=4) n, FILEI,d,barycorr
         endif
      if (n.eq.0) goto 4
      if (n.gt.NIMs) NIMs = n
      if (NIMs.gt.(_NDIM_)) then
         write(6,*) '                                      '
         write(6,*) ' The routine is set up to allow only  '
         write(6,*) ' _NDIM_ = ',_NDIM_,' exposures.       '
         write(6,*) ' The RJ file has too many.            '
         write(6,*) '                                      '
         write(6,*) ' You need to recompile after changing '
         write(6,*) ' the _NDIM_ parameter in the defaults.'
         write(6,*) '                                      '
         stop
         endif
      if (barycorr.lt.minBC) minBC = barycorr
      if (barycorr.gt.maxBC) maxBC = barycorr

c---------------------------------------------------
c 
c The 'observation list' input file contains a barycentric correction
c (BARYCORR) for each observation. BARYCORR is the quantity that must
c be added (relativistically) to the observed radial velocity to remove
c the effect of the Earth's motion around the Solar System barycenter.
c Values of BARYCORR in the observation list file come from from
c kbcvel.ascii and are calculated by the observers using only the
c target coordinates and the photon-weighted mean observing time. 
c
c To allow consistent signs in velocity addition formulae, we define
c VELEAR as the component of observed radial velocity due to the
c Earth's motion around the Solar System barycenter. This is minus the
c definition of BARYCORR, so that VELEAR = -BARYCORR.
c
      velear = -barycorr

c---------------------------------------------------
c
c VELORB is the component of observed radial velocity due to orbital
c motion of the source around the barycenter of the distant system.
c
c The 'null hypothesis' ignores orbital motion in the distant system
c when solving for the intrinsic stellar spectrum (aka the template).
c
c Adopting the null hypothesis is equivalent to setting VELORB = 0.0.
c
      velorb = 0.0
      if (DOVORB) velorb = vorb

c---------------------------------------------------
c
c Sum velocity components relativistically and calculate redshift Z.
c Let uppercase U, V, W be velocities and T be the relativistic sum.
c Define lowercase to be beta: u=U/c, v=V/c, w=W/c, and t=T/c.
c
c   T = (U+V+W + UVW/c^2) / (1 + (UV+UW+VW)/c^2)
c   t = (t+v+w + uvw) / (1 + (uv+uw+vw))
c
c   z = sqrt((1+T/c)/(1-T/c)) - 1 = sqrt((1+t)/(1-t)) - 1
c
c To avoid loss of precision, divide by the speed of light (_CMPS_)
c early in the calculation, obtaining beta = v / c for each component.
c
      betsys = velsys / _CMPS_
      betorb = velorb / _CMPS_
      betear = velear / _CMPS_
      betsum = (betsys+betorb+betear + betsys*betorb*betear)
     .       / (1 + betsys*betorb + betsys*betear + betorb*betear)
      z0(n) = dsqrt((1+betsum)/(1-betsum)) - 1.0d0
      velsum = betsum * _CMPS_

c     write(19,119) n, z0(n), barycorr, velsum
c 119 format(i4,1x,f16.12,1x,f16.4,1x,f16.12)
 
c---------------------------------------------------
c
c read in the data from a file; JJ sometimes used DP
c
      FILEU = RJDIR(1:nrjdir) // FILEI
      call rdsk(FILEU,pixarray)    ! single precision
c     call rdsk_DP(FILEU,pixarray) ! double precision
 
c---------------------------------------------------
c
c  go through and fill the arrays with data
c
      do m = MMIN, MMAX
         do i = 1, 4021
            pix_imn(i,m,n) = pixarray(i,m)
            nrm_imn(i,m,n) = 1.
            tem_imn(i,m,n) = 1.
            uuu_imn(i,m,n) = 0
            met_imn(i,m,n) = 0
            enddo
         enddo

      call block_meteor(FILEU,met_imn(1,1,n),METEORFILE)

      unim(n) = .true.  ! use this image
      call init_wavim(wav_imn(1,1,n),MMIN,MMAX) ! load in the initial, default WLS
 
c----------------------------------------------------
c
c report on what was read in...
c
      write(6,'(''RDSK:'',i3.3,1x,a11,1x,15f8.1)')
     .    n, FILEI(1:11), (pix_imn(2000,m,n),m=MMIN,MMAX)
      goto 3
 
   4  close(20)
      write(6,*) ' '

c----------------------------------------------------
c
c write information to .par file
c
      write(60,'(/,a)') '# OBSERVATION LIST FILE'
      write(60,'(a,i0)') 'nobs: ',NIMs
      write(60,'(a,f0.2)') 'velsys: ',velsys
      write(60,'(a,sp,f0.2)') 'minbc: ',minBC
      write(60,'(a,sp,f0.2)') 'maxbc: ',maxbc
      write(60,'(2a)') 'obsdir: ',trim(RJDIR)
      write(60,'(2a)') 'obslist_file: ',trim(RJSOURCE)
      write(60,'(2a)') 'meteor_file: ',trim(METEORFILE)
      call flush(60)

      return
      end
 
 
 
c-------------------------------------------------------
c
c this routine will take a template spectrum and the
c solar spectrum and will shift and scale the solar
c spectrum to make it agree as well as possible with
c the observed one; somewhat crude
c
      subroutine tem_vs_sun(tems,temz)
      implicit none
 
      real tems(_SDIM_) ! this is the extracted template
      real temz(_SDIM_) ! this is nearest solar template
 
      integer s, ss
      integer sa, sb
      integer s1, s2
      real temx(_SDIM_)
      real temy(_SDIM_)
      real*8 ssum, xsum, psum
      real*8 ww_s, ww
      real   sunw, valu, temw
      real*8 velu
 
c
c find s-domain over which the template is defined
c
      sa = 0
      sb = 0
      do s = 1, _SDIM_
         temx(s) = sunw(ww_s(s))
         if (tems(s).lt.0.9) sb = s
         if (tems(s).lt.0.9.and.sa.eq.0) sa = s
         enddo
 
      write(6,'(''sa: '',i6,1x,f8.6,1x,f12.5)')
     .                        sa,tems(sa),ww_s(sa)
      write(6,'(''sb: '',i6,1x,f8.6,1x,f12.5)')
     .                        sb,tems(sb),ww_s(sb)
 
c
c go thru domain 251 pixels at a time; compute ratio
c between sun and star
c
      do s = 00001, _SDIM_
         s1 = max( 0001 ,s-250)
         s2 = min(_SDIM_,s+250)
         ssum = 0.
         xsum = 0.
         do ss = s1, s2             ! compute this average over
            ssum = ssum + tems(ss)  ! 2x the region we'll fill it
            xsum = xsum + temx(ss)  ! over
            enddo
         do ss = s-125, s+125       ! scaled version of sun
            temy(ss) = temx(ss)*ssum/xsum
            enddo
         enddo
 
c
c now, go through and find the ideal offset that will
c match the sun to this spectrum
c
      psum = 0.
      ssum = 0.
      do s = sa+500, sb-500
         ww = ww_s(s)
         valu = (temw(ww*(1+1/3d8),temy)
     .          -temw(ww*(1-1/3d8),temy))/2
         psum = psum + valu*(temy(s)-tems(s))
         ssum = ssum + valu*valu
         enddo
      velu = psum/ssum
 
      write(6,*) 'tem_vs_sun: ',velu
 
c
c now, evaluate the shifted sun, and return it...
c
      if (.false.) then
         open(91,file='INFO_TEMvsSUN.00',status='unknown')
         do s = 0001, _SDIM_
            ww = ww_s(s)
            temz(s) = temw(ww*(1-velu/3d8),temy)
            write(91,191) s,tems(s),temx(s),temy(s),temz(s)
 191        format(i6,1x,8f10.6)
            enddo
         close(91)
         endif
 
      return
      end
 
 
      subroutine find_lines(tems,Ls,sal,sbl,temu,LABEL59)
      implicit none
 
      real    tems(_SDIM_)
      integer Ls
      integer sal(_LDIM_)
      integer sbl(_LDIM_)
      integer temu(_SDIM_)
      character*59 LABEL59
 
      integer L, s, ss
      real*8  ww_s
      character*80 FILENAME

      integer OUTPUT_LIN
      common /OUTPUT_LIN_/OUTPUT_LIN

      if (OUTPUT_LIN.ge.1) then
         FILENAME = TRIM(LABEL59) // '.lines_find'
         open(88,file=FILENAME,status='unknown')
         endif
      Ls = 0
      do s = 1+100, _SDIM_-100
         if (tems(s).gt.0.75) goto 1
         do ss = s-15, s+15
            if (tems(s).gt.tems(ss)) goto 1
            enddo
         Ls = Ls + 1
         sal(Ls) = s-25
         sbl(Ls) = s+25
         if (OUTPUT_LIN.ge.1) 
     .   write(88,188) Ls,sal(Ls),sbl(Ls),s,tems(s)
  188    format(i5,1x,i6,1x,i6,1x,i6,1x,f8.6)
    1    continue
         enddo
      if (OUTPUT_LIN.ge.1) close(88)
      write(6,*) '---> Ls: ',Ls
 
      do L = 1, Ls-1
         if (sbl(L).gt.sal(L+1)) then
            sbl(L  ) = (sbl(L)+sal(L+1))/2 -1
            sal(L+1) = sbl(L) + 2
            endif
         enddo
 
      do s = 1, _SDIM_
         temu(s) = 0
         enddo
 
      do L = 1, Ls
         do s = sal(L), sbl(L)
            temu(s) = L
            enddo
         enddo
 
      if (OUTPUT_LIN.ge.1)  then
         FILENAME = TRIM(LABEL59) // '.lines_id'
         open(89,file=FILENAME,status='unknown')
         do s = 1, _SDIM_
            write(89,189) s, ww_s(s), tems(s), temu(s)
  189       format(1x,i6.6,1x,f9.4,1x,f8.6,1x,i4)
            enddo
         close(89)
         endif
 
      return
      end
 
 
 
c-----------------------------------------------
c
c
c
      subroutine polyfit_w(wa,wb,AA,BB,CC,DD,EE)
      implicit none
 
      real*8 wa(4021)
      real*8 wb(4021)
      real*8 AA,BB,CC,DD,EE
 
      real*8 abcde_i(5,4021)
      common /wav_abcde_comm/abcde_i
      integer i
      real*8  rx
 
      AA = 0.0d0
      BB = 0.0d0
      CC = 0.0d0
      DD = 0.0d0
      EE = 0.0d0
      do i = 0021, 4001
         AA = AA + wa(i)*abcde_i(1,i)/1.0d11
         BB = BB + wa(i)*abcde_i(2,i)/1.0d11
         CC = CC + wa(i)*abcde_i(3,i)/1.0d11
         DD = DD + wa(i)*abcde_i(4,i)/1.0d11
         EE = EE + wa(i)*abcde_i(5,i)/1.0d11
         enddo
 
      do i = 0001, 4021
         rx = (i-2011)/1990.0d0
         wb(i) = AA + BB*rx    + CC*rx**2
     .              + DD*rx**3 + EE*rx**4
         enddo
 
      return
      end
 
 
 
 
      subroutine find_tem_init(pixi,wavi,nrmi,temi,gsni,NTEM)
      implicit none
 
      real    pixi(_IDIM_)
      real*8  wavi(_IDIM_)
      real    nrmi(_IDIM_)
      real    temi(_IDIM_)
      real    gsni(_IDIM_)
      integer NTEM
 
      real*8 normal(25,25)
      real*8 lamron(25,25)
      real*8 vector(_NTEMMX_)
      real*8 coeffs(_NTEMMX_)
      real   rcoeff(_NTEMMX_)
 
      integer n, n1, n2, nc
      real       x1, x2, xi
      integer nn1, nn2
 
      integer o, o1, o2, oj
      real    projn_o(_NTEMMX_,_ODIM_)
 
      integer i, i1, i2
      real    projn_i(_NTEMMX_,_IDIM_)
 
      integer j
 
      real*8  ri, rn
      integer ni, nj, nn
      real    modi(_IDIM_), modu
      integer nn1o(_ODIM_)
      integer nn2o(_ODIM_)
      integer NTEMoo
      common /NTEMo_o/NTEMoo
      data NTEMoo/0/
 
      real rspline, tv
      real dnci(_IDIM_), dnc
      integer ii
      real*8  w
      real    ftso(_ODIM_), ftsw
 
      real    lsfji(_JDIM_,_IDIM_)
      real    lsfj_ga
 
      do i = 001, _IDIM_
         do j = 001, 201
            lsfji(j,i) = lsfj_ga(j,gsni(i),0.)
            enddo
         enddo
 
      do o = 00001, _ODIM_
         ri = 1 + (o-101)*0.1
         ii = int(ri)
         if (ii.lt._IBEG_) ii = _IBEG_
         if (ii.gt._IEND_) ii = _IEND_
         w = wavi(ii) + (ri-ii)*(wavi(ii+1)-wavi(ii))
         ftso(o) = ftsw(w)
         enddo
 
      do i = 1, _IDIM_
         dnci(i) = 9999
         enddo
 
      if (NTEM.ne.NTEMoo) then
         do o = 00001, _ODIM_
            do n = 0001, _NTEMMX_
               projn_o(n,o) = 0.0
               enddo
            ri = 1.0 + (o-101)*0.1
            rn = 1.0 + (ri-_IBEG_)/(_IEND_-_IBEG_)*(NTEM-1)
            if (rn.lt.0001) rn = 0001
            if (rn.gt.NTEM) rn = NTEM
            call rspline_fp(rn,NTEM,projn_o(1,o),nn1o(o),nn2o(o))
            enddo
         NTEMoo = NTEM
         endif
 
      do n = 0001, NTEM, 15
         n1 = n
         n2 = n+24
         if (n2.gt.NTEM) then
            n1 = NTEM-24
            n2 = NTEM
            endif
         nc = (n1+n2)/2
         x1 = _IBEG_ + (n1-1.0)/(NTEM-1)*(_IEND_-_IBEG_)
         x2 = _IBEG_ + (n2-1.0)/(NTEM-1)*(_IEND_-_IBEG_)
         i1 = max(_IBEG_,int(x1-0.5))
         i2 = min(_IEND_,int(x2+0.5))
 
         do i = i1, i2
            o = 101 + (i-1)*10
            do nn = max(0001,n1-1), min(NTEM,n2+1)
               projn_i(nn,i) = 0.
               do j = -100, 100
                  oj = o + j
                  projn_i(nn,i) =
     .            projn_i(nn,i) + ftso(oj)*nrmi(i)*
     .                            lsfji(101+j,i)*projn_o(nn,oj)
                  enddo
               enddo
            if (n1.ne.   1) projn_i(n1,i)=projn_i(n1,i)+projn_i(n1-1,i)
            if (n2.ne.NTEM) projn_i(n2,i)=projn_i(n2,i)+projn_i(n2+1,i)
           enddo
 
         do ni = n1, n2
         do nj = n1, n2
            normal(ni-n1+1,nj-n1+1) = 0.0
            do i = i1, i2
               normal(ni-n1+1,nj-n1+1) =
     .         normal(ni-n1+1,nj-n1+1) + projn_i(ni,i)
     .                                  *projn_i(nj,i)
               enddo
            enddo
            enddo
 
         call dinvertNxN(normal,lamron,25,25)
 
         do nn = n1, n2
            vector(nn) = 0.
            do i = i1, i2
               vector(nn) = vector(nn) + pixi(i)*projn_i(nn,i)
               enddo
            enddo
 
         do nn = n1, n2
            coeffs(nn) = 0.
            do nj = n1, n2
               coeffs(nn) = coeffs(nn) + lamron(nn-n1+1,nj-n1+1)
     .                                  *vector(nj)
               enddo
            xi = _IBEG_ + (nn-1.0)/(NTEM-1)*(_IEND_-_IBEG_)
            rcoeff(nn) = coeffs(nn)
            enddo
 
 
         do i = i1, i2
            modu = 0.
            do ni = n1, n2
               modu = modu + coeffs(ni)*projn_i(ni,i)
               enddo
            ri = i
            rn = 1.0 + (ri-_IBEG_)/(_IEND_-_IBEG_)*(NTEM-1)
            tv = rspline(rn-n1+1,rcoeff(n1),n2-n1+1)
            o  = 101 + (i-1)*10
            dnc = abs(rn-nc)
            tv = rspline(rn-n1+1,rcoeff(n1),n2-n1+1)
            if (dnc.lt.dnci(i)) then
               dnci(i) = dnc
               temi(i) = tv
               modi(i) = modu
               endif
            enddo
 
         enddo
 
      return
      end
 
 
 
 
 
      subroutine find_tem_sngl(pixi,reji,wavi,modi,
     .                              nrmi,temi,lsfji,NTEM)
      implicit none
 
      real    pixi(_IDIM_)
      integer*1  reji(_IDIM_)
      real*8  wavi(_IDIM_)
      real    modi(_IDIM_)
      real    nrmi(_IDIM_)
      real    temi(_IDIM_)
      real    lsfji(_JDIM_,_IDIM_)
      integer NTEM
 
      real*8 normal(25,25)
      real*8 lamron(25,25)
      real*8 vector(_NTEMMX_)
      real*8 coeffs(_NTEMMX_)
      real   rcoeff(_NTEMMX_)
 
      integer n, n1, n2, nc
      real       x1, x2, xi
      integer nn1, nn2
 
      integer o, o1, o2, oj
      real    projn_o(_NTEMMX_,_ODIM_)
 
      integer i, i1, i2
      real    projn_i(_NTEMMX_,_IDIM_)
 
      integer j
 
      real*8  ri, rn
      integer ni, nj, nn
      integer nn1o(_ODIM_)
      integer nn2o(_ODIM_)
      integer NTEMo
      common /NTEMo_/NTEMo
      data NTEMo/0/
 
      real rspline, tv
      real dnci(_IDIM_), dnc
      integer ii
      real*8  w
      real    ftso(_ODIM_), ftsw
      real    nrm
 
 
      if (NTEM.ne.NTEMo) then
         do o = 00001, _ODIM_
            do n = 0001, _NTEMMX_
               projn_o(n,o) = 0.0
               enddo
            ri = 1.0 + (o-101)*0.1
            rn = 1.0 + (ri-_IBEG_)/(_IEND_-_IBEG_)*(NTEM-1)
            if (rn.lt.0001) rn = 0001
            if (rn.gt.NTEM) rn = NTEM
            call rspline_fp(rn,NTEM,projn_o(1,o),nn1o(o),nn2o(o))
            enddo
         NTEMo = NTEM
         endif
 
 
 
      do o = 00001, _ODIM_
         ri = 1 + (o-101)*0.1
         ii = int(ri)
         if (ii.lt._IBEG_) ii = _IBEG_
         if (ii.gt._IEND_) ii = _IEND_
         w   = wavi(ii) + (ri-ii)*(wavi(ii+1)-wavi(ii))
         nrm = nrmi(ii) + (ri-ii)*(nrmi(ii+1)-nrmi(ii))
         ftso(o) = ftsw(w)*nrm
         enddo
 
      do i = 1, _IDIM_
         dnci(i) = 9999
         enddo
 
      do n = 0001, NTEM, 15
         n1 = n
         n2 = n+24
         if (n2.gt.NTEM) then
            n1 = NTEM-24
            n2 = NTEM
            endif
         nc = (n1+n2)/2
         x1 = _IBEG_ + (n1-1.0)/(NTEM-1)*(_IEND_-_IBEG_)
         x2 = _IBEG_ + (n2-1.0)/(NTEM-1)*(_IEND_-_IBEG_)
         i1 = max(_IBEG_,int(x1-0.5))
         i2 = min(_IEND_,int(x2+0.5))
 
         do i = i1, i2
            o = 101 + (i-1)*10
            do nn = max(0001,n1-1), min(NTEM,n2+1)
               projn_i(nn,i) = 0.
               do j = -100, 100
                  oj = o + j
                  projn_i(nn,i) =
     .            projn_i(nn,i) + ftso(oj)*
     .                            lsfji(101+j,i)*projn_o(nn,oj)
                  enddo
               enddo
            if (n1.ne.   1) projn_i(n1,i)=projn_i(n1,i)+projn_i(n1-1,i)
            if (n2.ne.NTEM) projn_i(n2,i)=projn_i(n2,i)+projn_i(n2+1,i)
            enddo
 
         do ni = n1, n2
         do nj = n1, n2
            normal(ni-n1+1,nj-n1+1) = 0.0
            do i = i1, i2
               if (reji(i).ne.1)
     .         normal(ni-n1+1,nj-n1+1) =
     .         normal(ni-n1+1,nj-n1+1) + projn_i(ni,i)
     .                                  *projn_i(nj,i)
               enddo
            enddo
            enddo
 
         call dinvertNxN(normal,lamron,25,25)
 
         do nn = n1, n2
            vector(nn) = 0.
            do i = i1, i2
               if (reji(i).ne.1)
     .         vector(nn) = vector(nn) + pixi(i)*projn_i(nn,i)
               enddo
            enddo
 
         do nn = n1, n2
            coeffs(nn) = 0.
            do nj = n1, n2
               coeffs(nn) = coeffs(nn) + lamron(nn-n1+1,nj-n1+1)
     .                                  *vector(nj)
               enddo
            xi = _IBEG_ + (nn-1.0)/(NTEM-1)*(_IEND_-_IBEG_)
            rcoeff(nn) = coeffs(nn)
            enddo
 
 
         do i = i1, i2
            ri = i
            rn = 1.0 + (ri-_IBEG_)/(_IEND_-_IBEG_)*(NTEM-1)
            tv = rspline(rn-n1+1,rcoeff(n1),n2-n1+1)
            o  = 101 + (i-1)*10
            dnc = abs(rn-nc)
            tv = rspline(rn-n1+1,rcoeff(n1),n2-n1+1)
            if (dnc.lt.dnci(i)) then
               dnci(i) = dnc
               temi(i) = tv
               endif
            enddo
 
         enddo
 
      return
      end
 
 
c--------------------------------------------------------
c
c
c
      subroutine find_tem_snglQ(pixi,reji,modi,
     .                          nrmi,temi)
      implicit none
 
      real    pixi(_IDIM_)
      integer*1  reji(_IDIM_)
      real    modi(_IDIM_)
      real    nrmi(_IDIM_)
      real    temi(_IDIM_)
 
      integer N, Ns
      integer i, ii
      real*8 psum, msum
      real   ff
 
      real    tq(_IDIM_)
      real    rspline
      real*8  rn
 
 
      Ns = 0
      do i = 0001 + 2, _IDIM_ -2, 4
         psum = 0.
         msum = 0.
         do ii = -2, 2, 1
            ff = 0.250
            if (ii.eq.-2) ff = 0.125
            if (ii.eq. 2) ff = 0.125
            psum = psum + ff*pixi(i+ii)
            msum = msum + ff*temi(i+ii)*modi(i+ii)*nrmi(i+ii)
            enddo
         Ns = Ns + 1
         tq(Ns) = psum/msum
         enddo
 
      do i = 0001, _IDIM_
         rn = 1 + (i-3)*0.25
         temi(i) = temi(i)*rspline(rn,tq,Ns)
         if (temi(i).lt.0.05) temi(i) = 0.05
         if (temi(i).gt.1.00) temi(i) = 1.00
         enddo
 
      return
      end
 
 
 
c-----------------------------------------------
c
c
c
c
c-----------------------------------------------
 
      subroutine tweak_wls(pixi,reji,wavi,nrmi,temi,lsfji,coeffw)
      implicit none
 
      real   pixi(_IDIM_)
      integer*1 reji(_IDIM_)
      real*8 wavi(_IDIM_)
      real   nrmi(_IDIM_)
      real   temi(_IDIM_)
      real   lsfji(_JDIM_,_IDIM_)
      real   coeffw(5)
 
      real   modo(_ODIM_)
      real   fxto(_ODIM_)
 
      real*8 wavh(_IDIM_)
 
      real   modx(_IDIM_)
      real   modi(_IDIM_)
      real   modl(_IDIM_), modlu
      integer o
      real    ri, t, f, ftsw
      integer ii, j, i, ic
      real*8  w
      real    msuma
      real    msumb
      real    projn_i(5,_IDIM_)
      real*8  rn
      integer n, n1, n2
      integer    ni, nj
      integer    i1, i2
      real nrm
 
      real*8  normal(5,5)
      real*8  lamron(5,5)
      real*8  vector(5)
      real*8  coeffs(5)
 
      do o = 00001, _ODIM_
         ri = 1 + (o-101)*0.1
         ii = int(ri)
         if (ii.lt._IBEG_) ii = _IBEG_
         if (ii.gt._IEND_) ii = _IEND_
         w   = wavi(ii) + (ri-ii)*(wavi(ii+1)-wavi(ii))
         t   = temi(ii) + (ri-ii)*(temi(ii+1)-temi(ii))
         nrm = nrmi(ii) + (ri-ii)*(nrmi(ii+1)-nrmi(ii))
         f = ftsw(w)
         fxto(o) = f*t*nrm
         enddo
 
 
      do o = 00001, _ODIM_
         i = 1 + (o-101)*0.1
         if (i.le._IBEG_) i = _IBEG_
         if (i.ge._IEND_) i = _IEND_
         modo(o) = 0.
         do j = -100, 100
            if (o+j.gt.00001.and.o+j.lt._ODIM_)
     .      modo(o) = modo(o) + fxto(o+j)*lsfji(101+j,i)
            enddo
         enddo
 
      do i = 0001, _IDIM_
         o = 101 + (i-1)*10
         modx(i) = (modo(o+1) - modo(o-1))/2
         modi(i) = modo(o)
         rn = 1 + 4.*(i-_RATB_)/(_RATY_-_RATB_)
         call rspline_fp(rn,5,projn_i(1,i),n1,n2)
         enddo
 
      do ni = 1, 5
      do nj = 1, 5
         normal(ni,nj) = 0.
         do i = _RATB_, _RATY_
            if (reji(i).eq.0)
     .      normal(ni,nj) = normal(ni,nj) + projn_i(ni,i)*modx(i)
     .                                     *projn_i(nj,i)*modx(i)
            enddo
         enddo
         enddo
 
      call dinvertNxN(normal,lamron,5,5)
 
      do n = 1, 5
         vector(n) = 0
         do i = _RATB_, _RATY_
            if (reji(i).eq.0)
     .      vector(n) =
     .      vector(n) + (pixi(i)-modi(i))*projn_i(n,i)*modx(i)
            enddo
         enddo
 
      do n = 1, 5
         coeffs(n) = 0
         do nj = 1, 5
            coeffs(n) = coeffs(n) + lamron(n,nj)*vector(nj)
            enddo
         enddo
 
      do i = 0001, _IDIM_
         modl(i) = 0.
         do n = 1, 5
            modl(i) = modl(i) + coeffs(n)*projn_i(n,i)
            enddo
         enddo
 
      do i = 0001, _IBEG_
         wavh(i) = wavi(_IBEG_) + (i-_IBEG_)*(wavi(_IBEG_+1)-
     .                                        wavi(_IBEG_  ))
         enddo
      do i = _IEND_, _IDIM_
         wavh(i) = wavi(_IEND_) + (i-_IEND_)*(wavi(_IEND_  )-
     .                                        wavi(_IEND_-1))
         enddo
 
      do i = _IBEG_, _IEND_
         ic = i
         if (ic.lt._IBEG_+1) ic = _IBEG_+1
         if (ic.gt._IEND_-1) ic = _IEND_-1
         i1 = ic-1
         i2 = ic+1
         modlu = modl(i)
         if (modlu.gt. 1) modlu =  1
         if (modlu.lt.-1) modlu = -1
         wavh(i) = wavi(i) + modlu/06/(i2-i1)*(wavi(i2)-wavi(i1))
         enddo
 
      do i = 0001, _IDIM_
         wavi(i) = wavh(i)
         enddo
 
      do n = 1, 5
         coeffw(n) = coeffs(n)
         enddo
 
      return
      end
 
 
c-----------------------------------------------
c
c
c
c
c-----------------------------------------------
 
      subroutine tweak_gsn(pixi,reji,wavi,nrmi,gsni,lsfji,coeffg)
      implicit none
 
      real   pixi(_IDIM_)
      integer*1 reji(_IDIM_)
      real*8 wavi(_IDIM_)
      real   nrmi(_IDIM_)
      real   gsni(_IDIM_)
      real   lsfji(_JDIM_,_IDIM_)
      real   coeffg(_WAVN_)
 
      real   modo(_ODIM_)
      real   modd(_ODIM_)
      real   fxto(_ODIM_)
 
      real*8 wavh(_IDIM_)
 
      real   modx(_IDIM_)
      real   modi(_IDIM_)
      real   modl(_IDIM_)
      integer o
      real    ri, t, f, ftsw
      integer ii, j, i
      real*8  w
      real    msuma
      real    msumb
      real    projn_i(_WAVN_,_IDIM_)
      real*8  rn
      integer n, n1, n2
      integer    ni, nj
      integer    i1, i2
      real       rg
      integer    ig
      real       dgsni
 
      real*8  normal(_WAVN_,_WAVN_)
      real*8  lamron(_WAVN_,_WAVN_)
      real*8  vector(_WAVN_)
      real*8  coeffs(_WAVN_)
 
      real  dlsfji(_JDIM_,_IDIM_)
      real  lsfp, lsfm
      real  lsfj_ga
 
 
      do o = 00001, _ODIM_
         ri = 1 + (o-101)*0.1
         ii = int(ri)
         if (ii.lt._IBEG_) ii = _IBEG_
         if (ii.gt._IEND_) ii = _IEND_
         w = wavi(ii) + (ri-ii)*(wavi(ii+1)-wavi(ii))
         t = nrmi(ii) + (ri-ii)*(nrmi(ii+1)-nrmi(ii))
         f = ftsw(w)
         fxto(o) = f*t
         enddo
 
 
      do i = 0001, _IDIM_
         rg = gsni(i)
         do j = 001, _JDIM_
            lsfji(j,i)  = lsfj_ga(j,rg    ,0.0)
            lsfp        = lsfj_ga(j,rg+1.0,0.0)
            lsfm        = lsfj_ga(j,rg-1.0,0.0)
            dlsfji(j,i) = (lsfp-lsfm)/2.0
            enddo
         enddo
 
 
      do o = 00001, _ODIM_
         i = 1 + (o-101)*0.1
         if (i.le._IBEG_) i = _IBEG_
         if (i.ge._IEND_) i = _IEND_
         modo(o) = 0.
         modd(o) = 0.
         do j = -100, 100
            if (o+j.gt.00001.and.o+j.lt._ODIM_) then
               modo(o) = modo(o) + fxto(o+j)* lsfji(101+j,i)
               modd(o) = modd(o) + fxto(o+j)*dlsfji(101+j,i)
               endif
            enddo
         enddo
 
 
      do i = 0001, _IDIM_
         o = 101 + (i-1)*10
         modx(i) = modd(o)
         modi(i) = modo(o)
         rn = 1 + (_WAVN_-1.0)*(i-_WAVB_)/(_WAVY_-_WAVB_)
         call rspline_fp(rn,_WAVN_,projn_i(1,i),n1,n2)
         enddo
 
 
      do ni = 1, _WAVN_
      do nj = 1, _WAVN_
         normal(ni,nj) = 0.
         do i = _WAVB_, _WAVY_
            if (reji(i).eq.0)
     .      normal(ni,nj) = normal(ni,nj) + projn_i(ni,i)*modx(i)
     .                                     *projn_i(nj,i)*modx(i)
            enddo
         enddo
         enddo
 
 
      call dinvertNxN(normal,lamron,_WAVN_,_WAVN_)
 
 
      do n = 1, _WAVN_
         vector(n) = 0
         do i = _WAVB_, _WAVY_
            if (reji(i).eq.0)
     .      vector(n) =
     .      vector(n) + (pixi(i)-modi(i))*projn_i(n,i)*modx(i)
            enddo
         enddo
 
      do n = 1, _WAVN_
         coeffs(n) = 0
         do nj = 1, _WAVN_
            coeffs(n) = coeffs(n) + lamron(n,nj)*vector(nj)
            enddo
         enddo
 
 
      do i = 0001, _IDIM_
         do n = 1, _WAVN_
            dgsni = coeffs(n)*projn_i(n,i)
            if (dgsni.gt. 1) dgsni =  1
            if (dgsni.lt.-1) dgsni = -1
            gsni(i) = gsni(i) + dgsni*0.90
            enddo
         enddo
 
      do i = 0001, _IDIM_
         do j = 001, _JDIM_
            lsfji(j,i)  = lsfj_ga(j,gsni(i),0.0)
            enddo
         enddo
 
c     write(*,'(''tweak_gsn: '',40x,99f7.3)') (coeffs(n),n=1,_WAVN_)
 
      do n = 1, _WAVN_
         i = _WAVB_ + (_WAVY_-_WAVB_)*(n-1)/(_WAVN_-1.0)
         coeffg(n) = gsni(i)
         enddo
 
      return
      end
 
c#include "51.ROUTINES/find_vel4sunX.f"
 
 
 
c------------------------------------------------------
c
c
c
c
c
      subroutine vlinefit_LSQ(pixim,wavim, nrmim,rejim,
     .                        tems,vel,lsfjim,
     .                        Ls,sal,sbl,NIM,
     .                        vl,tl,il,ml,MMIN,MMAX)
      implicit none
 
      real    pixim(_IDIM_,_MDIM_)
      real*8  wavim(_IDIM_,_MDIM_)
      real    nrmim(_IDIM_,_MDIM_)
      integer*1  rejim(_IDIM_,_MDIM_)
      real    tems(_SDIM_)
      real*8  vel
      real    lsfjim(_JDIM_,_IDIM_,_MDIM_)
      integer Ls
      integer sal(_LDIM_)
      integer sbl(_LDIM_)
      integer NIM
      real*8  vl(_LDIM_)
      real*8  tl(_LDIM_)
      integer il(_LDIM_)
      integer ml(_LDIM_)
      integer MMIN, MMAX
 
      integer L
      integer s, sbar
      integer m, mmyn
      real*8  w, ww, w1, w2, ww_s
      integer o, o1, o2
      integer i, i1, i2
      integer iu, imin, i_w, ii, j
      real    ftso(_ODIM_), ftsw
      real*8  wavo(_ODIM_)
      real    temo(_ODIM_)
      real    temx(_ODIM_)
      real    temw
      real    modi(_IDIM_)
      real    modx(_IDIM_)
      real*8  psum, msum
      real*8  e, emin, vmin, rmin
      real    tmin
      real*8  v, vv, velu
      real*8  ri
 
 
 
      do L = 2, Ls-1
         tmin = 9e9
         do s = sal(L), sbl(L)
            if (tems(s).lt.tmin) tmin = tems(s)
            enddo
         sbar = (sal(L)+sbl(L))/2
         ww = ww_s(sbar)
         w  = ww*(1+vel)
         mmyn = 0
         imin = 100
         do m = MMIN, MMAX
            iu = i_w(w,wavim(1,m))
            if (abs(iu-2010).lt.abs(imin-2010)) then
               mmyn = m
               imin = iu
               endif
            enddo
         if (mmyn.eq.0) goto 333  !need logic to ignore missing lines
         m  = mmyn
         w1 = ww_s(sal(L))*(1+vel)
         w2 = ww_s(sbl(L))*(1+vel)
         i1 = i_w(w1,wavim(1,m))
         i2 = i_w(w2,wavim(1,m))
         o1 = 101 + (i1-1)*10 - 200
         o2 = 101 + (i2-1)*10 + 200
         do o = o1, o2
            ri = 1 + (o-101)*0.1
            ii = int(ri)
            if (ii.lt.0001) ii = 0001
            if (ii.gt.4020) ii = 4020
            wavo(o) = wavim(ii,m) + (ri-ii)*(wavim(ii+1,m)-wavim(ii,m))
            ftso(o) = ftsw(wavo(o))
            enddo
        velu = vel
        do o = o1, o2
           temo(o) =  temw(wavo(o)/(1+velu),tems)
           temx(o) = (temw(wavo(o)/(1+velu+0.5/3.000d8),tems)
     .               -temw(wavo(o)/(1+velu-0.5/3.000d8),tems))
           enddo
        psum = 0.
        msum = 0.
        do i = i1, i2
           if (rejim(i,m).eq.0) then
              modi(i) = 0.
              modx(i) = 0.
              o = 101 + (i-1)*10
              do j = -100, 100
                 modi(i) = modi(i) + nrmim(i,m)
     .                               *ftso(o+j)
     .                               *temo(o+j)
     .                               *lsfjim(101+j,i,m)
                 modx(i) = modx(i) + nrmim(i,m)
     .                               *ftso(o+j)
     .                               *temx(o+j)
     .                               *lsfjim(101+j,i,m)
                 enddo
              psum = psum + modx(i)*(pixim(i,m)-modi(i))
              msum = msum + modx(i)*modx(i)
              endif
           enddo!i
         v = psum/msum
         if (.not.(v.lt. 999.0)) v =  999.0
         if (.not.(v.gt.-999.0)) v = -999.0
         velu = velu + v/3.000d8
         vl(L) = velu
         tl(L) = tmin
         ml(L) = m
         il(L) = (i1+i2)/2
 333     continue
         enddo
 
      return
      end
 
 
 
c---------------------------------
c
c
      subroutine symm_lsfqk(lsfqk)
      implicit none
 
      real lsfqk(_QDIM_,_KDIM_)
      integer q, k, qq
 
      do k = 1, _KDIM_
      do q = 1, _QDIM_/2
         qq = _QDIM_-(q-1)
         lsfqk( q,k) = (lsfqk(q,k)+lsfqk(qq,k))/2
         lsfqk(qq,k) =  lsfqk(q,k)
         enddo
         enddo
 
      return
      end
 
 
c-----------------------------------
c
c
      subroutine lsfji2lsfqk(lsfji,lsfqk)
      implicit none
 
      real lsfji(_JDIM_,_IDIM_)
      real lsfqk(_QDIM_,_KDIM_)
 
      integer k, i
 
      do k = 1, _KDIM_
         i = _WAVB_ + (k-1.)/(_KDIM_-1)*(_WAVY_-_WAVB_)
         call lsfj2lsfq(lsfji(1,i),lsfqk(1,k))
         enddo
 
      return
      end
 
 
 
c------------------------------------------------------
c
c
c
      subroutine anal_vlines(vln,tln,iln,mln,
     .                       unim,Ns,Ls,NIT,
     .                       zbarn,zsign,
     .                       zbarl,zsigl,LABEL59)
      implicit none
 
      real*8  vln(_LDIM_,Ns)
      real*8  tln(_LDIM_,Ns)
      integer iln(_LDIM_,Ns)
      integer mln(_LDIM_,Ns)
      logical unim(Ns)
      integer Ns, Ls, NIT
      real*8  zbarn(Ns), zsign(Ns)
      real*8  zbarl(_LDIM_), zsigl(_LDIM_)
      character*59 LABEL59
 
      real*8 barn, sign
      real*8 barl, sigl
      real*8 bart, sigt
      real*8 bari, sigi
      real*8 barm, sigm
 
      real*8 xlist(999)
      real*8 bar
      real*8 sig
      integer NTOT, NUSE, Nu
      integer l, m, i, n
      real*8  v, dv, t
      integer Lu, Luu
      character*80 FILENAME
      real*8 vlist(_LDIM_)
      real*8 tlist(_LDIM_)
      real*8 ilist(_LDIM_)
      real*8 mlist(_LDIM_)
 
      write(6,*) '                         '
      write(6,*) 'enter anal_vlines... Ns: ',Ns
      write(6,*) '                     Ls: ',Ls
 
      write(FILENAME,'(''.'',i2.2,''.vln'')') NIT
      FILENAME = TRIM(LABEL59) // FILENAME
      open(41,file=FILENAME,status='unknown')
      do N = 1, Ns
         if (unim(n)) then
            do L = 1, Ls
               write(41,140) L,vln(l,n),tln(l,n),mln(l,n),iln(l,n),n
  140          format(i4,1x,f13.10,1x,1x,f8.5,1x,i2,1x,i4,1x,i3)
               enddo
            endif
         enddo
      close(41)
 
 
 
c----------------------------------------
c
c find initial velocity for each exposure
c
      do n = 1, Ns
         if (unim(n)) then
            Lu = 0
            do l = 2, Ls-1
               if (abs(zbarl(l)).lt.200.and.zsigl(l).lt.100) then
                  Lu = Lu + 1
                  xlist(Lu) = vln(l,n)
                  endif
               enddo
            call barsig(xlist,Lu,barn,sign,Luu)
            zbarn(n) = barn
            zsign(n) = sign
            endif
         enddo
 
 
 
      do l = 1, Ls
         Nu = 0
         do N = 1, Ns
            if (unim(n)) then
               Nu = Nu + 1
               vlist(Nu) = vln(l,n)-zbarn(n)
               tlist(Nu) = tln(l,n)
               ilist(Nu) = iln(l,n)
               mlist(Nu) = mln(l,n)
               endif
            enddo
         call barsig(vlist,Nu,barl,sigl,NUSE)
         call barsig(tlist,Nu,bart,sigt,NUSE)
         call barsig(ilist,Nu,bari,sigi,NUSE)
         call barsig(mlist,Nu,barm,sigm,NUSE)
         zbarl(l) = barl
         zsigl(l) = sigl
         enddo
 
 
 
 
      write(FILENAME,'(''.'',i2.2,''.vsn'')') NIT
      FILENAME = TRIM(LABEL59) // FILENAME
      open(43,file=FILENAME,status='unknown')
 
      do n = 1, Ns
         if (unim(n)) then
            Lu = 0
            do l = 2, Ls-1
               if (abs(zbarl(l))*3e8.lt.200.and.
     .                 zsigl(l) *3e8.lt.100) then
                  Lu = Lu + 1
                  xlist(Lu) = vln(l,n)-zbarl(l)
                  endif
               enddo
            call barsig(xlist,Lu,barn,sign,Luu)
            write(06,143) n,barn,sign,Lu,Luu
            write(43,143) n,barn,sign,Lu,Luu
  143       format(i4,1x,f12.9,1x,f12.9,1x,i4,1x,i4)
            zbarn(n) = barn
            zsign(n) = sign
            endif
         enddo
         close(43)
 
 
      write(FILENAME,'(''.'',i2.2,''.vsl'')') NIT
      FILENAME = TRIM(LABEL59) // FILENAME
      open(42,file=FILENAME,status='unknown')
 
      do l = 1, Ls
         Nu = 0
         do N = 1, Ns
            if (unim(N)) then
               Nu = Nu + 1
               vlist(Nu) = vln(l,n)-zbarn(n)
               tlist(Nu) = tln(l,n)
               ilist(Nu) = iln(l,n)
               mlist(Nu) = mln(l,n)
               endif
            enddo
         call barsig(vlist,Nu,barl,sigl,NUSE)
         call barsig(tlist,Nu,bart,sigt,NUSE)
         call barsig(ilist,Nu,bari,sigi,NUSE)
         call barsig(mlist,Nu,barm,sigm,NUSE)
         zbarl(l) = barl
         zsigl(l) = sigl
         write(42,347) l,barl,sigl,bart,barm,bari,Nu,NUSE
  347    format(i3,1x,f12.9,1x,f12.9,1x,f8.6,1x,2f8.2,1x,2i7)
         enddo
      close(42)
 
      write(6,*) '                         '
      write(6,*) 'exit  anal_vlines...'
      write(6,*) '                         '
      write(6,*) '                         '
 
      end
 
 
 
c-----------------------------------------------
c
c
c
      subroutine barsig(xlist,NTOT,bar,sig,NUSE)
      implicit none
 
      integer NTOT
      real*8 xlist(NTOT)
      real*8 bar
      real*8 sig
      integer NUSE
 
      integer n
      real*8  bsum, ssum
      integer nsum
      integer NIT
 
 
      bar = 0.e0
      sig = 9e9
      do NIT = 1, 30
         bsum = 0.
         ssum = 0.
         nsum = 0.
         do n = 1, NTOT
            if (abs(xlist(n)-bar).le.2.75*sig) then
               bsum = bsum + xlist(n)
               ssum = ssum + (xlist(n)-bar)**2
               nsum = nsum + 1
               endif
            enddo
         if (nsum.gt.0) bar = bsum / nsum
         if (nsum.gt.1) sig = sqrt(ssum/(nsum-1))
         enddo
      NUSE = nsum
      if (nsum.le.1) sig = 0.999
 
      return
      end
 
 
 
c----------------------------------------------
c
c 
      subroutine cnorm_tem(tems,temo,NIT,LABEL59)
      implicit none
 
      real   tems(_SDIM_)
      real   temo(_SDIM_)
      integer     NIT
      character*59 LABEL59

      real   tem2(_SDIM_)
      real   tem3(_SDIM_)
      real   tem95
 
      real   temh(_SDIM_)
      real   teml(_SDIM_)
 
      integer s, ss ,s1, s2
      integer l, Ls
      real    f, ftot, ptot
      character*80 FILENAME

      integer NITu
      integer OUTPUT_NRM
      common /OUTPUT_NRM_/OUTPUT_NRM
 
      write(6,*) 'enter cnorm_tem...'
 
      do s = 0001, _SDIM_
         s1 = max( 0001 ,s-150)
         s2 = min(_SDIM_,s+150)
         Ls = 0
         do ss = s1, s2
            Ls = Ls + 1
            teml(Ls) = tems(ss)
            enddo
         call rbubble(teml,Ls)
         tem95 = teml(int(0.90*Ls))
         ptot = 0.
         ftot = 0.
         do ss = s1, s2
            if (tems(ss).gt.tem95-0.01) then
               ftot = ftot + 1.0
               ptot = ptot + tems(ss)
               endif
            enddo
         tem2(s) = ptot/ftot
         if (s.eq.s/1000*1000)
     .   write(6,111) s,s1,s2,Ls,tems(s),tem2(s)
         enddo

      if (OUTPUT_NRM.ge.0) then 
         NITu = NIT
         if (OUTPUT_NRM.eq.1) NITu = 99
         write(FILENAME,'(''.'',i2.2,''.cnr'')') NITu
         FILENAME = TRIM(LABEL59) // FILENAME
         open(11,file=FILENAME,status='unknown')
         endif
      do s = 0001, _SDIM_
         s1 = max( 0001 ,s-250)
         s2 = min(_SDIM_,s+250)
         ftot = 0.
         ptot = 0.
         do ss = s1, s2
            f = 1.00-abs(ss-s)/200.
            ftot = ftot + f
            ptot = ptot + f*tem2(ss)
            enddo
         tem3(s) = ptot/ftot
         temh(s) = tems(s) - (tem3(s)-1.0000)
         if (OUTPUT_NRM.ge.0) 
     .   write(11,111) s,s1,s2,Ls,tems(s),tem2(s),tem3(s),temh(s)
         if (s.eq.s/1000*1000)
     .   write(06,111) s,s1,s2,Ls,tems(s),tem2(s),tem3(s),temo(s)
  111    format(i6,1x,i6,1x,i6,1x,i3,1x,f8.6,1x,f8.6,1x,f8.6,1x,f8.6)
         enddo
      if (OUTPUT_NRM.ge.0) close(11)
 
      do s = 0001, _SDIM_
         temo(s) = temh(s)
         enddo
 
      return
      end
 
 
 
 
c--------------------------------------
c
c bubble-sorts a real*4 list into ascending order
c
      subroutine rbubble(r1,NTOT)
      implicit none
      integer NTOT
      real r1(NTOT)
      real temp
      integer n
      logical change
 777  continue
      change = .false.
      do n = 1, NTOT-1
         if (r1(n).gt.r1(n+1)) then
            temp = r1(n)
            r1(n) = r1(n+1)
            r1(n+1) = temp
            change = .true.
            endif
         enddo
      if (change) goto 777
      end
 
c--------------------------------------
c
c bubble-sorts a real*8 list into ascending order
c
      subroutine dbubble(r1,NTOT)
      implicit none
      integer NTOT
      real*8 r1(NTOT)
      integer n
      real*8 temp
      logical change
 777  continue
      change = .false.
      do n = 1, NTOT-1
         if (r1(n).gt.r1(n+1)) then
            temp = r1(n)
            r1(n) = r1(n+1)
            r1(n+1) = temp
            change = .true.
            endif
         enddo
      if (change) goto 777
      end
 
 
 
c---------------------------------------------------- 
c 
c this routine will report whether a given wavelength 
c is impacted by telluric lines:  it will report the 
c telluric line depth; 0.0 means none 
c 
c it will tabulate Jeff's list between 5000 and 6000  
c angstroms 
c 
 
      real function tell_depth(w) 
      implicit none 
      real*8 w 
 
      real rvtell_bd(4,406)    
      common /rvtell_/rvtell_bd 
      integer n 
 
      logical first 
      common /td_first/first 
      data first/.true./ 
      integer s 
      real*8 ww 
 
      real tdepth_s(1000001) 
      common /tdepth_s_/tdepth_s 
     
      if (first) then  
         write(6,*) 'SET UP TELLURIC-DEPTH ARRAY...' 
         do s = 0000001, 1000001  
            ww = 5000.0d0 + (s-1)*0.01 
            tdepth_s(s) = 0. 
            do n = 1, 406 
               if (ww.ge.rvtell_bd(1,n).and. 
     .             ww.le.rvtell_bd(3,n)) 
     .             tdepth_s(s) = rvtell_bd(4,n) 
               enddo 
            enddo 
         first = .false. 
         endif 
         
      tell_depth = 0. 
      s = 1 + (w-5000.0d0)*100.0 + 0.5  
      if (s.lt.0000001) return 
      if (s.gt.1000001) return 
      tell_depth = tdepth_s(s) 
  
      return 
      end 
 
 
c---------------------------------------- 
c 
c this is the list that JValenti gave me, 
c from his analysis of the NSO atlas. 
c 
      block data rvtell 
      real rvtell_bd(4,406)    
      common /rvtell_/rvtell_bd 
     
      data rvtell_bd / 
c     wvac_beg,wvac_min,wvac_end,depth, 
     .5026.158,5026.168,5026.178,0.026, 
     .5028.544,5028.559,5028.575,0.031, 
     .5028.828,5028.843,5028.863,0.029, 
     .5029.460,5029.480,5029.501,0.044, 
     .5029.460,5029.480,5029.501,0.044, 
     .5033.086,5033.086,5033.086,0.020, 
     .5034.095,5034.105,5034.115,0.025, 
     .5034.323,5034.338,5034.353,0.031, 
     .5034.495,5034.536,5034.566,0.049, 
     .5034.653,5034.673,5034.693,0.053, 
     .5035.246,5035.266,5035.286,0.043, 
     .5035.388,5035.408,5035.429,0.045, 
     .5035.677,5035.723,5035.753,0.042, 
     .5036.281,5036.296,5036.316,0.036, 
     .5037.286,5037.321,5037.352,0.111, 
     .5037.697,5037.738,5037.778,0.125, 
     .5038.261,5038.291,5038.327,0.044, 
     .5038.845,5038.860,5038.880,0.039, 
     .5039.536,5039.551,5039.566,0.029, 
     .5039.713,5039.734,5039.754,0.038, 
     .5041.381,5041.396,5041.416,0.034, 
     .5041.879,5041.904,5041.930,0.046, 
     .5042.113,5042.143,5042.164,0.057, 
     .5042.886,5042.917,5042.952,0.066, 
     .5043.054,5043.085,5043.120,0.077, 
     .5043.939,5043.970,5044.001,0.079, 
     .5052.357,5052.377,5052.403,0.044, 
     .5053.807,5053.823,5053.838,0.033, 
     .5056.404,5056.435,5056.466,0.076, 
     .5057.566,5057.596,5057.632,0.106, 
     .5058.139,5058.139,5058.139,0.020, 
     .5058.195,5058.246,5058.267,0.040, 
     .5059.982,5060.013,5060.044,0.066, 
     .5061.099,5061.115,5061.135,0.037, 
     .5065.895,5065.946,5065.977,0.081, 
     .5066.357,5066.378,5066.398,0.039, 
     .5067.061,5067.086,5067.143,0.050, 
     .5067.174,5067.194,5067.225,0.076, 
     .5068.263,5068.278,5068.289,0.026, 
     .5068.684,5068.700,5068.710,0.026, 
     .5069.090,5069.106,5069.121,0.031, 
     .5070.150,5070.160,5070.170,0.026, 
     .5069.090,5069.106,5069.121,0.031, 
     .5070.150,5070.160,5070.170,0.026, 
     .5071.220,5071.230,5071.235,0.024, 
     .5072.213,5072.223,5072.233,0.026, 
     .5072.846,5072.882,5072.918,0.106, 
     .5075.823,5075.844,5075.865,0.035, 
     .5076.236,5076.256,5076.282,0.036, 
     .5076.421,5076.463,5076.499,0.088, 
     .5077.422,5077.442,5077.463,0.046, 
     .5079.104,5079.140,5079.196,0.087, 
     .5079.238,5079.238,5079.238,0.021, 
     .5079.542,5079.558,5079.578,0.033, 
     .5083.794,5083.820,5083.841,0.051, 
     .5084.172,5084.182,5084.192,0.026, 
     .5084.508,5084.513,5084.523,0.024, 
     .5086.060,5086.070,5086.086,0.026, 
     .5086.350,5086.350,5086.350,0.020, 
     .5086.603,5086.629,5086.655,0.056, 
     .5092.251,5092.262,5092.267,0.024, 
     .5092.506,5092.526,5092.547,0.039, 
     .5094.032,5094.042,5094.057,0.027, 
     .5094.400,5094.410,5094.421,0.027, 
     .5095.870,5095.885,5095.906,0.032, 
     .5102.450,5102.461,5102.476,0.028, 
     .5408.010,5408.022,5408.034,0.026, 
     .5408.238,5408.238,5408.244,0.021, 
     .5409.473,5409.491,5409.514,0.038, 
     .5409.743,5409.760,5409.772,0.034, 
     .5412.179,5412.185,5412.197,0.022, 
     .5413.082,5413.100,5413.123,0.035, 
     .5414.337,5414.366,5414.395,0.058, 
     .5419.367,5419.391,5419.409,0.038, 
     .5420.319,5420.331,5420.343,0.023, 
     .5420.631,5420.649,5420.666,0.031, 
     .5419.367,5419.391,5419.409,0.038, 
     .5420.319,5420.331,5420.343,0.023, 
     .5420.631,5420.649,5420.666,0.031, 
     .5421.184,5421.189,5421.195,0.022, 
     .5422.925,5422.948,5422.977,0.048, 
     .5429.657,5429.663,5429.669,0.023, 
     .5434.823,5434.835,5434.847,0.026, 
     .5435.568,5435.597,5435.621,0.049, 
     .5437.081,5437.099,5437.117,0.033, 
     .5438.903,5438.903,5438.903,0.020, 
     .5442.273,5442.297,5442.314,0.032, 
     .5448.080,5448.098,5448.116,0.035, 
     .5456.358,5456.388,5456.418,0.055, 
     .5457.472,5457.478,5457.484,0.024, 
     .5459.404,5459.427,5459.445,0.037, 
     .5459.404,5459.427,5459.445,0.037, 
     .5462.028,5462.040,5462.046,0.022, 
     .5462.398,5462.428,5462.452,0.049, 
     .5465.122,5465.128,5465.140,0.024, 
     .5470.217,5470.235,5470.277,0.033, 
     .5476.756,5476.762,5476.768,0.022, 
     .5478.149,5478.149,5478.155,0.022, 
     .5676.761,5676.780,5676.793,0.028, 
     .5681.730,5681.755,5681.788,0.043, 
     .5683.758,5683.874,5683.913,0.059, 
     .5685.407,5685.433,5685.458,0.049, 
     .5687.425,5687.483,5687.548,0.125, 
     .5687.593,5687.619,5687.632,0.031, 
     .5688.584,5688.596,5688.609,0.026, 
     .5689.548,5689.594,5689.633,0.060, 
     .5690.202,5690.222,5690.248,0.035, 
     .5690.390,5690.429,5690.468,0.094, 
     .5689.548,5689.594,5689.633,0.060, 
     .5690.202,5690.222,5690.248,0.035, 
     .5690.390,5690.429,5690.468,0.094, 
     .5692.367,5692.419,5692.496,0.115, 
     .5692.730,5692.749,5692.775,0.035, 
     .5695.934,5695.960,5695.986,0.041, 
     .5698.149,5698.187,5698.226,0.094, 
     .5698.668,5698.688,5698.714,0.038, 
     .5699.266,5699.318,5699.357,0.087, 
     .5700.676,5700.722,5700.767,0.102, 
     .5714.029,5714.088,5714.101,0.025, 
     .5717.473,5717.506,5717.545,0.069, 
     .5718.919,5718.945,5718.965,0.030, 
     .5719.541,5719.580,5719.626,0.097, 
     .5720.490,5720.509,5720.536,0.032, 
     .5721.917,5721.950,5721.983,0.053, 
     .5723.870,5723.896,5723.922,0.044, 
     .5726.696,5726.703,5726.716,0.023, 
     .5726.841,5726.873,5726.906,0.069, 
     .5728.475,5728.494,5728.514,0.033, 
     .5729.637,5729.670,5729.703,0.059, 
     .5729.821,5729.841,5729.860,0.027, 
     .5729.637,5729.670,5729.703,0.059, 
     .5729.821,5729.841,5729.860,0.027, 
     .5733.100,5733.107,5733.113,0.021, 
     .5735.547,5735.574,5735.593,0.038, 
     .5737.615,5737.680,5737.727,0.089, 
     .5742.169,5742.202,5742.228,0.053, 
     .5744.948,5744.974,5744.994,0.036, 
     .5745.694,5745.786,5745.826,0.073, 
     .5752.928,5752.948,5752.974,0.039, 
     .5754.180,5754.213,5754.239,0.054, 
     .5761.582,5761.602,5761.622,0.032, 
     .5763.375,5763.395,5763.415,0.037, 
     .5767.132,5767.152,5767.172,0.030, 
     .5788.101,5788.108,5788.108,0.023, 
     .5788.175,5788.182,5788.189,0.024, 
     .5788.537,5788.550,5788.564,0.028, 
     .5788.785,5788.799,5788.805,0.023, 
     .5788.866,5788.872,5788.886,0.024, 
     .5789.174,5789.187,5789.241,0.027, 
     .5789.851,5789.858,5789.871,0.024, 
     .5790.093,5790.099,5790.113,0.026, 
     .5790.763,5790.763,5790.770,0.021, 
     .5789.174,5789.187,5789.241,0.027, 
     .5789.851,5789.858,5789.871,0.024, 
     .5790.093,5790.099,5790.113,0.026, 
     .5790.763,5790.763,5790.770,0.021, 
     .5797.527,5797.527,5797.527,0.020, 
     .5799.955,5799.962,5799.975,0.024, 
     .5800.628,5800.642,5800.648,0.025, 
     .5802.655,5802.662,5802.669,0.024, 
     .5803.315,5803.329,5803.342,0.028, 
     .5805.627,5805.634,5805.641,0.024, 
     .5806.288,5806.295,5806.302,0.024, 
     .5854.840,5854.840,5854.847,0.021, 
     .5859.582,5859.610,5859.679,0.038, 
     .5859.582,5859.610,5859.679,0.038, 
     .5861.609,5861.637,5861.658,0.035, 
     .5865.620,5865.640,5865.661,0.032, 
     .5867.575,5867.603,5867.623,0.039, 
     .5868.753,5868.760,5868.774,0.026, 
     .5870.600,5870.635,5870.669,0.068, 
     .5871.159,5871.159,5871.166,0.021, 
     .5871.269,5871.297,5871.324,0.045, 
     .5873.560,5873.580,5873.601,0.044, 
     .5875.120,5875.140,5875.161,0.034, 
     .5875.451,5875.472,5875.500,0.031, 
     .5875.555,5875.589,5875.624,0.066, 
     .5876.080,5876.128,5876.176,0.101, 
     .5876.342,5876.439,5876.480,0.077, 
     .5877.310,5877.330,5877.351,0.031, 
     .5877.551,5877.558,5877.565,0.023, 
     .5879.169,5879.190,5879.217,0.035, 
     .5879.556,5879.605,5879.771,0.103, 
     .5880.497,5880.511,5880.525,0.028, 
     .5880.705,5880.739,5880.767,0.062, 
     .5880.864,5880.940,5880.981,0.099, 
     .5879.169,5879.190,5879.217,0.035, 
     .5879.556,5879.605,5879.771,0.103, 
     .5880.497,5880.511,5880.525,0.028, 
     .5880.705,5880.739,5880.767,0.062, 
     .5880.864,5880.940,5880.981,0.099, 
     .5881.057,5881.099,5881.134,0.070, 
     .5881.715,5881.729,5881.743,0.028, 
     .5881.826,5881.867,5882.012,0.116, 
     .5882.365,5882.379,5882.393,0.026, 
     .5882.469,5882.497,5882.525,0.044, 
     .5882.746,5882.815,5882.864,0.078, 
     .5882.968,5883.009,5883.099,0.066, 
     .5883.764,5883.910,5883.993,0.253, 
     .5884.173,5884.207,5884.242,0.067, 
     .5885.482,5885.642,5885.676,0.057, 
     .5885.884,5885.981,5886.085,0.371, 
     .5886.286,5886.342,5886.453,0.136, 
     .5886.654,5886.688,5886.723,0.075, 
     .5887.132,5887.222,5887.326,0.301, 
     .5887.527,5887.666,5887.756,0.242, 
     .5887.818,5887.853,5887.888,0.048, 
     .5888.651,5888.706,5888.769,0.159, 
     .5889.067,5889.081,5889.102,0.029, 
     .5889.560,5889.643,5889.733,0.287, 
     .5889.782,5889.865,5889.900,0.059, 
     .5889.941,5889.955,5889.962,0.045, 
     .5890.004,5890.087,5890.351,0.263, 
     .5890.726,5890.740,5890.760,0.028, 
     .5891.121,5891.163,5891.198,0.068, 
     .5891.448,5891.656,5891.767,0.354, 
     .5892.253,5892.399,5892.490,0.260, 
     .5892.635,5892.635,5892.635,0.020, 
     .5892.962,5893.045,5893.108,0.153, 
     .5893.455,5893.511,5893.580,0.141, 
     .5894.352,5894.393,5894.435,0.097, 
     .5894.908,5894.950,5895.040,0.095, 
     .5895.096,5895.144,5895.179,0.077, 
     .5896.382,5896.500,5896.577,0.144, 
     .5896.626,5896.626,5896.632,0.020, 
     .5896.751,5896.834,5896.918,0.195, 
     .5897.064,5897.085,5897.106,0.034, 
     .5897.419,5897.467,5897.558,0.105, 
     .5897.732,5897.760,5897.787,0.042, 
     .5897.906,5897.940,5897.975,0.058, 
     .5898.017,5898.163,5898.442,0.360, 
     .5898.755,5898.776,5898.797,0.028, 
     .5898.894,5899.005,5899.082,0.184, 
     .5899.792,5900.043,5900.182,0.363, 
     .5900.405,5900.426,5900.440,0.026, 
     .5900.893,5900.921,5900.955,0.063, 
     .5898.998,5899.005,5899.082,0.184, 
     .5899.792,5900.043,5900.182,0.363, 
     .5900.405,5900.426,5900.440,0.026, 
     .5900.893,5900.921,5900.955,0.063, 
     .5901.213,5901.464,5901.596,0.455, 
     .5902.000,5902.140,5902.203,0.128, 
     .5902.251,5902.272,5902.286,0.030, 
     .5903.478,5903.541,5903.597,0.122, 
     .5903.673,5903.708,5903.736,0.042, 
     .5903.820,5903.848,5903.904,0.049, 
     .5903.952,5903.973,5903.987,0.030, 
     .5905.089,5905.124,5905.159,0.054, 
     .5905.257,5905.292,5905.410,0.090, 
     .5905.536,5905.550,5905.564,0.024, 
     .5906.296,5906.303,5906.310,0.022, 
     .5906.980,5907.001,5907.029,0.036, 
     .5907.232,5907.267,5907.399,0.101, 
     .5907.448,5907.476,5907.504,0.051, 
     .5907.797,5907.846,5907.888,0.129, 
     .5908.153,5908.209,5908.293,0.118, 
     .5908.915,5908.998,5909.089,0.248, 
     .5909.424,5909.452,5909.480,0.037, 
     .5909.990,5909.997,5910.011,0.023, 
     .5910.137,5910.193,5910.361,0.126, 
     .5910.459,5910.494,5910.522,0.043, 
     .5910.591,5910.766,5910.934,0.197, 
     .5911.857,5911.878,5911.899,0.037, 
     .5911.990,5912.018,5912.038,0.047, 
     .5912.542,5912.696,5912.738,0.070, 
     .5912.913,5912.997,5913.088,0.208, 
     .5914.137,5914.200,5914.256,0.158, 
     .5914.865,5914.935,5914.991,0.115, 
     .5915.369,5915.439,5915.516,0.152, 
     .5915.558,5915.635,5915.684,0.135, 
     .5915.999,5916.027,5916.055,0.041, 
     .5916.552,5916.587,5916.622,0.070, 
     .5917.351,5917.393,5917.428,0.077, 
     .5917.932,5917.953,5918.037,0.033, 
     .5918.303,5918.408,5918.612,0.250, 
     .5918.899,5919.053,5919.214,0.374, 
     .5919.502,5919.649,5919.796,0.470, 
     .5920.147,5920.168,5920.189,0.034, 
     .5920.490,5920.561,5920.638,0.163, 
     .5918.997,5919.053,5919.214,0.374, 
     .5919.502,5919.649,5919.796,0.470, 
     .5920.147,5920.168,5920.189,0.034, 
     .5920.490,5920.561,5920.638,0.163, 
     .5921.129,5921.157,5921.178,0.040, 
     .5921.620,5921.655,5921.683,0.053, 
     .5922.335,5922.518,5922.602,0.240, 
     .5922.672,5922.714,5922.756,0.076, 
     .5923.149,5923.184,5923.220,0.046, 
     .5923.563,5923.830,5923.929,0.238, 
     .5924.146,5924.273,5924.385,0.319, 
     .5924.729,5924.750,5924.771,0.029, 
     .5924.883,5924.996,5925.094,0.191, 
     .5926.590,5926.619,5926.640,0.030, 
     .5928.200,5928.291,5928.376,0.218, 
     .5928.791,5928.812,5928.840,0.033, 
     .5929.093,5929.121,5929.156,0.057, 
     .5929.403,5929.410,5929.417,0.022, 
     .5930.570,5930.599,5930.627,0.048, 
     .5930.979,5931.007,5931.035,0.042, 
     .5931.985,5932.098,5932.210,0.342, 
     .5932.696,5932.788,5932.872,0.200, 
     .5933.936,5933.936,5933.936,0.020, 
     .5934.062,5934.091,5934.119,0.048, 
     .5935.169,5935.197,5935.218,0.043, 
     .5935.768,5935.824,5935.881,0.121, 
     .5937.432,5937.432,5937.439,0.021, 
     .5937.911,5938.059,5938.102,0.098, 
     .5938.582,5938.596,5938.617,0.027, 
     .5939.936,5939.972,5940.014,0.074, 
     .5940.346,5940.424,5940.508,0.166, 
     .5940.833,5940.882,5940.995,0.080, 
     .5939.936,5939.972,5940.014,0.074, 
     .5940.346,5940.424,5940.508,0.166, 
     .5940.833,5941.087,5941.299,0.414, 
     .5941.433,5941.447,5941.461,0.023, 
     .5941.511,5941.631,5941.730,0.329, 
     .5942.245,5942.570,5942.747,0.414, 
     .5944.224,5944.316,5944.394,0.192, 
     .5944.655,5944.740,5944.818,0.171, 
     .5945.214,5945.320,5945.362,0.081, 
     .5945.575,5945.652,5945.758,0.151, 
     .5945.850,5946.013,5946.133,0.367, 
     .5946.614,5946.636,5946.657,0.030, 
     .5946.770,5947.067,5947.159,0.342, 
     .5947.449,5947.449,5947.449,0.020, 
     .5948.185,5948.228,5948.277,0.049, 
     .5948.730,5948.766,5948.801,0.055, 
     .5948.985,5949.176,5949.283,0.322, 
     .5949.424,5949.431,5949.445,0.023, 
     .5949.573,5949.828,5949.892,0.136, 
     .5950.125,5950.154,5950.175,0.044, 
     .5950.302,5950.345,5950.387,0.112, 
     .5950.862,5950.862,5950.869,0.021, 
     .5950.933,5950.968,5950.997,0.050, 
     .5951.216,5951.230,5951.237,0.023, 
     .5951.301,5951.308,5951.315,0.022, 
     .5951.443,5951.507,5951.563,0.126, 
     .5952.825,5952.832,5952.839,0.022, 
     .5953.456,5953.477,5953.499,0.040, 
     .5954.875,5954.953,5955.031,0.192, 
     .5956.308,5956.358,5956.407,0.109, 
     .5956.741,5956.783,5956.840,0.083, 
     .5957.792,5957.884,5957.991,0.258, 
     .5958.154,5958.253,5958.324,0.158, 
     .5958.544,5958.630,5958.722,0.226, 
     .5958.814,5958.836,5958.857,0.031, 
     .5959.148,5959.162,5959.177,0.025, 
     .5959.290,5959.312,5959.326,0.028, 
     .5959.958,5959.994,5960.036,0.069, 
     .5959.148,5959.162,5959.177,0.025, 
     .5959.290,5959.312,5959.326,0.028, 
     .5959.958,5959.994,5960.036,0.069, 
     .5961.408,5961.444,5961.472,0.053, 
     .5962.432,5962.468,5962.511,0.124, 
     .5964.915,5964.944,5964.972,0.047, 
     .5966.183,5966.332,5966.368,0.075, 
     .5966.582,5966.667,5966.745,0.185, 
     .5967.237,5967.237,5967.237,0.020, 
     .5967.280,5967.315,5967.344,0.056, 
     .5967.629,5967.835,5967.928,0.214, 
     .5968.028,5968.056,5968.092,0.039, 
     .5968.177,5968.277,5968.441,0.242, 
     .5968.990,5969.040,5969.090,0.104, 
     .5970.009,5970.059,5970.116,0.135, 
     .5971.264,5971.350,5971.436,0.194, 
     .5974.248,5974.276,5974.298,0.049, 
     .5975.026,5975.105,5975.240,0.192, 
     .5975.912,5975.926,5975.947,0.030, 
     .5976.490,5976.526,5976.562,0.064, 
     .5976.848,5977.012,5977.119,0.170, 
     .5977.370,5977.455,5977.484,0.041, 
     .5977.763,5977.806,5977.977,0.127, 
     .5980.123,5980.158,5980.187,0.055, 
     .5980.502,5980.516,5980.531,0.024, 
     .5980.123,5980.158,5980.187,0.055, 
     .5980.502,5980.516,5980.531,0.024, 
     .5981.196,5981.232,5981.260,0.059, 
     .5981.676,5981.683,5981.690,0.022, 
     .5981.890,5981.890,5981.890,0.020, 
     .5984.275,5984.282,5984.282,0.021, 
     .5984.949,5984.963,5984.977,0.027, 
     .5985.171,5985.228,5985.278,0.125, 
     .5986.984,5987.028,5987.071,0.060, 
     .5988.075,5988.096,5988.118,0.034, 
     .5988.519,5988.555,5988.627,0.084, 
     .5989.244,5989.294,5989.337,0.102, 
     .5990.543,5990.622,5990.651,0.048, 
     .5990.802,5990.852,5990.902,0.118, 
     .5991.958,5992.001,5992.051,0.118, 
     .5993.028,5993.050,5993.079,0.036, 
     .5993.912,5993.920,5993.927,0.023, 
     .5994.495,5994.523,5994.610,0.071, 
     .5997.314,5997.350,5997.386,0.084, 
     .5999.668,5999.733,5999.769,0.085, 
     .6000.143,6000.158,6000.172,0.025, 
     .5999.668,5999.733,5999.769,0.085, 
     .6000.143,6000.158,6000.172,0.025, 
     .6002.132,6002.139,6002.146,0.021, 
     .6002.622,6002.643,6002.665,0.036, 
     .6003.862,6003.876,6003.884,0.024, 
     .6004.641,6004.662,6004.684,0.040, 
     .6009.333,6009.362,6009.384,0.035, 
     .6015.589,6015.603,6015.618,0.027, 
     .6015.806,6015.842,6015.878,0.039, 
     .6016.001,6016.001,6016.001,0.020, 
     .6019.116,6019.124,6019.138,0.023/ 
      end 
 
 
c---------------------------------------------------------------------
c
c will flag the rej array for pixels that
c are impacted by telluric activity
c
 
      subroutine tell_rej(wav_imn,rej_imn,Ns)
      implicit none
 
      real*8 wav_imn(_IDIM_,_MDIM_,Ns)
      integer*1 rej_imn(_IDIM_,_MDIM_,Ns)
      integer Ns
 
      integer i,m,n
      real tell_depth
 
      do n = 1, Ns
         do m = 01, _MDIM_
            do i = 0001, _IDIM_
               if (tell_depth(wav_imn(i,m,n)).gt.0.001)
     .             rej_imn(i,m,n) = 9
               enddo
            enddo
         enddo
 
      return
      end
 
 
 
c---------------------------------------------------
c
c For each observed spectrum and each echelle order, compute
c the distribution of residuals (hist) and the corresponding
c cumulative distribution (cist). Print residuals at 11
c fiducial percentiles:
c
c   0.01% 0.1% 1% 10% 25% 50% 75% 90% 99% 99.9% 99.99%
c
      subroutine info_nim(pix_imn,mod_imn,wav_imn,
     .                    nrm_imn,bar_imn,rej_imn,unim,
     .                    NIT,NIMs,LABEL59,MMIN,MMAX)
      implicit none
 
      real    pix_imn(_IDIM_,_MDIM_,NIMs)
      real    mod_imn(_IDIM_,_MDIM_,NIMs)
      real*8  wav_imn(_IDIM_,_MDIM_,NIMs)
      integer*1  rej_imn(_IDIM_,_MDIM_,NIMs)
      real    nrm_imn(_IDIM_,_MDIM_,NIMs)
      real    bar_imn(_IDIM_,_MDIM_,NIMs)
      logical unim(NIMs)
      character*59 LABEL59
      integer MMIN, MMAX
 
      real*8  ww
      integer NIT, NIMs
 
      integer m, i, n
      character*80 FILENAME
      integer mu, m1, m2, mo
 
      integer hist(20001), h, hv
      integer cist(20001), ctot
      integer v, nrej
      real resv(11), res
      real vlist(11), vu
      common/vlist_/vlist
      data vlist/0.0001,0.0010,0.0100,0.1000,0.2500,0.5000,
     .           0.7500,0.9000,0.9900,0.9990,0.9999/
 
      real*8  wav0(_IDIM_,_MDIM_)
      real    dwav, ddwav
      real    di, dis
      real    dw, dws
      integer nns, imid
 
      call init_wavim(wav0,MMIN,MMAX)
 
      write(FILENAME,'(''.'',i2.2,''.nim'')') NIT
      FILENAME = TRIM(LABEL59) // FILENAME
 
      write(06,'(''INFO: '',80a)') FILENAME
 
      open(44,file=FILENAME,status='unknown')
      do n = 1, NIMs
         if (.not.unim(n)) goto 555
         dis  = 0.0
	 dws = 0.0
	 nns = 0.0
         do mu = MMIN, MMAX+1
            m1 = mu
            m2 = mu
            mo = mu
            if (mu.eq.MMAX+1) then
               m1 = MMIN
               m2 = MMAX
               mo = 0
               endif
            do h = 1, 20001
               hist(h) = 0
               enddo
            nrej = 0
            do m = m1, m2
               do i = 0001, _IDIM_
                  if (rej_imn(i,m,n).ne.0) nrej = nrej+1
                  if (rej_imn(i,m,n).eq.0) then
                     res = (pix_imn(i,m,n)-nrm_imn(i,m,n)*
     .                                     mod_imn(i,m,n))/
     .                      bar_imn(i,m,n)
                     h = 10001 + res*10000
                     if (h.lt.00002) h = 00002
                     if (h.gt.20000) h = 20000
                     hist(h) = hist(h) + 1
                     endif
                  enddo
               imid = (_WAVY_-_WAVB_)/2
               di = 2*(wav_imn(imid  ,m,n)-wav0(imid  ,m))/
     .                (   wav0(imid+1,m)  -wav0(imid-1,m))
               dw =  ((wav_imn(_WAVY_,m,n)-wav_imn(_WAVB_,m,n))-
     .                (   wav0(_WAVY_,m)  -wav0(_WAVB_,m)))/
     .                (   wav0(_WAVY_,m)  -wav0(_WAVB_,m))*
     .                (_WAVY_-_WAVB_)/2
               dis  = dis + di
               dws  = dws + dw
               nns  = nns + 1
               enddo!m1...m2
            ctot = 0
            do h = 1, 20001
               ctot = ctot + hist(h)
               cist(h) = ctot
               enddo
            do v = 1, 11
               vu = ctot*vlist(v)
               do h = 1, 20000
                  if (cist(h  ).lt.vu.and.
     .                cist(h+1).ge.vu) then
                      hv = 10001
                      if (h.lt.10001) hv = h+1
                      if (h.gt.10001) hv = h-1
                      endif
                  enddo
               resv(v) = (hv-10001)*0.0001
               enddo
            write(44,144) n,mo,(resv(v),v=1,11),nrej,ctot,
     .                    dis/nns,dws/nns
            if (mo.eq.0)
     .      write(06,144) n,mo,(resv(v),v=1,11),nrej,ctot,
     .                    dis/nns,dws/nns
            enddo!m
  144       format(i3.3,1x,i2.2,1x,
     .             5(1x,f7.4),4x,f7.4,4x,5(f7.4,1x),i5,1x,i8,1x,
     .             f8.3,1x,f8.3)
  555    continue!not unim
         enddo!n
      close(44)
 
      return
      end
 
 
 
 
      subroutine fill_mod_imn(mod_imn,wav_imn,lsfqkmn,tems,z0,unim,
     .                        Ns,MMIN,MMAX)
      implicit none
 
      real    mod_imn(_IDIM_,_MDIM_,Ns)
      real*8  wav_imn(_IDIM_,_MDIM_,Ns)
      real*8  z0(Ns)
      real    lsfqkmn(_QDIM_,_KDIM_,_MDIM_,Ns)
      real    tems(_SDIM_)
      logical unim(Ns)
      integer Ns
      integer MMIN, MMAX
 
      integer n, m
      real    fxto(_ODIM_)
      real    lsfji(_JDIM_,_IDIM_)
 
      do n = 1, Ns
         if (unim(n)) then
            if (n.eq.n/10*10) write(6,*) '   n: ',n
            do m = MMIN, MMAX
               call lsfqk2lsfji(lsfqkmn(1,1,m,n),lsfji)          ! make a fine-scale LSF
               call find_fxto(fxto,wav_imn(1,m,n),z0(n),tems)    ! fxto is the super sam x10
               call find_modi(mod_imn(1,m,n),fxto,lsfji)         ! find the model of the obsns
               enddo
            endif
         enddo
 
      return
      end
 
 
 
c
c this will determine which pixels can contribute to
c the template
c
c
      subroutine find_uuuimn(pix_imn,mod_imn,rat_imn,
     .                       bar_imn,rej_imn,uuu_imn,met_imn,
     .                       wav_imn,z0,unim,Ns,MMIN,MMAX)
      implicit none
 
      real    pix_imn(_IDIM_,_MDIM_,Ns)
      real    mod_imn(_IDIM_,_MDIM_,Ns)
      real    rat_imn(_IDIM_,_MDIM_,Ns)
      real    bar_imn(_IDIM_,_MDIM_,Ns)
      integer*1  rej_imn(_IDIM_,_MDIM_,Ns)
      integer*1  uuu_imn(_IDIM_,_MDIM_,Ns)
      integer*1  met_imn(_IDIM_,_MDIM_,Ns)
      real*8  wav_imn(_IDIM_,_MDIM_,Ns)
      real*8  z0(Ns)
      logical unim(Ns)
      integer Ns
      integer MMIN, MMAX
 
      integer s, n, m, i, i_w
      real*8  w, ww_s
 
      integer U, Us, Uu, Ux
      real*8 ru(9999), rbar, rsig
      integer iu(9999)
      integer mu(9999)
      integer nu(9999)
 
      do n = 01, Ns
      do m = 01, _MDIM_
      do i = 01, _IDIM_
         uuu_imn(i,m,n) = rej_imn(i,m,n)
         if (.not.unim(n)) uuu_imn(i,m,n) = 9
         enddo
         enddo
         enddo
 
 
c
c could make it so that only one pixel in each
c exposure contributes ; not that way now
c
      write(6,*) '---> find_uuuimn...'
      open(91,file='LOG.UNIM',status='unknown')
      do s = 00001, _SDIM_
         w = ww_s(s)
         Us = 0
         do n = 01, Ns
            if (unim(n)) then
               do m = MMIN, MMAX
                  if (wav_imn(_RATB_,m,n)/(1+z0(n)).lt.w.and.
     .                wav_imn(_RATY_,m,n)/(1+z0(n)).gt.w) then
                      i = i_w(w/(1+z0(n)),wav_imn(1,m,n))
                      if (rej_imn(i,m,n).eq.0) then
                         Us = Us + 1
                         ru(Us) = (pix_imn(i,m,n)-
     .                             mod_imn(i,m,n)*rat_imn(i,m,n))/
     .                             bar_imn(i,m,n)
                         iu(Us) = i
                         mu(Us) = m
                         nu(Us) = n
                         if (.not.(ru(Us).gt.-9e9)) then
                            Us = Us-1
                            uuu_imn(i,m,n) = 5
                            endif
                         endif
                      endif
                  enddo!m
               endif
            enddo!n
         call barsig(ru,Us,rbar,rsig,Uu)
         if (s.eq.s/1000*1000) write(06,111) s,w,rbar,rsig,Us,Uu
         Ux = 0
         do u = 1, Us
            if (abs(ru(u)-rbar).gt.2.25*rsig.and.
     .          iu(u).gt._RATB_.and.
     .          iu(u).lt._RATY_) then
                uuu_imn(iu(u),mu(u),nu(u)) = 1
                Ux = Ux + 1
                endif
            if (abs(ru(u)-rbar).gt.7.50*rsig.and.           ! higher threshold for
     .          (iu(u).le._RATB_.or.iu(u).ge._RATY_)) then  ! rejecting in vicinity
                uuu_imn(iu(u),mu(u),nu(u)) = 1              ! of the ends
                Ux = Ux + 1                              
                endif
            enddo
         write(91,111) s,w,rbar,rsig,Us,Uu,Ux
  111    format(i6.6,1x,f10.5,1x,f8.4,1x,f8.4,1x,3i5)
         enddo!s
      close(91)

      do n = 01, Ns
      do m = 01, _MDIM_
      do i = 01, _IDIM_
         if (met_imn(i,m,n).ge.2) uuu_imn(i,m,n) = 1
         enddo
         enddo
         enddo
 
      return
      end
 
 
 
c--------------------------------------------------------------------
c
c
c
      subroutine fill_rat_imn(pix_imn,mod_imn,
     .                        bar_imn,rat_imn,
     .                        rej_imn,met_imn,unim,Ns,
     .                        NIT,MMIN,MMAX)
      implicit none
 
      real    pix_imn(_IDIM_,_MDIM_,Ns)
      real    mod_imn(_IDIM_,_MDIM_,Ns)
      real    bar_imn(_IDIM_,_MDIM_,Ns)
      real    rat_imn(_IDIM_,_MDIM_,Ns)
      integer*1  rej_imn(_IDIM_,_MDIM_,Ns)
      integer*1  met_imn(_IDIM_,_MDIM_,Ns)
      logical unim(Ns)
      integer Ns, NIT
      integer MMIN, MMAX
 
      integer n, m
 
      do n = 1, Ns
         if (unim(n)) then
            do m = MMIN, MMAX
c              call find_rati(pix_imn(1,m,n),mod_imn(1,m,n),     ! find the scale that takes
c    .                        rat_imn(1,m,n),rej_imn(1,m,n),NIT) !    model*rat ---> pix
               call find_rati_bar(pix_imn(1,m,n),mod_imn(1,m,n),   ! find the scale that takes
     .                            bar_imn(1,m,n),rat_imn(1,m,n),   !    model*rat ---> pix
     .                            rej_imn(1,m,n),met_imn(1,m,n),NIT)                  
               call find_bari(pix_imn(1,m,n),bar_imn(1,m,n),     ! average pixel value
     .                                       rej_imn(1,m,n))     !    to help scaling
               enddo
            endif
         enddo
 
      return
      end
 
 
      logical function all_ok(array,Ns)
      implicit none
 
      integer Ns
      real array(Ns)
 
      integer n
 
      all_ok = .true.
      do n = 1, Ns
         if (.not.(array(n).ge.-9e9)) all_ok = .false.
         enddo
 
      return
      end
 
 
 
c--------------------------------------------------------------------
c
c block data stuff
c
c
c
c To use this, take w(i) and set:
c
c   A = sum[w(i)*abcde_i(1,i)/1.0d9] i=0021,4001
c   B = sum[w(i)*abcde_i(2,i)/1.0d9] i=0021,4001
c   C = sum[w(i)*abcde_i(3,i)/1.0d9] i=0021,4001
c   D = sum[w(i)*abcde_i(4,i)/1.0d9] i=0021,4001
c   E = sum[w(i)*abcde_i(5,i)/1.0d9] i=0021,4001
c
c Then:  rx = (i-2011)/1990.0d0
c
c  wnew(i) = A + B*rx    + C*rx**2
c              + D*rx**3 + E*xr**4
c
c
      block data wav_abcde_bd
      real*8 abcde_i(5,4021)
      common /wav_abcde_comm/abcde_i
      data    abcde_i /   ! * 1e9
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0001
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0002
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0003
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0004
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0005
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0006
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0007
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0008
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0009
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0010
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0011
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0012
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0013
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0014
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0015
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0016
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0017
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0018
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0019
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   0020
     . 46933356d0, 187969431d0,-657067409d0,-438595524d0, 985601736d0,   !   0021  1.00000
     . 46603006d0, 187213442d0,-653103699d0,-437273125d0, 980648963d0,   !   0022  1.00000
     . 46273569d0, 186458451d0,-649149453d0,-435952388d0, 975707393d0,   !   0023  1.00000
     . 45945045d0, 185704457d0,-645204661d0,-434633312d0, 970777013d0,   !   0024  1.00000
     . 45617432d0, 184951460d0,-641269310d0,-433315896d0, 965857808d0,   !   0025  1.00000
     . 45290729d0, 184199459d0,-637343390d0,-432000140d0, 960949767d0,   !   0026  1.00000
     . 44964936d0, 183448453d0,-633426890d0,-430686041d0, 956052876d0,   !   0027  1.00000
     . 44640050d0, 182698442d0,-629519797d0,-429373601d0, 951167122d0,   !   0028  1.00000
     . 44316072d0, 181949426d0,-625622102d0,-428062817d0, 946292491d0,   !   0029  1.00000
     . 43992999d0, 181201403d0,-621733792d0,-426753688d0, 941428972d0,   !   0030  1.00000
     . 43670830d0, 180454375d0,-617854856d0,-425446215d0, 936576550d0,   !   0031  1.00000
     . 43349566d0, 179708339d0,-613985284d0,-424140396d0, 931735212d0,   !   0032  1.00000
     . 43029203d0, 178963296d0,-610125064d0,-422836231d0, 926904946d0,   !   0033  1.00000
     . 42709742d0, 178219245d0,-606274184d0,-421533718d0, 922085738d0,   !   0034  1.00000
     . 42391181d0, 177476185d0,-602432635d0,-420232856d0, 917277575d0,   !   0035  1.00000
     . 42073519d0, 176734116d0,-598600403d0,-418933646d0, 912480445d0,   !   0036  1.00000
     . 41756755d0, 175993038d0,-594777479d0,-417636085d0, 907694334d0,   !   0037  1.00000
     . 41440888d0, 175252950d0,-590963850d0,-416340174d0, 902919229d0,   !   0038  1.00000
     . 41125916d0, 174513851d0,-587159507d0,-415045911d0, 898155117d0,   !   0039  1.00000
     . 40811839d0, 173775742d0,-583364437d0,-413753296d0, 893401985d0,   !   0040  1.00000
     . 40498656d0, 173038621d0,-579578629d0,-412462327d0, 888659820d0,   !   0041  1.00000
     . 40186364d0, 172302488d0,-575802073d0,-411173004d0, 883928609d0,   !   0042  1.00000
     . 39874964d0, 171567342d0,-572034757d0,-409885327d0, 879208339d0,   !   0043  1.00000
     . 39564454d0, 170833184d0,-568276669d0,-408599293d0, 874498997d0,   !   0044  1.00000
     . 39254833d0, 170100011d0,-564527800d0,-407314903d0, 869800570d0,   !   0045  1.00000
     . 38946100d0, 169367825d0,-560788137d0,-406032155d0, 865113045d0,   !   0046  1.00000
     . 38638253d0, 168636625d0,-557057670d0,-404751049d0, 860436409d0,   !   0047  1.00000
     . 38331292d0, 167906409d0,-553336387d0,-403471583d0, 855770649d0,   !   0048  1.00000
     . 38025216d0, 167177178d0,-549624278d0,-402193758d0, 851115751d0,   !   0049  1.00000
     . 37720023d0, 166448931d0,-545921330d0,-400917572d0, 846471704d0,   !   0050  1.00000
     . 37415712d0, 165721667d0,-542227534d0,-399643024d0, 841838494d0,   !   0051  1.00000
     . 37112282d0, 164995387d0,-538542877d0,-398370114d0, 837216108d0,   !   0052  1.00000
     . 36809733d0, 164270089d0,-534867349d0,-397098840d0, 832604534d0,   !   0053  1.00000
     . 36508062d0, 163545772d0,-531200939d0,-395829202d0, 828003757d0,   !   0054  1.00000
     . 36207269d0, 162822438d0,-527543636d0,-394561199d0, 823413766d0,   !   0055  1.00000
     . 35907353d0, 162100084d0,-523895428d0,-393294831d0, 818834547d0,   !   0056  1.00000
     . 35608313d0, 161378711d0,-520256304d0,-392030095d0, 814266088d0,   !   0057  1.00000
     . 35310147d0, 160658318d0,-516626254d0,-390766992d0, 809708375d0,   !   0058  1.00000
     . 35012854d0, 159938904d0,-513005266d0,-389505521d0, 805161396d0,   !   0059  1.00000
     . 34716434d0, 159220469d0,-509393329d0,-388245681d0, 800625137d0,   !   0060  1.00000
     . 34420885d0, 158503013d0,-505790433d0,-386987470d0, 796099586d0,   !   0061  1.00000
     . 34126206d0, 157786535d0,-502196565d0,-385730889d0, 791584730d0,   !   0062  1.00000
     . 33832396d0, 157071034d0,-498611716d0,-384475936d0, 787080556d0,   !   0063  1.00000
     . 33539454d0, 156356510d0,-495035873d0,-383222610d0, 782587052d0,   !   0064  1.00000
     . 33247379d0, 155642962d0,-491469027d0,-381970911d0, 778104203d0,   !   0065  1.00000
     . 32956169d0, 154930391d0,-487911165d0,-380720838d0, 773631998d0,   !   0066  1.00000
     . 32665824d0, 154218795d0,-484362278d0,-379472390d0, 769170424d0,   !   0067  1.00000
     . 32376342d0, 153508174d0,-480822353d0,-378225566d0, 764719467d0,   !   0068  1.00000
     . 32087723d0, 152798528d0,-477291380d0,-376980365d0, 760279115d0,   !   0069  1.00000
     . 31799965d0, 152089855d0,-473769348d0,-375736787d0, 755849355d0,   !   0070  1.00000
     . 31513067d0, 151382156d0,-470256246d0,-374494830d0, 751430175d0,   !   0071  1.00000
     . 31227028d0, 150675430d0,-466752063d0,-373254494d0, 747021560d0,   !   0072  1.00000
     . 30941847d0, 149969676d0,-463256788d0,-372015778d0, 742623500d0,   !   0073  1.00000
     . 30657523d0, 149264894d0,-459770409d0,-370778682d0, 738235980d0,   !   0074  1.00000
     . 30374055d0, 148561084d0,-456292916d0,-369543203d0, 733858988d0,   !   0075  1.00000
     . 30091441d0, 147858245d0,-452824299d0,-368309342d0, 729492511d0,   !   0076  1.00000
     . 29809681d0, 147156376d0,-449364545d0,-367077097d0, 725136536d0,   !   0077  1.00000
     . 29528773d0, 146455477d0,-445913644d0,-365846469d0, 720791051d0,   !   0078  1.00000
     . 29248717d0, 145755548d0,-442471586d0,-364617455d0, 716456043d0,   !   0079  1.00000
     . 28969511d0, 145056587d0,-439038358d0,-363390055d0, 712131499d0,   !   0080  1.00000
     . 28691154d0, 144358595d0,-435613951d0,-362164269d0, 707817406d0,   !   0081  1.00000
     . 28413645d0, 143661571d0,-432198352d0,-360940095d0, 703513751d0,   !   0082  1.00000
     . 28136983d0, 142965514d0,-428791552d0,-359717532d0, 699220522d0,   !   0083  1.00000
     . 27861167d0, 142270424d0,-425393539d0,-358496581d0, 694937707d0,   !   0084  1.00000
     . 27586196d0, 141576301d0,-422004303d0,-357277239d0, 690665291d0,   !   0085  1.00000
     . 27312069d0, 140883144d0,-418623832d0,-356059507d0, 686403263d0,   !   0086  1.00000
     . 27038784d0, 140190952d0,-415252116d0,-354843382d0, 682151610d0,   !   0087  1.00000
     . 26766341d0, 139499725d0,-411889143d0,-353628865d0, 677910319d0,   !   0088  1.00000
     . 26494738d0, 138809462d0,-408534904d0,-352415955d0, 673679377d0,   !   0089  1.00000
     . 26223974d0, 138120164d0,-405189386d0,-351204650d0, 669458772d0,   !   0090  1.00000
     . 25954049d0, 137431828d0,-401852579d0,-349994950d0, 665248491d0,   !   0091  1.00000
     . 25684961d0, 136744456d0,-398524472d0,-348786855d0, 661048521d0,   !   0092  1.00000
     . 25416709d0, 136058046d0,-395205054d0,-347580362d0, 656858849d0,   !   0093  1.00000
     . 25149292d0, 135372598d0,-391894315d0,-346375472d0, 652679464d0,   !   0094  1.00000
     . 24882709d0, 134688112d0,-388592243d0,-345172184d0, 648510352d0,   !   0095  1.00000
     . 24616958d0, 134004586d0,-385298827d0,-343970496d0, 644351500d0,   !   0096  1.00000
     . 24352040d0, 133322021d0,-382014058d0,-342770409d0, 640202896d0,   !   0097  1.00000
     . 24087952d0, 132640416d0,-378737923d0,-341571920d0, 636064528d0,   !   0098  1.00000
     . 23824693d0, 131959770d0,-375470413d0,-340375029d0, 631936382d0,   !   0099  1.00000
     . 23562263d0, 131280083d0,-372211515d0,-339179736d0, 627818446d0,   !   0100  1.00000
     . 23300661d0, 130601355d0,-368961220d0,-337986040d0, 623710708d0,   !   0101  1.00000
     . 23039885d0, 129923584d0,-365719517d0,-336793939d0, 619613154d0,   !   0102  1.00000
     . 22779934d0, 129246771d0,-362486394d0,-335603433d0, 615525772d0,   !   0103  1.00000
     . 22520807d0, 128570915d0,-359261841d0,-334414521d0, 611448549d0,   !   0104  1.00000
     . 22262503d0, 127896015d0,-356045848d0,-333227202d0, 607381474d0,   !   0105  1.00000
     . 22005021d0, 127222071d0,-352838402d0,-332041476d0, 603324533d0,   !   0106  1.00000
     . 21748360d0, 126549082d0,-349639494d0,-330857341d0, 599277713d0,   !   0107  1.00000
     . 21492519d0, 125877049d0,-346449113d0,-329674797d0, 595241003d0,   !   0108  1.00000
     . 21237496d0, 125205970d0,-343267248d0,-328493843d0, 591214389d0,   !   0109  1.00000
     . 20983291d0, 124535844d0,-340093888d0,-327314478d0, 587197859d0,   !   0110  1.00000
     . 20729903d0, 123866672d0,-336929022d0,-326136701d0, 583191400d0,   !   0111  1.00000
     . 20477331d0, 123198453d0,-333772640d0,-324960511d0, 579195000d0,   !   0112  1.00000
     . 20225572d0, 122531187d0,-330624730d0,-323785908d0, 575208646d0,   !   0113  1.00000
     . 19974627d0, 121864872d0,-327485283d0,-322612891d0, 571232327d0,   !   0114  1.00000
     . 19724495d0, 121199509d0,-324354287d0,-321441459d0, 567266028d0,   !   0115  1.00000
     . 19475173d0, 120535096d0,-321231731d0,-320271610d0, 563309738d0,   !   0116  1.00000
     . 19226662d0, 119871634d0,-318117606d0,-319103345d0, 559363444d0,   !   0117  1.00000
     . 18978960d0, 119209122d0,-315011899d0,-317936662d0, 555427133d0,   !   0118  1.00000
     . 18732065d0, 118547559d0,-311914601d0,-316771561d0, 551500794d0,   !   0119  1.00000
     . 18485978d0, 117886945d0,-308825700d0,-315608040d0, 547584413d0,   !   0120  1.00000
     . 18240696d0, 117227280d0,-305745187d0,-314446099d0, 543677978d0,   !   0121  1.00000
     . 17996219d0, 116568562d0,-302673049d0,-313285737d0, 539781477d0,   !   0122  1.00000
     . 17752546d0, 115910792d0,-299609277d0,-312126954d0, 535894897d0,   !   0123  1.00000
     . 17509676d0, 115253969d0,-296553860d0,-310969747d0, 532018225d0,   !   0124  1.00000
     . 17267607d0, 114598092d0,-293506786d0,-309814117d0, 528151450d0,   !   0125  1.00000
     . 17026338d0, 113943161d0,-290468047d0,-308660063d0, 524294558d0,   !   0126  1.00000
     . 16785869d0, 113289175d0,-287437629d0,-307507584d0, 520447537d0,   !   0127  1.00000
     . 16546198d0, 112636134d0,-284415524d0,-306356679d0, 516610375d0,   !   0128  1.00000
     . 16307325d0, 111984038d0,-281401720d0,-305207346d0, 512783060d0,   !   0129  1.00000
     . 16069247d0, 111332885d0,-278396207d0,-304059587d0, 508965578d0,   !   0130  1.00000
     . 15831965d0, 110682676d0,-275398974d0,-302913398d0, 505157918d0,   !   0131  1.00000
     . 15595477d0, 110033409d0,-272410010d0,-301768781d0, 501360066d0,   !   0132  1.00000
     . 15359782d0, 109385085d0,-269429305d0,-300625733d0, 497572011d0,   !   0133  1.00000
     . 15124879d0, 108737703d0,-266456848d0,-299484254d0, 493793741d0,   !   0134  1.00000
     . 14890766d0, 108091263d0,-263492628d0,-298344343d0, 490025242d0,   !   0135  1.00000
     . 14657444d0, 107445763d0,-260536636d0,-297206000d0, 486266502d0,   !   0136  1.00000
     . 14424910d0, 106801203d0,-257588859d0,-296069223d0, 482517510d0,   !   0137  1.00000
     . 14193164d0, 106157583d0,-254649288d0,-294934012d0, 478778252d0,   !   0138  1.00000
     . 13962204d0, 105514903d0,-251717911d0,-293800365d0, 475048716d0,   !   0139  1.00000
     . 13732030d0, 104873162d0,-248794719d0,-292668283d0, 471328890d0,   !   0140  1.00000
     . 13502641d0, 104232359d0,-245879701d0,-291537763d0, 467618761d0,   !   0141  1.00000
     . 13274035d0, 103592493d0,-242972846d0,-290408806d0, 463918317d0,   !   0142  1.00000
     . 13046211d0, 102953565d0,-240074143d0,-289281411d0, 460227546d0,   !   0143  1.00000
     . 12819169d0, 102315574d0,-237183582d0,-288155576d0, 456546436d0,   !   0144  1.00000
     . 12592907d0, 101678520d0,-234301152d0,-287031300d0, 452874973d0,   !   0145  1.00000
     . 12367424d0, 101042401d0,-231426843d0,-285908584d0, 449213146d0,   !   0146  1.00000
     . 12142719d0, 100407217d0,-228560645d0,-284787426d0, 445560943d0,   !   0147  1.00000
     . 11918792d0,  99772968d0,-225702545d0,-283667825d0, 441918350d0,   !   0148  1.00000
     . 11695640d0,  99139654d0,-222852535d0,-282549781d0, 438285357d0,   !   0149  1.00000
     . 11473263d0,  98507274d0,-220010604d0,-281433292d0, 434661949d0,   !   0150  1.00000
     . 11251660d0,  97875826d0,-217176740d0,-280318358d0, 431048116d0,   !   0151  1.00000
     . 11030831d0,  97245312d0,-214350934d0,-279204978d0, 427443845d0,   !   0152  1.00000
     . 10810772d0,  96615730d0,-211533174d0,-278093152d0, 423849123d0,   !   0153  1.00000
     . 10591485d0,  95987080d0,-208723451d0,-276982877d0, 420263938d0,   !   0154  1.00000
     . 10372967d0,  95359361d0,-205921753d0,-275874154d0, 416688278d0,   !   0155  1.00000
     . 10155218d0,  94732573d0,-203128071d0,-274766981d0, 413122131d0,   !   0156  1.00000
     .  9938237d0,  94106716d0,-200342394d0,-273661359d0, 409565485d0,   !   0157  1.00000
     .  9722022d0,  93481788d0,-197564711d0,-272557285d0, 406018326d0,   !   0158  1.00000
     .  9506573d0,  92857789d0,-194795011d0,-271454759d0, 402480644d0,   !   0159  1.00000
     .  9291888d0,  92234720d0,-192033285d0,-270353781d0, 398952425d0,   !   0160  1.00000
     .  9077966d0,  91612578d0,-189279521d0,-269254349d0, 395433657d0,   !   0161  1.00000
     .  8864807d0,  90991365d0,-186533710d0,-268156462d0, 391924329d0,   !   0162  1.00000
     .  8652409d0,  90371078d0,-183795840d0,-267060120d0, 388424428d0,   !   0163  1.00000
     .  8440771d0,  89751719d0,-181065902d0,-265965323d0, 384933942d0,   !   0164  1.00000
     .  8229893d0,  89133286d0,-178343885d0,-264872068d0, 381452858d0,   !   0165  1.00000
     .  8019773d0,  88515779d0,-175629778d0,-263780356d0, 377981164d0,   !   0166  1.00000
     .  7810409d0,  87899197d0,-172923570d0,-262690185d0, 374518849d0,   !   0167  1.00000
     .  7601802d0,  87283540d0,-170225253d0,-261601554d0, 371065899d0,   !   0168  1.00000
     .  7393950d0,  86668807d0,-167534814d0,-260514464d0, 367622304d0,   !   0169  1.00000
     .  7186852d0,  86054998d0,-164852243d0,-259428912d0, 364188050d0,   !   0170  1.00000
     .  6980506d0,  85442112d0,-162177531d0,-258344898d0, 360763125d0,   !   0171  1.00000
     .  6774913d0,  84830149d0,-159510667d0,-257262422d0, 357347518d0,   !   0172  1.00000
     .  6570070d0,  84219109d0,-156851639d0,-256181482d0, 353941216d0,   !   0173  1.00000
     .  6365978d0,  83608990d0,-154200439d0,-255102077d0, 350544207d0,   !   0174  1.00000
     .  6162634d0,  82999792d0,-151557055d0,-254024207d0, 347156478d0,   !   0175  1.00000
     .  5960037d0,  82391516d0,-148921477d0,-252947871d0, 343778018d0,   !   0176  1.00000
     .  5758188d0,  81784159d0,-146293695d0,-251873069d0, 340408815d0,   !   0177  1.00000
     .  5557084d0,  81177723d0,-143673697d0,-250799798d0, 337048857d0,   !   0178  1.00000
     .  5356724d0,  80572205d0,-141061475d0,-249728059d0, 333698130d0,   !   0179  1.00000
     .  5157108d0,  79967607d0,-138457017d0,-248657850d0, 330356624d0,   !   0180  1.00000
     .  4958235d0,  79363927d0,-135860313d0,-247589171d0, 327024326d0,   !   0181  1.00000
     .  4760103d0,  78761164d0,-133271353d0,-246522021d0, 323701224d0,   !   0182  1.00000
     .  4562711d0,  78159319d0,-130690126d0,-245456400d0, 320387306d0,   !   0183  1.00000
     .  4366059d0,  77558391d0,-128116622d0,-244392305d0, 317082559d0,   !   0184  1.00000
     .  4170145d0,  76958379d0,-125550830d0,-243329736d0, 313786973d0,   !   0185  1.00000
     .  3974969d0,  76359283d0,-122992741d0,-242268694d0, 310500534d0,   !   0186  1.00000
     .  3780529d0,  75761102d0,-120442344d0,-241209175d0, 307223231d0,   !   0187  1.00000
     .  3586824d0,  75163835d0,-117899628d0,-240151181d0, 303955051d0,   !   0188  1.00000
     .  3393853d0,  74567484d0,-115364583d0,-239094710d0, 300695983d0,   !   0189  1.00000
     .  3201616d0,  73972045d0,-112837199d0,-238039761d0, 297446014d0,   !   0190  1.00000
     .  3010111d0,  73377521d0,-110317466d0,-236986333d0, 294205132d0,   !   0191  1.00000
     .  2819336d0,  72783908d0,-107805373d0,-235934426d0, 290973326d0,   !   0192  1.00000
     .  2629292d0,  72191209d0,-105300910d0,-234884039d0, 287750584d0,   !   0193  1.00000
     .  2439977d0,  71599421d0,-102804067d0,-233835170d0, 284536892d0,   !   0194  1.00000
     .  2251390d0,  71008544d0,-100314833d0,-232787820d0, 281332240d0,   !   0195  1.00000
     .  2063530d0,  70418578d0, -97833198d0,-231741986d0, 278136616d0,   !   0196  1.00000
     .  1876396d0,  69829522d0, -95359151d0,-230697669d0, 274950007d0,   !   0197  1.00000
     .  1689987d0,  69241376d0, -92892683d0,-229654868d0, 271772401d0,   !   0198  1.00000
     .  1504302d0,  68654140d0, -90433784d0,-228613581d0, 268603786d0,   !   0199  1.00000
     .  1319340d0,  68067812d0, -87982442d0,-227573808d0, 265444151d0,   !   0200  1.00000
     .  1135099d0,  67482392d0, -85538648d0,-226535548d0, 262293483d0,   !   0201  1.00000
     .   951580d0,  66897880d0, -83102391d0,-225498800d0, 259151770d0,   !   0202  1.00000
     .   768780d0,  66314275d0, -80673661d0,-224463564d0, 256019002d0,   !   0203  1.00000
     .   586699d0,  65731578d0, -78252448d0,-223429838d0, 252895164d0,   !   0204  1.00000
     .   405336d0,  65149786d0, -75838742d0,-222397622d0, 249780246d0,   !   0205  1.00000
     .   224689d0,  64568900d0, -73432532d0,-221366915d0, 246674236d0,   !   0206  1.00000
     .    44758d0,  63988919d0, -71033808d0,-220337715d0, 243577122d0,   !   0207  1.00000
     .  -134457d0,  63409844d0, -68642560d0,-219310023d0, 240488891d0,   !   0208  1.00000
     .  -312960d0,  62831672d0, -66258778d0,-218283837d0, 237409532d0,   !   0209  1.00000
     .  -490749d0,  62254404d0, -63882452d0,-217259156d0, 234339033d0,   !   0210  1.00000
     .  -667828d0,  61678040d0, -61513570d0,-216235981d0, 231277382d0,   !   0211  1.00000
     .  -844195d0,  61102578d0, -59152124d0,-215214309d0, 228224568d0,   !   0212  1.00000
     . -1019853d0,  60528018d0, -56798103d0,-214194140d0, 225180577d0,   !   0213  1.00000
     . -1194803d0,  59954360d0, -54451496d0,-213175473d0, 222145399d0,   !   0214  1.00000
     . -1369045d0,  59381604d0, -52112293d0,-212158307d0, 219119021d0,   !   0215  1.00000
     . -1542580d0,  58809748d0, -49780486d0,-211142642d0, 216101432d0,   !   0216  1.00000
     . -1715410d0,  58238792d0, -47456062d0,-210128477d0, 213092620d0,   !   0217  1.00000
     . -1887536d0,  57668736d0, -45139012d0,-209115811d0, 210092573d0,   !   0218  1.00000
     . -2058958d0,  57099579d0, -42829326d0,-208104642d0, 207101278d0,   !   0219  1.00000
     . -2229677d0,  56531321d0, -40526993d0,-207094970d0, 204118725d0,   !   0220  1.00000
     . -2399695d0,  55963962d0, -38232004d0,-206086795d0, 201144902d0,   !   0221  1.00000
     . -2569013d0,  55397500d0, -35944348d0,-205080116d0, 198179795d0,   !   0222  1.00000
     . -2737632d0,  54831935d0, -33664016d0,-204074931d0, 195223395d0,   !   0223  1.00000
     . -2905551d0,  54267267d0, -31390997d0,-203071239d0, 192275688d0,   !   0224  1.00000
     . -3072774d0,  53703495d0, -29125280d0,-202069041d0, 189336663d0,   !   0225  1.00000
     . -3239300d0,  53140619d0, -26866856d0,-201068335d0, 186406309d0,   !   0226  1.00000
     . -3405131d0,  52578638d0, -24615715d0,-200069120d0, 183484613d0,   !   0227  1.00000
     . -3570268d0,  52017552d0, -22371847d0,-199071395d0, 180571563d0,   !   0228  1.00000
     . -3734711d0,  51457360d0, -20135241d0,-198075161d0, 177667149d0,   !   0229  1.00000
     . -3898462d0,  50898062d0, -17905887d0,-197080415d0, 174771357d0,   !   0230  1.00000
     . -4061521d0,  50339657d0, -15683776d0,-196087156d0, 171884177d0,   !   0231  1.00000
     . -4223890d0,  49782145d0, -13468896d0,-195095386d0, 169005597d0,   !   0232  1.00000
     . -4385570d0,  49225525d0, -11261239d0,-194105101d0, 166135604d0,   !   0233  1.00000
     . -4546562d0,  48669796d0,  -9060794d0,-193116302d0, 163274187d0,   !   0234  1.00000
     . -4706867d0,  48114959d0,  -6867551d0,-192128987d0, 160421335d0,   !   0235  1.00000
     . -4866485d0,  47561013d0,  -4681499d0,-191143156d0, 157577035d0,   !   0236  1.00000
     . -5025418d0,  47007957d0,  -2502630d0,-190158808d0, 154741276d0,   !   0237  1.00000
     . -5183666d0,  46455790d0,   -330932d0,-189175942d0, 151914046d0,   !   0238  1.00000
     . -5341232d0,  45904513d0,   1833603d0,-188194558d0, 149095333d0,   !   0239  1.00000
     . -5498115d0,  45354124d0,   3990988d0,-187214653d0, 146285126d0,   !   0240  1.00000
     . -5654317d0,  44804624d0,   6141231d0,-186236229d0, 143483413d0,   !   0241  1.00000
     . -5809839d0,  44256011d0,   8284342d0,-185259283d0, 140690182d0,   !   0242  1.00000
     . -5964681d0,  43708286d0,  10420332d0,-184283815d0, 137905422d0,   !   0243  1.00000
     . -6118845d0,  43161447d0,  12549211d0,-183309824d0, 135129120d0,   !   0244  1.00000
     . -6272333d0,  42615494d0,  14670988d0,-182337310d0, 132361266d0,   !   0245  1.00000
     . -6425143d0,  42070427d0,  16785673d0,-181366271d0, 129601847d0,   !   0246  1.00000
     . -6577279d0,  41526246d0,  18893277d0,-180396706d0, 126850852d0,   !   0247  1.00000
     . -6728740d0,  40982948d0,  20993809d0,-179428615d0, 124108269d0,   !   0248  1.00000
     . -6879529d0,  40440536d0,  23087280d0,-178461997d0, 121374087d0,   !   0249  1.00000
     . -7029645d0,  39899006d0,  25173700d0,-177496852d0, 118648293d0,   !   0250  1.00000
     . -7179089d0,  39358361d0,  27253078d0,-176533177d0, 115930877d0,   !   0251  1.00000
     . -7327864d0,  38818597d0,  29325424d0,-175570973d0, 113221826d0,   !   0252  1.00000
     . -7475969d0,  38279716d0,  31390749d0,-174610238d0, 110521130d0,   !   0253  1.00000
     . -7623406d0,  37741717d0,  33449063d0,-173650972d0, 107828775d0,   !   0254  1.00000
     . -7770176d0,  37204599d0,  35500374d0,-172693174d0, 105144752d0,   !   0255  1.00000
     . -7916280d0,  36668361d0,  37544694d0,-171736843d0, 102469047d0,   !   0256  1.00000
     . -8061718d0,  36133004d0,  39582033d0,-170781979d0,  99801650d0,   !   0257  1.00000
     . -8206492d0,  35598526d0,  41612399d0,-169828579d0,  97142549d0,   !   0258  1.00000
     . -8350603d0,  35064927d0,  43635804d0,-168876645d0,  94491732d0,   !   0259  1.00000
     . -8494052d0,  34532208d0,  45652257d0,-167926174d0,  91849188d0,   !   0260  1.00000
     . -8636839d0,  34000366d0,  47661768d0,-166977165d0,  89214905d0,   !   0261  1.00000
     . -8778966d0,  33469402d0,  49664346d0,-166029619d0,  86588872d0,   !   0262  1.00000
     . -8920434d0,  32939315d0,  51660003d0,-165083534d0,  83971077d0,   !   0263  1.00000
     . -9061244d0,  32410105d0,  53648747d0,-164138910d0,  81361508d0,   !   0264  1.00000
     . -9201396d0,  31881771d0,  55630589d0,-163195745d0,  78760154d0,   !   0265  1.00000
     . -9340892d0,  31354312d0,  57605539d0,-162254039d0,  76167003d0,   !   0266  1.00000
     . -9479732d0,  30827729d0,  59573606d0,-161313790d0,  73582045d0,   !   0267  1.00000
     . -9617919d0,  30302021d0,  61534800d0,-160374999d0,  71005267d0,   !   0268  1.00000
     . -9755452d0,  29777186d0,  63489132d0,-159437664d0,  68436657d0,   !   0269  1.00000
     . -9892332d0,  29253226d0,  65436611d0,-158501784d0,  65876205d0,   !   0270  1.00000
     .-10028561d0,  28730138d0,  67377247d0,-157567359d0,  63323898d0,   !   0271  1.00000
     .-10164140d0,  28207923d0,  69311049d0,-156634387d0,  60779726d0,   !   0272  1.00000
     .-10299070d0,  27686580d0,  71238028d0,-155702868d0,  58243677d0,   !   0273  1.00000
     .-10433351d0,  27166109d0,  73158194d0,-154772802d0,  55715739d0,   !   0274  1.00000
     .-10566985d0,  26646509d0,  75071556d0,-153844186d0,  53195900d0,   !   0275  1.00000
     .-10699972d0,  26127780d0,  76978125d0,-152917021d0,  50684150d0,   !   0276  1.00000
     .-10832314d0,  25609921d0,  78877910d0,-151991305d0,  48180477d0,   !   0277  1.00000
     .-10964011d0,  25092931d0,  80770920d0,-151067038d0,  45684870d0,   !   0278  1.00000
     .-11095065d0,  24576811d0,  82657166d0,-150144219d0,  43197316d0,   !   0279  1.00000
     .-11225476d0,  24061559d0,  84536658d0,-149222847d0,  40717805d0,   !   0280  1.00000
     .-11355246d0,  23547176d0,  86409405d0,-148302921d0,  38246324d0,   !   0281  1.00000
     .-11484375d0,  23033660d0,  88275417d0,-147384441d0,  35782864d0,   !   0282  1.00000
     .-11612865d0,  22521011d0,  90134705d0,-146467405d0,  33327411d0,   !   0283  1.00000
     .-11740716d0,  22009229d0,  91987277d0,-145551812d0,  30879955d0,   !   0284  1.00000
     .-11867930d0,  21498313d0,  93833143d0,-144637663d0,  28440485d0,   !   0285  1.00000
     .-11994507d0,  20988262d0,  95672314d0,-143724955d0,  26008989d0,   !   0286  1.00000
     .-12120448d0,  20479077d0,  97504800d0,-142813689d0,  23585455d0,   !   0287  1.00000
     .-12245754d0,  19970757d0,  99330609d0,-141903863d0,  21169872d0,   !   0288  1.00000
     .-12370427d0,  19463300d0, 101149751d0,-140995476d0,  18762229d0,   !   0289  1.00000
     .-12494468d0,  18956707d0, 102962237d0,-140088528d0,  16362514d0,   !   0290  1.00000
     .-12617876d0,  18450977d0, 104768077d0,-139183018d0,  13970716d0,   !   0291  1.00000
     .-12740654d0,  17946110d0, 106567279d0,-138278945d0,  11586824d0,   !   0292  1.00000
     .-12862802d0,  17442105d0, 108359854d0,-137376308d0,   9210826d0,   !   0293  1.00000
     .-12984321d0,  16938962d0, 110145811d0,-136475106d0,   6842711d0,   !   0294  1.00000
     .-13105212d0,  16436680d0, 111925161d0,-135575339d0,   4482467d0,   !   0295  1.00000
     .-13225477d0,  15935258d0, 113697912d0,-134677005d0,   2130084d0,   !   0296  1.00000
     .-13345116d0,  15434696d0, 115464075d0,-133780104d0,   -214449d0,   !   0297  1.00000
     .-13464129d0,  14934995d0, 117223660d0,-132884636d0,  -2551146d0,   !   0298  1.00000
     .-13582519d0,  14436152d0, 118976675d0,-131990598d0,  -4880017d0,   !   0299  1.00000
     .-13700286d0,  13938168d0, 120723131d0,-131097991d0,  -7201073d0,   !   0300  1.00000
     .-13817431d0,  13441041d0, 122463038d0,-130206813d0,  -9514326d0,   !   0301  1.00000
     .-13933954d0,  12944773d0, 124196405d0,-129317064d0, -11819786d0,   !   0302  1.00000
     .-14049858d0,  12449362d0, 125923241d0,-128428742d0, -14117465d0,   !   0303  1.00000
     .-14165143d0,  11954807d0, 127643557d0,-127541848d0, -16407374d0,   !   0304  1.00000
     .-14279809d0,  11461108d0, 129357363d0,-126656379d0, -18689525d0,   !   0305  1.00000
     .-14393859d0,  10968265d0, 131064667d0,-125772336d0, -20963929d0,   !   0306  1.00000
     .-14507292d0,  10476277d0, 132765480d0,-124889718d0, -23230598d0,   !   0307  1.00000
     .-14620110d0,   9985144d0, 134459810d0,-124008522d0, -25489541d0,   !   0308  1.00000
     .-14732313d0,   9494865d0, 136147669d0,-123128750d0, -27740771d0,   !   0309  1.00000
     .-14843904d0,   9005439d0, 137829065d0,-122250400d0, -29984300d0,   !   0310  1.00000
     .-14954882d0,   8516866d0, 139504009d0,-121373471d0, -32220137d0,   !   0311  1.00000
     .-15065248d0,   8029146d0, 141172509d0,-120497962d0, -34448295d0,   !   0312  1.00000
     .-15175005d0,   7542278d0, 142834576d0,-119623872d0, -36668784d0,   !   0313  1.00000
     .-15284152d0,   7056262d0, 144490218d0,-118751202d0, -38881617d0,   !   0314  1.00000
     .-15392690d0,   6571097d0, 146139446d0,-117879948d0, -41086803d0,   !   0315  1.00000
     .-15500621d0,   6086782d0, 147782270d0,-117010112d0, -43284355d0,   !   0316  1.00000
     .-15607945d0,   5603318d0, 149418698d0,-116141692d0, -45474283d0,   !   0317  1.00000
     .-15714664d0,   5120703d0, 151048741d0,-115274687d0, -47656599d0,   !   0318  1.00000
     .-15820778d0,   4638937d0, 152672408d0,-114409097d0, -49831314d0,   !   0319  1.00000
     .-15926288d0,   4158020d0, 154289709d0,-113544920d0, -51998439d0,   !   0320  1.00000
     .-16031196d0,   3677951d0, 155900653d0,-112682156d0, -54157986d0,   !   0321  1.00000
     .-16135502d0,   3198729d0, 157505250d0,-111820804d0, -56309964d0,   !   0322  1.00000
     .-16239207d0,   2720354d0, 159103509d0,-110960863d0, -58454387d0,   !   0323  1.00000
     .-16342313d0,   2242826d0, 160695440d0,-110102332d0, -60591264d0,   !   0324  1.00000
     .-16444819d0,   1766144d0, 162281053d0,-109245211d0, -62720607d0,   !   0325  1.00000
     .-16546728d0,   1290308d0, 163860357d0,-108389498d0, -64842427d0,   !   0326  1.00000
     .-16648039d0,    815317d0, 165433361d0,-107535193d0, -66956736d0,   !   0327  1.00000
     .-16748755d0,    341170d0, 167000076d0,-106682295d0, -69063543d0,   !   0328  1.00000
     .-16848875d0,   -132131d0, 168560510d0,-105830803d0, -71162861d0,   !   0329  1.00000
     .-16948402d0,   -604590d0, 170114674d0,-104980716d0, -73254701d0,   !   0330  1.00000
     .-17047335d0,  -1076206d0, 171662577d0,-104132034d0, -75339074d0,   !   0331  1.00000
     .-17145676d0,  -1546980d0, 173204228d0,-103284755d0, -77415990d0,   !   0332  1.00000
     .-17243426d0,  -2016911d0, 174739637d0,-102438879d0, -79485461d0,   !   0333  1.00000
     .-17340586d0,  -2486001d0, 176268813d0,-101594405d0, -81547498d0,   !   0334  1.00000
     .-17437156d0,  -2954250d0, 177791766d0,-100751332d0, -83602113d0,   !   0335  1.00000
     .-17533138d0,  -3421658d0, 179308505d0, -99909659d0, -85649315d0,   !   0336  1.00000
     .-17628532d0,  -3888226d0, 180819041d0, -99069385d0, -87689117d0,   !   0337  1.00000
     .-17723340d0,  -4353955d0, 182323381d0, -98230510d0, -89721529d0,   !   0338  1.00000
     .-17817562d0,  -4818844d0, 183821537d0, -97393033d0, -91746562d0,   !   0339  1.00000
     .-17911200d0,  -5282895d0, 185313517d0, -96556952d0, -93764228d0,   !   0340  1.00000
     .-18004254d0,  -5746108d0, 186799330d0, -95722268d0, -95774537d0,   !   0341  1.00000
     .-18096726d0,  -6208483d0, 188278987d0, -94888979d0, -97777501d0,   !   0342  1.00000
     .-18188615d0,  -6670022d0, 189752497d0, -94057084d0, -99773130d0,   !   0343  1.00000
     .-18279924d0,  -7130724d0, 191219869d0, -93226583d0,-101761436d0,   !   0344  1.00000
     .-18370653d0,  -7590589d0, 192681113d0, -92397475d0,-103742429d0,   !   0345  1.00000
     .-18460803d0,  -8049620d0, 194136237d0, -91569758d0,-105716121d0,   !   0346  1.00000
     .-18550375d0,  -8507815d0, 195585253d0, -90743432d0,-107682523d0,   !   0347  1.00000
     .-18639370d0,  -8965175d0, 197028168d0, -89918497d0,-109641645d0,   !   0348  1.00000
     .-18727790d0,  -9421702d0, 198464992d0, -89094951d0,-111593499d0,   !   0349  1.00000
     .-18815633d0,  -9877395d0, 199895736d0, -88272793d0,-113538096d0,   !   0350  1.00000
     .-18902903d0, -10332255d0, 201320407d0, -87452024d0,-115475446d0,   !   0351  1.00000
     .-18989600d0, -10786282d0, 202739017d0, -86632641d0,-117405560d0,   !   0352  1.00000
     .-19075724d0, -11239478d0, 204151573d0, -85814644d0,-119328451d0,   !   0353  1.00000
     .-19161276d0, -11691842d0, 205558086d0, -84998032d0,-121244127d0,   !   0354  1.00000
     .-19246259d0, -12143375d0, 206958564d0, -84182805d0,-123152601d0,   !   0355  1.00000
     .-19330672d0, -12594077d0, 208353018d0, -83368961d0,-125053884d0,   !   0356  1.00000
     .-19414516d0, -13043949d0, 209741457d0, -82556499d0,-126947985d0,   !   0357  1.00000
     .-19497792d0, -13492992d0, 211123889d0, -81745420d0,-128834917d0,   !   0358  1.00000
     .-19580502d0, -13941206d0, 212500325d0, -80935722d0,-130714690d0,   !   0359  1.00000
     .-19662647d0, -14388591d0, 213870774d0, -80127404d0,-132587316d0,   !   0360  1.00000
     .-19744226d0, -14835149d0, 215235244d0, -79320465d0,-134452804d0,   !   0361  1.00000
     .-19825242d0, -15280878d0, 216593746d0, -78514904d0,-136311166d0,   !   0362  1.00000
     .-19905694d0, -15725781d0, 217946289d0, -77710722d0,-138162413d0,   !   0363  1.00000
     .-19985585d0, -16169857d0, 219292882d0, -76907916d0,-140006556d0,   !   0364  1.00000
     .-20064914d0, -16613108d0, 220633535d0, -76106486d0,-141843605d0,   !   0365  1.00000
     .-20143684d0, -17055532d0, 221968256d0, -75306431d0,-143673572d0,   !   0366  1.00000
     .-20221894d0, -17497132d0, 223297055d0, -74507751d0,-145496468d0,   !   0367  1.00000
     .-20299546d0, -17937907d0, 224619942d0, -73710444d0,-147312302d0,   !   0368  1.00000
     .-20376640d0, -18377858d0, 225936926d0, -72914510d0,-149121087d0,   !   0369  1.00000
     .-20453179d0, -18816986d0, 227248016d0, -72119947d0,-150922832d0,   !   0370  1.00000
     .-20529161d0, -19255291d0, 228553221d0, -71326756d0,-152717549d0,   !   0371  1.00000
     .-20604590d0, -19692773d0, 229852551d0, -70534935d0,-154505249d0,   !   0372  1.00000
     .-20679464d0, -20129433d0, 231146014d0, -69744483d0,-156285943d0,   !   0373  1.00000
     .-20753786d0, -20565272d0, 232433621d0, -68955399d0,-158059641d0,   !   0374  1.00000
     .-20827556d0, -21000289d0, 233715381d0, -68167684d0,-159826353d0,   !   0375  1.00000
     .-20900775d0, -21434486d0, 234991302d0, -67381335d0,-161586092d0,   !   0376  1.00000
     .-20973445d0, -21867864d0, 236261395d0, -66596351d0,-163338868d0,   !   0377  1.00000
     .-21045565d0, -22300421d0, 237525668d0, -65812734d0,-165084691d0,   !   0378  1.00000
     .-21117138d0, -22732160d0, 238784130d0, -65030480d0,-166823573d0,   !   0379  1.00000
     .-21188163d0, -23163080d0, 240036791d0, -64249590d0,-168555523d0,   !   0380  1.00000
     .-21258642d0, -23593182d0, 241283661d0, -63470062d0,-170280554d0,   !   0381  1.00000
     .-21328576d0, -24022467d0, 242524747d0, -62691897d0,-171998675d0,   !   0382  1.00000
     .-21397965d0, -24450935d0, 243760061d0, -61915092d0,-173709898d0,   !   0383  1.00000
     .-21466811d0, -24878586d0, 244989610d0, -61139647d0,-175414233d0,   !   0384  1.00000
     .-21535115d0, -25305421d0, 246213404d0, -60365562d0,-177111691d0,   !   0385  1.00000
     .-21602877d0, -25731441d0, 247431452d0, -59592835d0,-178802282d0,   !   0386  1.00000
     .-21670098d0, -26156646d0, 248643764d0, -58821466d0,-180486019d0,   !   0387  1.00000
     .-21736780d0, -26581036d0, 249850349d0, -58051453d0,-182162910d0,   !   0388  1.00000
     .-21802922d0, -27004613d0, 251051215d0, -57282796d0,-183832967d0,   !   0389  1.00000
     .-21868527d0, -27427376d0, 252246373d0, -56515495d0,-185496201d0,   !   0390  1.00000
     .-21933595d0, -27849326d0, 253435830d0, -55749548d0,-187152623d0,   !   0391  1.00000
     .-21998127d0, -28270463d0, 254619598d0, -54984954d0,-188802242d0,   !   0392  1.00000
     .-22062124d0, -28690788d0, 255797683d0, -54221713d0,-190445070d0,   !   0393  1.00000
     .-22125586d0, -29110302d0, 256970097d0, -53459824d0,-192081118d0,   !   0394  1.00000
     .-22188515d0, -29529005d0, 258136848d0, -52699285d0,-193710396d0,   !   0395  1.00000
     .-22250912d0, -29946898d0, 259297944d0, -51940097d0,-195332915d0,   !   0396  1.00000
     .-22312777d0, -30363981d0, 260453396d0, -51182258d0,-196948685d0,   !   0397  1.00000
     .-22374112d0, -30780254d0, 261603213d0, -50425767d0,-198557718d0,   !   0398  1.00000
     .-22434917d0, -31195718d0, 262747403d0, -49670624d0,-200160024d0,   !   0399  1.00000
     .-22495193d0, -31610374d0, 263885975d0, -48916828d0,-201755613d0,   !   0400  1.00000
     .-22554942d0, -32024221d0, 265018940d0, -48164377d0,-203344496d0,   !   0401  1.00000
     .-22614164d0, -32437262d0, 266146305d0, -47413272d0,-204926684d0,   !   0402  1.00000
     .-22672859d0, -32849495d0, 267268081d0, -46663511d0,-206502187d0,   !   0403  1.00000
     .-22731030d0, -33260922d0, 268384276d0, -45915093d0,-208071017d0,   !   0404  1.00000
     .-22788676d0, -33671543d0, 269494899d0, -45168018d0,-209633183d0,   !   0405  1.00000
     .-22845799d0, -34081359d0, 270599960d0, -44422285d0,-211188697d0,   !   0406  1.00000
     .-22902400d0, -34490369d0, 271699467d0, -43677892d0,-212737569d0,   !   0407  1.00000
     .-22958480d0, -34898575d0, 272793430d0, -42934840d0,-214279809d0,   !   0408  1.00000
     .-23014039d0, -35305977d0, 273881857d0, -42193127d0,-215815428d0,   !   0409  1.00000
     .-23069078d0, -35712576d0, 274964758d0, -41452752d0,-217344437d0,   !   0410  1.00000
     .-23123598d0, -36118372d0, 276042142d0, -40713715d0,-218866847d0,   !   0411  1.00000
     .-23177601d0, -36523365d0, 277114019d0, -39976015d0,-220382667d0,   !   0412  1.00000
     .-23231087d0, -36927557d0, 278180396d0, -39239650d0,-221891909d0,   !   0413  1.00000
     .-23284057d0, -37330947d0, 279241283d0, -38504621d0,-223394582d0,   !   0414  1.00000
     .-23336512d0, -37733536d0, 280296689d0, -37770925d0,-224890699d0,   !   0415  1.00000
     .-23388453d0, -38135324d0, 281346624d0, -37038564d0,-226380268d0,   !   0416  1.00000
     .-23439880d0, -38536313d0, 282391095d0, -36307534d0,-227863302d0,   !   0417  1.00000
     .-23490795d0, -38936502d0, 283430113d0, -35577837d0,-229339809d0,   !   0418  1.00000
     .-23541198d0, -39335892d0, 284463686d0, -34849470d0,-230809801d0,   !   0419  1.00000
     .-23591091d0, -39734484d0, 285491824d0, -34122433d0,-232273289d0,   !   0420  1.00000
     .-23640475d0, -40132278d0, 286514535d0, -33396726d0,-233730282d0,   !   0421  1.00000
     .-23689349d0, -40529275d0, 287531828d0, -32672347d0,-235180792d0,   !   0422  1.00000
     .-23737716d0, -40925474d0, 288543712d0, -31949295d0,-236624828d0,   !   0423  1.00000
     .-23785576d0, -41320877d0, 289550197d0, -31227571d0,-238062402d0,   !   0424  1.00000
     .-23832930d0, -41715485d0, 290551292d0, -30507172d0,-239493524d0,   !   0425  1.00000
     .-23879778d0, -42109296d0, 291547004d0, -29788098d0,-240918204d0,   !   0426  1.00000
     .-23926123d0, -42502313d0, 292537344d0, -29070348d0,-242336453d0,   !   0427  1.00000
     .-23971964d0, -42894536d0, 293522320d0, -28353922d0,-243748280d0,   !   0428  1.00000
     .-24017302d0, -43285964d0, 294501942d0, -27638818d0,-245153698d0,   !   0429  1.00000
     .-24062139d0, -43676600d0, 295476218d0, -26925036d0,-246552716d0,   !   0430  1.00000
     .-24106475d0, -44066442d0, 296445157d0, -26212574d0,-247945344d0,   !   0431  1.00000
     .-24150312d0, -44455492d0, 297408768d0, -25501433d0,-249331594d0,   !   0432  1.00000
     .-24193650d0, -44843750d0, 298367060d0, -24791611d0,-250711475d0,   !   0433  1.00000
     .-24236489d0, -45231216d0, 299320042d0, -24083107d0,-252084998d0,   !   0434  1.00000
     .-24278832d0, -45617892d0, 300267724d0, -23375921d0,-253452173d0,   !   0435  1.00000
     .-24320678d0, -46003777d0, 301210113d0, -22670052d0,-254813011d0,   !   0436  1.00000
     .-24362029d0, -46388872d0, 302147219d0, -21965498d0,-256167523d0,   !   0437  1.00000
     .-24402886d0, -46773178d0, 303079051d0, -21262259d0,-257515718d0,   !   0438  1.00000
     .-24443249d0, -47156695d0, 304005618d0, -20560335d0,-258857607d0,   !   0439  1.00000
     .-24483120d0, -47539424d0, 304926928d0, -19859723d0,-260193200d0,   !   0440  1.00000
     .-24522499d0, -47921365d0, 305842991d0, -19160425d0,-261522508d0,   !   0441  1.00000
     .-24561387d0, -48302518d0, 306753816d0, -18462438d0,-262845541d0,   !   0442  1.00000
     .-24599785d0, -48682884d0, 307659410d0, -17765762d0,-264162310d0,   !   0443  1.00000
     .-24637694d0, -49062464d0, 308559784d0, -17070395d0,-265472825d0,   !   0444  1.00000
     .-24675115d0, -49441258d0, 309454946d0, -16376338d0,-266777096d0,   !   0445  1.00000
     .-24712049d0, -49819267d0, 310344905d0, -15683590d0,-268075133d0,   !   0446  1.00000
     .-24748496d0, -50196491d0, 311229670d0, -14992149d0,-269366948d0,   !   0447  1.00000
     .-24784457d0, -50572931d0, 312109250d0, -14302014d0,-270652550d0,   !   0448  1.00000
     .-24819935d0, -50948586d0, 312983653d0, -13613185d0,-271931949d0,   !   0449  1.00000
     .-24854928d0, -51323458d0, 313852889d0, -12925661d0,-273205157d0,   !   0450  1.00000
     .-24889438d0, -51697548d0, 314716966d0, -12239442d0,-274472182d0,   !   0451  1.00000
     .-24923467d0, -52070855d0, 315575893d0, -11554526d0,-275733037d0,   !   0452  1.00000
     .-24957014d0, -52443380d0, 316429679d0, -10870912d0,-276987730d0,   !   0453  1.00000
     .-24990082d0, -52815123d0, 317278333d0, -10188600d0,-278236272d0,   !   0454  1.00000
     .-25022670d0, -53186086d0, 318121863d0,  -9507588d0,-279478674d0,   !   0455  1.00000
     .-25054779d0, -53556269d0, 318960279d0,  -8827877d0,-280714946d0,   !   0456  1.00000
     .-25086411d0, -53925671d0, 319793589d0,  -8149465d0,-281945097d0,   !   0457  1.00000
     .-25117567d0, -54294294d0, 320621803d0,  -7472351d0,-283169140d0,   !   0458  1.00000
     .-25148246d0, -54662139d0, 321444928d0,  -6796534d0,-284387082d0,   !   0459  1.00000
     .-25178451d0, -55029205d0, 322262974d0,  -6122014d0,-285598936d0,   !   0460  1.00000
     .-25208182d0, -55395493d0, 323075949d0,  -5448790d0,-286804710d0,   !   0461  1.00000
     .-25237439d0, -55761004d0, 323883863d0,  -4776861d0,-288004416d0,   !   0462  1.00000
     .-25266224d0, -56125737d0, 324686724d0,  -4106226d0,-289198064d0,   !   0463  1.00000
     .-25294538d0, -56489695d0, 325484540d0,  -3436884d0,-290385663d0,   !   0464  1.00000
     .-25322381d0, -56852877d0, 326277321d0,  -2768835d0,-291567225d0,   !   0465  1.00000
     .-25349755d0, -57215283d0, 327065076d0,  -2102077d0,-292742759d0,   !   0466  1.00000
     .-25376660d0, -57576914d0, 327847813d0,  -1436610d0,-293912275d0,   !   0467  1.00000
     .-25403097d0, -57937771d0, 328625540d0,   -772434d0,-295075784d0,   !   0468  1.00000
     .-25429067d0, -58297855d0, 329398267d0,   -109546d0,-296233296d0,   !   0469  1.00000
     .-25454571d0, -58657165d0, 330166003d0,    552053d0,-297384821d0,   !   0470  1.00000
     .-25479610d0, -59015702d0, 330928755d0,   1212365d0,-298530369d0,   !   0471  1.00000
     .-25504184d0, -59373466d0, 331686534d0,   1871390d0,-299669951d0,   !   0472  1.00000
     .-25528295d0, -59730459d0, 332439347d0,   2529130d0,-300803576d0,   !   0473  1.00000
     .-25551943d0, -60086681d0, 333187203d0,   3185585d0,-301931255d0,   !   0474  1.00000
     .-25575130d0, -60442131d0, 333930112d0,   3840756d0,-303052998d0,   !   0475  1.00000
     .-25597855d0, -60796812d0, 334668081d0,   4494643d0,-304168816d0,   !   0476  1.00000
     .-25620121d0, -61150722d0, 335401119d0,   5147248d0,-305278717d0,   !   0477  1.00000
     .-25641927d0, -61503863d0, 336129236d0,   5798572d0,-306382713d0,   !   0478  1.00000
     .-25663275d0, -61856235d0, 336852439d0,   6448614d0,-307480814d0,   !   0479  1.00000
     .-25684166d0, -62207839d0, 337570738d0,   7097377d0,-308573029d0,   !   0480  1.00000
     .-25704600d0, -62558675d0, 338284141d0,   7744861d0,-309659369d0,   !   0481  1.00000
     .-25724578d0, -62908744d0, 338992657d0,   8391066d0,-310739844d0,   !   0482  1.00000
     .-25744102d0, -63258046d0, 339696294d0,   9035994d0,-311814464d0,   !   0483  1.00000
     .-25763171d0, -63606582d0, 340395061d0,   9679646d0,-312883239d0,   !   0484  1.00000
     .-25781788d0, -63954351d0, 341088968d0,  10322021d0,-313946179d0,   !   0485  1.00000
     .-25799952d0, -64301356d0, 341778022d0,  10963122d0,-315003295d0,   !   0486  1.00000
     .-25817665d0, -64647596d0, 342462232d0,  11602949d0,-316054596d0,   !   0487  1.00000
     .-25834928d0, -64993071d0, 343141606d0,  12241502d0,-317100092d0,   !   0488  1.00000
     .-25851741d0, -65337782d0, 343816154d0,  12878784d0,-318139794d0,   !   0489  1.00000
     .-25868105d0, -65681730d0, 344485884d0,  13514793d0,-319173712d0,   !   0490  1.00000
     .-25884021d0, -66024916d0, 345150805d0,  14149532d0,-320201855d0,   !   0491  1.00000
     .-25899491d0, -66367339d0, 345810925d0,  14783001d0,-321224234d0,   !   0492  1.00000
     .-25914514d0, -66709000d0, 346466253d0,  15415201d0,-322240859d0,   !   0493  1.00000
     .-25929092d0, -67049900d0, 347116798d0,  16046132d0,-323251739d0,   !   0494  1.00000
     .-25943225d0, -67390039d0, 347762567d0,  16675796d0,-324256886d0,   !   0495  1.00000
     .-25956915d0, -67729418d0, 348403570d0,  17304194d0,-325256308d0,   !   0496  1.00000
     .-25970162d0, -68068037d0, 349039816d0,  17931326d0,-326250016d0,   !   0497  1.00000
     .-25982967d0, -68405896d0, 349671312d0,  18557193d0,-327238020d0,   !   0498  1.00000
     .-25995331d0, -68742997d0, 350298068d0,  19181796d0,-328220329d0,   !   0499  1.00000
     .-26007255d0, -69079340d0, 350920092d0,  19805136d0,-329196955d0,   !   0500  1.00000
     .-26018740d0, -69414924d0, 351537392d0,  20427214d0,-330167906d0,   !   0501  1.00000
     .-26029786d0, -69749752d0, 352149977d0,  21048030d0,-331133194d0,   !   0502  1.00000
     .-26040395d0, -70083822d0, 352757857d0,  21667585d0,-332092827d0,   !   0503  1.00000
     .-26050567d0, -70417137d0, 353361038d0,  22285880d0,-333046816d0,   !   0504  1.00000
     .-26060303d0, -70749695d0, 353959530d0,  22902916d0,-333995170d0,   !   0505  1.00000
     .-26069604d0, -71081498d0, 354553342d0,  23518695d0,-334937900d0,   !   0506  1.00000
     .-26078471d0, -71412547d0, 355142481d0,  24133215d0,-335875016d0,   !   0507  1.00000
     .-26086904d0, -71742841d0, 355726957d0,  24746480d0,-336806528d0,   !   0508  1.00000
     .-26094905d0, -72072381d0, 356306777d0,  25358488d0,-337732445d0,   !   0509  1.00000
     .-26102475d0, -72401168d0, 356881951d0,  25969242d0,-338652777d0,   !   0510  1.00000
     .-26109613d0, -72729202d0, 357452487d0,  26578741d0,-339567535d0,   !   0511  1.00000
     .-26116322d0, -73056484d0, 358018394d0,  27186988d0,-340476728d0,   !   0512  1.00000
     .-26122602d0, -73383015d0, 358579679d0,  27793982d0,-341380367d0,   !   0513  1.00000
     .-26128454d0, -73708793d0, 359136352d0,  28399725d0,-342278460d0,   !   0514  1.00000
     .-26133878d0, -74033822d0, 359688421d0,  29004217d0,-343171018d0,   !   0515  1.00000
     .-26138876d0, -74358099d0, 360235894d0,  29607459d0,-344058052d0,   !   0516  1.00000
     .-26143448d0, -74681627d0, 360778780d0,  30209452d0,-344939570d0,   !   0517  1.00000
     .-26147595d0, -75004406d0, 361317087d0,  30810197d0,-345815583d0,   !   0518  1.00000
     .-26151318d0, -75326436d0, 361850824d0,  31409694d0,-346686100d0,   !   0519  1.00000
     .-26154619d0, -75647718d0, 362379999d0,  32007945d0,-347551132d0,   !   0520  1.00000
     .-26157496d0, -75968252d0, 362904621d0,  32604951d0,-348410688d0,   !   0521  1.00000
     .-26159953d0, -76288038d0, 363424699d0,  33200711d0,-349264778d0,   !   0522  1.00000
     .-26161989d0, -76607078d0, 363940239d0,  33795228d0,-350113412d0,   !   0523  1.00000
     .-26163605d0, -76925372d0, 364451252d0,  34388501d0,-350956599d0,   !   0524  1.00000
     .-26164802d0, -77242920d0, 364957746d0,  34980533d0,-351794351d0,   !   0525  1.00000
     .-26165581d0, -77559722d0, 365459728d0,  35571322d0,-352626676d0,   !   0526  1.00000
     .-26165943d0, -77875780d0, 365957207d0,  36160871d0,-353453584d0,   !   0527  1.00000
     .-26165889d0, -78191094d0, 366450193d0,  36749180d0,-354275085d0,   !   0528  1.00000
     .-26165419d0, -78505664d0, 366938692d0,  37336251d0,-355091189d0,   !   0529  1.00000
     .-26164535d0, -78819490d0, 367422714d0,  37922083d0,-355901906d0,   !   0530  1.00000
     .-26163236d0, -79132574d0, 367902267d0,  38506677d0,-356707245d0,   !   0531  1.00000
     .-26161525d0, -79444915d0, 368377359d0,  39090036d0,-357507216d0,   !   0532  1.00000
     .-26159401d0, -79756515d0, 368847999d0,  39672158d0,-358301829d0,   !   0533  1.00000
     .-26156866d0, -80067374d0, 369314195d0,  40253046d0,-359091095d0,   !   0534  1.00000
     .-26153921d0, -80377491d0, 369775956d0,  40832700d0,-359875021d0,   !   0535  1.00000
     .-26150566d0, -80686869d0, 370233289d0,  41411120d0,-360653619d0,   !   0536  1.00000
     .-26146802d0, -80995507d0, 370686204d0,  41988309d0,-361426898d0,   !   0537  1.00000
     .-26142630d0, -81303405d0, 371134708d0,  42564266d0,-362194868d0,   !   0538  1.00000
     .-26138051d0, -81610565d0, 371578810d0,  43138992d0,-362957539d0,   !   0539  1.00000
     .-26133066d0, -81916986d0, 372018518d0,  43712489d0,-363714920d0,   !   0540  1.00000
     .-26127675d0, -82222670d0, 372453841d0,  44284756d0,-364467020d0,   !   0541  1.00000
     .-26121879d0, -82527617d0, 372884787d0,  44855796d0,-365213851d0,   !   0542  1.00000
     .-26115680d0, -82831826d0, 373311364d0,  45425608d0,-365955421d0,   !   0543  1.00000
     .-26109078d0, -83135300d0, 373733581d0,  45994194d0,-366691740d0,   !   0544  1.00000
     .-26102073d0, -83438038d0, 374151445d0,  46561554d0,-367422818d0,   !   0545  1.00000
     .-26094668d0, -83740040d0, 374564966d0,  47127689d0,-368148664d0,   !   0546  1.00000
     .-26086862d0, -84041308d0, 374974151d0,  47692600d0,-368869289d0,   !   0547  1.00000
     .-26078656d0, -84341842d0, 375379009d0,  48256289d0,-369584701d0,   !   0548  1.00000
     .-26070051d0, -84641641d0, 375779549d0,  48818755d0,-370294911d0,   !   0549  1.00000
     .-26061048d0, -84940708d0, 376175777d0,  49379999d0,-370999929d0,   !   0550  1.00000
     .-26051649d0, -85239042d0, 376567703d0,  49940023d0,-371699763d0,   !   0551  1.00000
     .-26041853d0, -85536644d0, 376955336d0,  50498827d0,-372394424d0,   !   0552  1.00000
     .-26031661d0, -85833514d0, 377338682d0,  51056413d0,-373083921d0,   !   0553  1.00000
     .-26021075d0, -86129652d0, 377717751d0,  51612780d0,-373768263d0,   !   0554  1.00000
     .-26010094d0, -86425060d0, 378092551d0,  52167930d0,-374447462d0,   !   0555  1.00000
     .-25998721d0, -86719738d0, 378463089d0,  52721863d0,-375121525d0,   !   0556  1.00000
     .-25986956d0, -87013686d0, 378829375d0,  53274581d0,-375790463d0,   !   0557  1.00000
     .-25974799d0, -87306905d0, 379191417d0,  53826084d0,-376454285d0,   !   0558  1.00000
     .-25962252d0, -87599395d0, 379549222d0,  54376373d0,-377113002d0,   !   0559  1.00000
     .-25949315d0, -87891157d0, 379902799d0,  54925450d0,-377766621d0,   !   0560  1.00000
     .-25935989d0, -88182191d0, 380252157d0,  55473313d0,-378415154d0,   !   0561  1.00000
     .-25922275d0, -88472498d0, 380597303d0,  56019966d0,-379058610d0,   !   0562  1.00000
     .-25908173d0, -88762079d0, 380938245d0,  56565408d0,-379696998d0,   !   0563  1.00000
     .-25893685d0, -89050933d0, 381274992d0,  57109640d0,-380330328d0,   !   0564  1.00000
     .-25878812d0, -89339061d0, 381607553d0,  57652663d0,-380958609d0,   !   0565  1.00000
     .-25863554d0, -89626464d0, 381935934d0,  58194478d0,-381581851d0,   !   0566  1.00000
     .-25847912d0, -89913143d0, 382260145d0,  58735086d0,-382200064d0,   !   0567  1.00000
     .-25831886d0, -90199097d0, 382580194d0,  59274488d0,-382813257d0,   !   0568  1.00000
     .-25815479d0, -90484327d0, 382896089d0,  59812684d0,-383421439d0,   !   0569  1.00000
     .-25798689d0, -90768835d0, 383207837d0,  60349675d0,-384024621d0,   !   0570  1.00000
     .-25781520d0, -91052619d0, 383515448d0,  60885463d0,-384622811d0,   !   0571  1.00000
     .-25763970d0, -91335682d0, 383818929d0,  61420047d0,-385216020d0,   !   0572  1.00000
     .-25746041d0, -91618023d0, 384118288d0,  61953429d0,-385804256d0,   !   0573  1.00000
     .-25727734d0, -91899642d0, 384413534d0,  62485610d0,-386387529d0,   !   0574  1.00000
     .-25709049d0, -92180541d0, 384704675d0,  63016590d0,-386965849d0,   !   0575  1.00000
     .-25689988d0, -92460720d0, 384991719d0,  63546370d0,-387539226d0,   !   0576  1.00000
     .-25670551d0, -92740179d0, 385274674d0,  64074952d0,-388107668d0,   !   0577  1.00000
     .-25650739d0, -93018918d0, 385553548d0,  64602335d0,-388671185d0,   !   0578  1.00000
     .-25630553d0, -93296940d0, 385828349d0,  65128522d0,-389229787d0,   !   0579  1.00000
     .-25609993d0, -93574242d0, 386099086d0,  65653512d0,-389783483d0,   !   0580  1.00000
     .-25589061d0, -93850828d0, 386365766d0,  66177306d0,-390332283d0,   !   0581  1.00000
     .-25567757d0, -94126696d0, 386628398d0,  66699906d0,-390876195d0,   !   0582  1.00000
     .-25546082d0, -94401847d0, 386886990d0,  67221312d0,-391415231d0,   !   0583  1.00000
     .-25524037d0, -94676282d0, 387141550d0,  67741524d0,-391949398d0,   !   0584  1.00000
     .-25501623d0, -94950002d0, 387392085d0,  68260545d0,-392478706d0,   !   0585  1.00000
     .-25478840d0, -95223006d0, 387638605d0,  68778374d0,-393003166d0,   !   0586  1.00000
     .-25455690d0, -95495296d0, 387881116d0,  69295013d0,-393522785d0,   !   0587  1.00000
     .-25432172d0, -95766871d0, 388119628d0,  69810462d0,-394037574d0,   !   0588  1.00000
     .-25408289d0, -96037733d0, 388354148d0,  70324722d0,-394547542d0,   !   0589  1.00000
     .-25384040d0, -96307882d0, 388584685d0,  70837794d0,-395052699d0,   !   0590  1.00000
     .-25359427d0, -96577318d0, 388811245d0,  71349679d0,-395553053d0,   !   0591  1.00000
     .-25334450d0, -96846042d0, 389033838d0,  71860378d0,-396048615d0,   !   0592  1.00000
     .-25309110d0, -97114054d0, 389252472d0,  72369891d0,-396539393d0,   !   0593  1.00000
     .-25283408d0, -97381355d0, 389467154d0,  72878219d0,-397025397d0,   !   0594  1.00000
     .-25257345d0, -97647945d0, 389677893d0,  73385364d0,-397506637d0,   !   0595  1.00000
     .-25230921d0, -97913825d0, 389884696d0,  73891326d0,-397983121d0,   !   0596  1.00000
     .-25204138d0, -98178996d0, 390087572d0,  74396105d0,-398454860d0,   !   0597  1.00000
     .-25176996d0, -98443457d0, 390286528d0,  74899703d0,-398921861d0,   !   0598  1.00000
     .-25149496d0, -98707210d0, 390481573d0,  75402121d0,-399384136d0,   !   0599  1.00000
     .-25121638d0, -98970255d0, 390672715d0,  75903359d0,-399841693d0,   !   0600  1.00000
     .-25093424d0, -99232592d0, 390859961d0,  76403418d0,-400294541d0,   !   0601  1.00000
     .-25064855d0, -99494221d0, 391043320d0,  76902299d0,-400742689d0,   !   0602  1.00000
     .-25035931d0, -99755144d0, 391222800d0,  77400003d0,-401186148d0,   !   0603  1.00000
     .-25006652d0,-100015361d0, 391398408d0,  77896530d0,-401624926d0,   !   0604  1.00000
     .-24977021d0,-100274873d0, 391570153d0,  78391883d0,-402059033d0,   !   0605  1.00000
     .-24947036d0,-100533679d0, 391738042d0,  78886060d0,-402488478d0,   !   0606  1.00000
     .-24916701d0,-100791780d0, 391902083d0,  79379064d0,-402913270d0,   !   0607  1.00000
     .-24886014d0,-101049178d0, 392062286d0,  79870894d0,-403333419d0,   !   0608  1.00000
     .-24854977d0,-101305871d0, 392218656d0,  80361553d0,-403748933d0,   !   0609  1.00000
     .-24823591d0,-101561862d0, 392371203d0,  80851040d0,-404159822d0,   !   0610  1.00000
     .-24791857d0,-101817150d0, 392519934d0,  81339357d0,-404566096d0,   !   0611  1.00000
     .-24759775d0,-102071736d0, 392664858d0,  81826504d0,-404967764d0,   !   0612  1.00000
     .-24727346d0,-102325620d0, 392805981d0,  82312482d0,-405364834d0,   !   0613  1.00000
     .-24694571d0,-102578803d0, 392943313d0,  82797293d0,-405757316d0,   !   0614  1.00000
     .-24661450d0,-102831286d0, 393076861d0,  83280936d0,-406145220d0,   !   0615  1.00000
     .-24627985d0,-103083068d0, 393206633d0,  83763413d0,-406528554d0,   !   0616  1.00000
     .-24594177d0,-103334151d0, 393332636d0,  84244724d0,-406907328d0,   !   0617  1.00000
     .-24560025d0,-103584534d0, 393454879d0,  84724871d0,-407281551d0,   !   0618  1.00000
     .-24525531d0,-103834219d0, 393573370d0,  85203854d0,-407651232d0,   !   0619  1.00000
     .-24490696d0,-104083206d0, 393688117d0,  85681674d0,-408016381d0,   !   0620  1.00000
     .-24455521d0,-104331495d0, 393799127d0,  86158332d0,-408377006d0,   !   0621  1.00000
     .-24420005d0,-104579087d0, 393906408d0,  86633829d0,-408733117d0,   !   0622  1.00000
     .-24384151d0,-104825982d0, 394009969d0,  87108165d0,-409084723d0,   !   0623  1.00000
     .-24347958d0,-105072182d0, 394109816d0,  87581342d0,-409431833d0,   !   0624  1.00000
     .-24311428d0,-105317685d0, 394205958d0,  88053360d0,-409774456d0,   !   0625  1.00000
     .-24274562d0,-105562494d0, 394298404d0,  88524220d0,-410112602d0,   !   0626  1.00000
     .-24237359d0,-105806608d0, 394387159d0,  88993922d0,-410446279d0,   !   0627  1.00000
     .-24199822d0,-106050028d0, 394472234d0,  89462469d0,-410775498d0,   !   0628  1.00000
     .-24161950d0,-106292754d0, 394553634d0,  89929860d0,-411100266d0,   !   0629  1.00000
     .-24123745d0,-106534787d0, 394631369d0,  90396096d0,-411420593d0,   !   0630  1.00000
     .-24085207d0,-106776128d0, 394705446d0,  90861179d0,-411736489d0,   !   0631  1.00000
     .-24046337d0,-107016776d0, 394775872d0,  91325109d0,-412047961d0,   !   0632  1.00000
     .-24007136d0,-107256733d0, 394842657d0,  91787887d0,-412355021d0,   !   0633  1.00000
     .-23967605d0,-107495999d0, 394905806d0,  92249513d0,-412657676d0,   !   0634  1.00000
     .-23927744d0,-107734574d0, 394965329d0,  92709989d0,-412955935d0,   !   0635  1.00000
     .-23887554d0,-107972459d0, 395021233d0,  93169316d0,-413249808d0,   !   0636  1.00000
     .-23847036d0,-108209654d0, 395073526d0,  93627493d0,-413539304d0,   !   0637  1.00000
     .-23806191d0,-108446161d0, 395122215d0,  94084523d0,-413824432d0,   !   0638  1.00000
     .-23765020d0,-108681979d0, 395167309d0,  94540406d0,-414105201d0,   !   0639  1.00000
     .-23723523d0,-108917108d0, 395208815d0,  94995142d0,-414381620d0,   !   0640  1.00000
     .-23681701d0,-109151550d0, 395246741d0,  95448733d0,-414653698d0,   !   0641  1.00000
     .-23639555d0,-109385305d0, 395281094d0,  95901180d0,-414921445d0,   !   0642  1.00000
     .-23597085d0,-109618374d0, 395311884d0,  96352482d0,-415184868d0,   !   0643  1.00000
     .-23554293d0,-109850756d0, 395339116d0,  96802642d0,-415443978d0,   !   0644  1.00000
     .-23511179d0,-110082453d0, 395362799d0,  97251660d0,-415698784d0,   !   0645  1.00000
     .-23467744d0,-110313464d0, 395382942d0,  97699536d0,-415949293d0,   !   0646  1.00000
     .-23423989d0,-110543791d0, 395399550d0,  98146272d0,-416195516d0,   !   0647  1.00000
     .-23379915d0,-110773434d0, 395412633d0,  98591869d0,-416437462d0,   !   0648  1.00000
     .-23335522d0,-111002394d0, 395422198d0,  99036327d0,-416675139d0,   !   0649  1.00000
     .-23290810d0,-111230670d0, 395428253d0,  99479647d0,-416908556d0,   !   0650  1.00000
     .-23245782d0,-111458264d0, 395430805d0,  99921830d0,-417137722d0,   !   0651  1.00000
     .-23200437d0,-111685175d0, 395429862d0, 100362876d0,-417362648d0,   !   0652  1.00000
     .-23154777d0,-111911405d0, 395425431d0, 100802788d0,-417583340d0,   !   0653  1.00000
     .-23108802d0,-112136954d0, 395417522d0, 101241564d0,-417799809d0,   !   0654  1.00000
     .-23062513d0,-112361823d0, 395406140d0, 101679207d0,-418012063d0,   !   0655  1.00000
     .-23015910d0,-112586011d0, 395391295d0, 102115717d0,-418220111d0,   !   0656  1.00000
     .-22968995d0,-112809520d0, 395372993d0, 102551095d0,-418423963d0,   !   0657  1.00000
     .-22921768d0,-113032350d0, 395351242d0, 102985342d0,-418623627d0,   !   0658  1.00000
     .-22874230d0,-113254501d0, 395326050d0, 103418459d0,-418819112d0,   !   0659  1.00000
     .-22826382d0,-113475974d0, 395297425d0, 103850446d0,-419010428d0,   !   0660  1.00000
     .-22778225d0,-113696769d0, 395265374d0, 104281304d0,-419197582d0,   !   0661  1.00000
     .-22729759d0,-113916888d0, 395229905d0, 104711034d0,-419380584d0,   !   0662  1.00000
     .-22680984d0,-114136330d0, 395191025d0, 105139638d0,-419559443d0,   !   0663  1.00000
     .-22631903d0,-114355096d0, 395148743d0, 105567115d0,-419734168d0,   !   0664  1.00000
     .-22582515d0,-114573186d0, 395103065d0, 105993467d0,-419904768d0,   !   0665  1.00000
     .-22532822d0,-114790601d0, 395054000d0, 106418694d0,-420071251d0,   !   0666  1.00000
     .-22482824d0,-115007342d0, 395001555d0, 106842798d0,-420233627d0,   !   0667  1.00000
     .-22432521d0,-115223409d0, 394945738d0, 107265778d0,-420391904d0,   !   0668  1.00000
     .-22381916d0,-115438802d0, 394886556d0, 107687637d0,-420546091d0,   !   0669  1.00000
     .-22331007d0,-115653522d0, 394824017d0, 108108375d0,-420696197d0,   !   0670  1.00000
     .-22279797d0,-115867569d0, 394758128d0, 108527992d0,-420842232d0,   !   0671  1.00000
     .-22228286d0,-116080945d0, 394688898d0, 108946490d0,-420984203d0,   !   0672  1.00000
     .-22176475d0,-116293649d0, 394616333d0, 109363869d0,-421122120d0,   !   0673  1.00000
     .-22124363d0,-116505682d0, 394540442d0, 109780130d0,-421255992d0,   !   0674  1.00000
     .-22071954d0,-116717044d0, 394461231d0, 110195274d0,-421385827d0,   !   0675  1.00000
     .-22019246d0,-116927736d0, 394378709d0, 110609302d0,-421511634d0,   !   0676  1.00000
     .-21966241d0,-117137759d0, 394292883d0, 111022215d0,-421633422d0,   !   0677  1.00000
     .-21912939d0,-117347113d0, 394203761d0, 111434013d0,-421751200d0,   !   0678  1.00000
     .-21859341d0,-117555798d0, 394111350d0, 111844697d0,-421864977d0,   !   0679  1.00000
     .-21805449d0,-117763815d0, 394015657d0, 112254269d0,-421974761d0,   !   0680  1.00000
     .-21751262d0,-117971165d0, 393916691d0, 112662729d0,-422080561d0,   !   0681  1.00000
     .-21696782d0,-118177848d0, 393814458d0, 113070077d0,-422182387d0,   !   0682  1.00000
     .-21642009d0,-118383864d0, 393708967d0, 113476315d0,-422280246d0,   !   0683  1.00000
     .-21586944d0,-118589214d0, 393600225d0, 113881444d0,-422374149d0,   !   0684  1.00000
     .-21531588d0,-118793898d0, 393488239d0, 114285464d0,-422464102d0,   !   0685  1.00000
     .-21475941d0,-118997918d0, 393373017d0, 114688376d0,-422550116d0,   !   0686  1.00000
     .-21420005d0,-119201273d0, 393254566d0, 115090182d0,-422632199d0,   !   0687  1.00000
     .-21363779d0,-119403964d0, 393132894d0, 115490881d0,-422710359d0,   !   0688  1.00000
     .-21307266d0,-119605991d0, 393008009d0, 115890475d0,-422784606d0,   !   0689  1.00000
     .-21250465d0,-119807356d0, 392879918d0, 116288964d0,-422854948d0,   !   0690  1.00000
     .-21193377d0,-120008058d0, 392748628d0, 116686350d0,-422921394d0,   !   0691  1.00000
     .-21136003d0,-120208097d0, 392614146d0, 117082633d0,-422983952d0,   !   0692  1.00000
     .-21078344d0,-120407476d0, 392476482d0, 117477814d0,-423042632d0,   !   0693  1.00000
     .-21020400d0,-120606193d0, 392335641d0, 117871893d0,-423097442d0,   !   0694  1.00000
     .-20962173d0,-120804250d0, 392191631d0, 118264873d0,-423148391d0,   !   0695  1.00000
     .-20903663d0,-121001647d0, 392044460d0, 118656753d0,-423195487d0,   !   0696  1.00000
     .-20844870d0,-121198384d0, 391894136d0, 119047534d0,-423238740d0,   !   0697  1.00000
     .-20785796d0,-121394463d0, 391740665d0, 119437218d0,-423278157d0,   !   0698  1.00000
     .-20726441d0,-121589883d0, 391584055d0, 119825805d0,-423313748d0,   !   0699  1.00000
     .-20666806d0,-121784644d0, 391424314d0, 120213295d0,-423345521d0,   !   0700  1.00000
     .-20606892d0,-121978749d0, 391261449d0, 120599691d0,-423373485d0,   !   0701  1.00000
     .-20546699d0,-122172196d0, 391095467d0, 120984991d0,-423397649d0,   !   0702  1.00000
     .-20486229d0,-122364987d0, 390926376d0, 121369199d0,-423418021d0,   !   0703  1.00000
     .-20425481d0,-122557122d0, 390754183d0, 121752313d0,-423434609d0,   !   0704  1.00000
     .-20364458d0,-122748601d0, 390578896d0, 122134336d0,-423447424d0,   !   0705  1.00000
     .-20303158d0,-122939426d0, 390400522d0, 122515267d0,-423456472d0,   !   0706  1.00000
     .-20241584d0,-123129596d0, 390219068d0, 122895108d0,-423461763d0,   !   0707  1.00000
     .-20179736d0,-123319111d0, 390034543d0, 123273860d0,-423463305d0,   !   0708  1.00000
     .-20117614d0,-123507974d0, 389846952d0, 123651523d0,-423461108d0,   !   0709  1.00000
     .-20055220d0,-123696183d0, 389656304d0, 124028098d0,-423455179d0,   !   0710  1.00000
     .-19992554d0,-123883740d0, 389462607d0, 124403587d0,-423445527d0,   !   0711  1.00000
     .-19929617d0,-124070645d0, 389265866d0, 124777989d0,-423432161d0,   !   0712  1.00000
     .-19866409d0,-124256899d0, 389066091d0, 125151306d0,-423415090d0,   !   0713  1.00000
     .-19802932d0,-124442501d0, 388863287d0, 125523539d0,-423394321d0,   !   0714  1.00000
     .-19739186d0,-124627453d0, 388657463d0, 125894688d0,-423369864d0,   !   0715  1.00000
     .-19675171d0,-124811755d0, 388448626d0, 126264754d0,-423341727d0,   !   0716  1.00000
     .-19610889d0,-124995408d0, 388236783d0, 126633738d0,-423309919d0,   !   0717  1.00000
     .-19546341d0,-125178411d0, 388021941d0, 127001642d0,-423274448d0,   !   0718  1.00000
     .-19481526d0,-125360766d0, 387804108d0, 127368465d0,-423235323d0,   !   0719  1.00000
     .-19416447d0,-125542474d0, 387583291d0, 127734208d0,-423192552d0,   !   0720  1.00000
     .-19351102d0,-125723533d0, 387359498d0, 128098873d0,-423146144d0,   !   0721  1.00000
     .-19285494d0,-125903946d0, 387132736d0, 128462461d0,-423096107d0,   !   0722  1.00000
     .-19219623d0,-126083712d0, 386903011d0, 128824971d0,-423042451d0,   !   0723  1.00000
     .-19153490d0,-126262832d0, 386670332d0, 129186406d0,-422985182d0,   !   0724  1.00000
     .-19087095d0,-126441307d0, 386434706d0, 129546765d0,-422924311d0,   !   0725  1.00000
     .-19020439d0,-126619137d0, 386196140d0, 129906049d0,-422859845d0,   !   0726  1.00000
     .-18953523d0,-126796323d0, 385954641d0, 130264261d0,-422791793d0,   !   0727  1.00000
     .-18886347d0,-126972864d0, 385710216d0, 130621399d0,-422720164d0,   !   0728  1.00000
     .-18818913d0,-127148762d0, 385462874d0, 130977466d0,-422644965d0,   !   0729  1.00000
     .-18751221d0,-127324017d0, 385212620d0, 131332461d0,-422566206d0,   !   0730  1.00000
     .-18683272d0,-127498630d0, 384959463d0, 131686386d0,-422483894d0,   !   0731  1.00000
     .-18615066d0,-127672600d0, 384703409d0, 132039242d0,-422398039d0,   !   0732  1.00000
     .-18546604d0,-127845929d0, 384444466d0, 132391030d0,-422308649d0,   !   0733  1.00000
     .-18477887d0,-128018618d0, 384182642d0, 132741749d0,-422215732d0,   !   0734  1.00000
     .-18408916d0,-128190665d0, 383917942d0, 133091402d0,-422119296d0,   !   0735  1.00000
     .-18339692d0,-128362073d0, 383650375d0, 133439989d0,-422019351d0,   !   0736  1.00000
     .-18270214d0,-128532842d0, 383379949d0, 133787510d0,-421915904d0,   !   0737  1.00000
     .-18200485d0,-128702971d0, 383106669d0, 134133967d0,-421808965d0,   !   0738  1.00000
     .-18130504d0,-128872462d0, 382830543d0, 134479361d0,-421698540d0,   !   0739  1.00000
     .-18060272d0,-129041315d0, 382551579d0, 134823692d0,-421584639d0,   !   0740  1.00000
     .-17989790d0,-129209531d0, 382269783d0, 135166961d0,-421467271d0,   !   0741  1.00000
     .-17919059d0,-129377109d0, 381985164d0, 135509169d0,-421346443d0,   !   0742  1.00000
     .-17848080d0,-129544051d0, 381697727d0, 135850317d0,-421222164d0,   !   0743  1.00000
     .-17776852d0,-129710358d0, 381407481d0, 136190405d0,-421094442d0,   !   0744  1.00000
     .-17705378d0,-129876029d0, 381114432d0, 136529435d0,-420963286d0,   !   0745  1.00000
     .-17633657d0,-130041065d0, 380818587d0, 136867407d0,-420828705d0,   !   0746  1.00000
     .-17561690d0,-130205466d0, 380519955d0, 137204323d0,-420690705d0,   !   0747  1.00000
     .-17489479d0,-130369234d0, 380218541d0, 137540182d0,-420549297d0,   !   0748  1.00000
     .-17417023d0,-130532368d0, 379914354d0, 137874986d0,-420404488d0,   !   0749  1.00000
     .-17344324d0,-130694869d0, 379607400d0, 138208736d0,-420256286d0,   !   0750  1.00000
     .-17271381d0,-130856738d0, 379297686d0, 138541433d0,-420104700d0,   !   0751  1.00000
     .-17198197d0,-131017975d0, 378985220d0, 138873076d0,-419949738d0,   !   0752  1.00000
     .-17124772d0,-131178580d0, 378670008d0, 139203668d0,-419791409d0,   !   0753  1.00000
     .-17051105d0,-131338555d0, 378352059d0, 139533209d0,-419629720d0,   !   0754  1.00000
     .-16977199d0,-131497899d0, 378031378d0, 139861700d0,-419464681d0,   !   0755  1.00000
     .-16903054d0,-131656613d0, 377707973d0, 140189141d0,-419296299d0,   !   0756  1.00000
     .-16828670d0,-131814698d0, 377381852d0, 140515534d0,-419124583d0,   !   0757  1.00000
     .-16754048d0,-131972154d0, 377053021d0, 140840879d0,-418949541d0,   !   0758  1.00000
     .-16679189d0,-132128981d0, 376721487d0, 141165177d0,-418771181d0,   !   0759  1.00000
     .-16604094d0,-132285181d0, 376387258d0, 141488430d0,-418589512d0,   !   0760  1.00000
     .-16528763d0,-132440753d0, 376050341d0, 141810637d0,-418404542d0,   !   0761  1.00000
     .-16453197d0,-132595698d0, 375710742d0, 142131800d0,-418216278d0,   !   0762  1.00000
     .-16377396d0,-132750017d0, 375368469d0, 142451919d0,-418024731d0,   !   0763  1.00000
     .-16301362d0,-132903710d0, 375023529d0, 142770996d0,-417829906d0,   !   0764  1.00000
     .-16225096d0,-133056777d0, 374675928d0, 143089031d0,-417631814d0,   !   0765  1.00000
     .-16148597d0,-133209220d0, 374325675d0, 143406025d0,-417430462d0,   !   0766  1.00000
     .-16071867d0,-133361038d0, 373972776d0, 143721978d0,-417225858d0,   !   0767  1.00000
     .-15994906d0,-133512233d0, 373617238d0, 144036893d0,-417018011d0,   !   0768  1.00000
     .-15917714d0,-133662804d0, 373259068d0, 144350769d0,-416806928d0,   !   0769  1.00000
     .-15840294d0,-133812752d0, 372898273d0, 144663607d0,-416592619d0,   !   0770  1.00000
     .-15762645d0,-133962077d0, 372534860d0, 144975409d0,-416375091d0,   !   0771  1.00000
     .-15684768d0,-134110781d0, 372168837d0, 145286175d0,-416154352d0,   !   0772  1.00000
     .-15606664d0,-134258863d0, 371800210d0, 145595905d0,-415930410d0,   !   0773  1.00000
     .-15528333d0,-134406325d0, 371428986d0, 145904601d0,-415703275d0,   !   0774  1.00000
     .-15449777d0,-134553166d0, 371055173d0, 146212264d0,-415472953d0,   !   0775  1.00000
     .-15370995d0,-134699387d0, 370678776d0, 146518894d0,-415239454d0,   !   0776  1.00000
     .-15291989d0,-134844989d0, 370299805d0, 146824492d0,-415002785d0,   !   0777  1.00000
     .-15212759d0,-134989972d0, 369918264d0, 147129060d0,-414762954d0,   !   0778  1.00000
     .-15133307d0,-135134337d0, 369534162d0, 147432597d0,-414519970d0,   !   0779  1.00000
     .-15053632d0,-135278083d0, 369147505d0, 147735105d0,-414273841d0,   !   0780  1.00000
     .-14973735d0,-135421212d0, 368758301d0, 148036585d0,-414024575d0,   !   0781  1.00000
     .-14893617d0,-135563725d0, 368366555d0, 148337036d0,-413772179d0,   !   0782  1.00000
     .-14813280d0,-135705621d0, 367972276d0, 148636462d0,-413516663d0,   !   0783  1.00000
     .-14732722d0,-135846901d0, 367575470d0, 148934861d0,-413258034d0,   !   0784  1.00000
     .-14651946d0,-135987566d0, 367176145d0, 149232235d0,-412996301d0,   !   0785  1.00000
     .-14570952d0,-136127616d0, 366774306d0, 149528585d0,-412731471d0,   !   0786  1.00000
     .-14489740d0,-136267051d0, 366369961d0, 149823911d0,-412463552d0,   !   0787  1.00000
     .-14408312d0,-136405873d0, 365963118d0, 150118215d0,-412192554d0,   !   0788  1.00000
     .-14326667d0,-136544081d0, 365553782d0, 150411497d0,-411918483d0,   !   0789  1.00000
     .-14244807d0,-136681677d0, 365141961d0, 150703758d0,-411641347d0,   !   0790  1.00000
     .-14162733d0,-136818660d0, 364727662d0, 150994999d0,-411361156d0,   !   0791  1.00000
     .-14080444d0,-136955031d0, 364310891d0, 151285221d0,-411077917d0,   !   0792  1.00000
     .-13997942d0,-137090791d0, 363891657d0, 151574424d0,-410791638d0,   !   0793  1.00000
     .-13915228d0,-137225940d0, 363469964d0, 151862610d0,-410502327d0,   !   0794  1.00000
     .-13832302d0,-137360479d0, 363045821d0, 152149779d0,-410209992d0,   !   0795  1.00000
     .-13749164d0,-137494408d0, 362619235d0, 152435932d0,-409914641d0,   !   0796  1.00000
     .-13665816d0,-137627727d0, 362190211d0, 152721070d0,-409616282d0,   !   0797  1.00000
     .-13582259d0,-137760438d0, 361758758d0, 153005194d0,-409314924d0,   !   0798  1.00000
     .-13498492d0,-137892540d0, 361324881d0, 153288304d0,-409010574d0,   !   0799  1.00000
     .-13414516d0,-138024034d0, 360888589d0, 153570402d0,-408703240d0,   !   0800  1.00000
     .-13330333d0,-138154922d0, 360449887d0, 153851488d0,-408392931d0,   !   0801  1.00000
     .-13245943d0,-138285202d0, 360008782d0, 154131563d0,-408079654d0,   !   0802  1.00000
     .-13161346d0,-138414876d0, 359565282d0, 154410628d0,-407763417d0,   !   0803  1.00000
     .-13076544d0,-138543944d0, 359119393d0, 154688684d0,-407444228d0,   !   0804  1.00000
     .-12991536d0,-138672406d0, 358671123d0, 154965731d0,-407122096d0,   !   0805  1.00000
     .-12906324d0,-138800264d0, 358220477d0, 155241771d0,-406797028d0,   !   0806  1.00000
     .-12820909d0,-138927518d0, 357767463d0, 155516804d0,-406469032d0,   !   0807  1.00000
     .-12735291d0,-139054167d0, 357312087d0, 155790832d0,-406138116d0,   !   0808  1.00000
     .-12649470d0,-139180214d0, 356854358d0, 156063854d0,-405804288d0,   !   0809  1.00000
     .-12563448d0,-139305657d0, 356394280d0, 156335872d0,-405467557d0,   !   0810  1.00000
     .-12477224d0,-139430498d0, 355931861d0, 156606887d0,-405127929d0,   !   0811  1.00000
     .-12390801d0,-139554738d0, 355467108d0, 156876899d0,-404785413d0,   !   0812  1.00000
     .-12304178d0,-139678376d0, 355000028d0, 157145909d0,-404440017d0,   !   0813  1.00000
     .-12217356d0,-139801413d0, 354530627d0, 157413918d0,-404091749d0,   !   0814  1.00000
     .-12130336d0,-139923850d0, 354058912d0, 157680928d0,-403740617d0,   !   0815  1.00000
     .-12043118d0,-140045687d0, 353584891d0, 157946938d0,-403386628d0,   !   0816  1.00000
     .-11955704d0,-140166925d0, 353108569d0, 158211950d0,-403029790d0,   !   0817  1.00000
     .-11868093d0,-140287563d0, 352629954d0, 158475964d0,-402670113d0,   !   0818  1.00000
     .-11780287d0,-140407604d0, 352149052d0, 158738982d0,-402307602d0,   !   0819  1.00000
     .-11692286d0,-140527047d0, 351665870d0, 159001004d0,-401942266d0,   !   0820  1.00000
     .-11604091d0,-140645892d0, 351180415d0, 159262031d0,-401574114d0,   !   0821  1.00000
     .-11515702d0,-140764141d0, 350692693d0, 159522063d0,-401203152d0,   !   0822  1.00000
     .-11427120d0,-140881793d0, 350202712d0, 159781102d0,-400829389d0,   !   0823  1.00000
     .-11338347d0,-140998850d0, 349710478d0, 160039149d0,-400452833d0,   !   0824  1.00000
     .-11249382d0,-141115311d0, 349215997d0, 160296204d0,-400073491d0,   !   0825  1.00000
     .-11160226d0,-141231177d0, 348719278d0, 160552268d0,-399691372d0,   !   0826  1.00000
     .-11070880d0,-141346449d0, 348220325d0, 160807342d0,-399306483d0,   !   0827  1.00000
     .-10981345d0,-141461127d0, 347719146d0, 161061427d0,-398918831d0,   !   0828  1.00000
     .-10891621d0,-141575212d0, 347215748d0, 161314524d0,-398528426d0,   !   0829  1.00000
     .-10801709d0,-141688705d0, 346710137d0, 161566633d0,-398135274d0,   !   0830  1.00000
     .-10711609d0,-141801604d0, 346202321d0, 161817755d0,-397739383d0,   !   0831  1.00000
     .-10621323d0,-141913912d0, 345692305d0, 162067892d0,-397340762d0,   !   0832  1.00000
     .-10530850d0,-142025629d0, 345180096d0, 162317043d0,-396939418d0,   !   0833  1.00000
     .-10440192d0,-142136755d0, 344665702d0, 162565211d0,-396535358d0,   !   0834  1.00000
     .-10349350d0,-142247291d0, 344149128d0, 162812395d0,-396128592d0,   !   0835  1.00000
     .-10258323d0,-142357237d0, 343630382d0, 163058596d0,-395719125d0,   !   0836  1.00000
     .-10167113d0,-142466593d0, 343109470d0, 163303816d0,-395306967d0,   !   0837  1.00000
     .-10075720d0,-142575361d0, 342586398d0, 163548055d0,-394892124d0,   !   0838  1.00000
     . -9984145d0,-142683541d0, 342061174d0, 163791314d0,-394474605d0,   !   0839  1.00000
     . -9892389d0,-142791132d0, 341533804d0, 164033594d0,-394054418d0,   !   0840  1.00000
     . -9800452d0,-142898137d0, 341004295d0, 164274896d0,-393631569d0,   !   0841  1.00000
     . -9708335d0,-143004554d0, 340472653d0, 164515220d0,-393206068d0,   !   0842  1.00000
     . -9616038d0,-143110386d0, 339938885d0, 164754568d0,-392777921d0,   !   0843  1.00000
     . -9523562d0,-143215631d0, 339402998d0, 164992939d0,-392347137d0,   !   0844  1.00000
     . -9430909d0,-143320291d0, 338864997d0, 165230336d0,-391913722d0,   !   0845  1.00000
     . -9338078d0,-143424367d0, 338324891d0, 165466759d0,-391477685d0,   !   0846  1.00000
     . -9245070d0,-143527858d0, 337782685d0, 165702208d0,-391039034d0,   !   0847  1.00000
     . -9151887d0,-143630766d0, 337238386d0, 165936685d0,-390597776d0,   !   0848  1.00000
     . -9058527d0,-143733090d0, 336692001d0, 166170190d0,-390153919d0,   !   0849  1.00000
     . -8964993d0,-143834831d0, 336143535d0, 166402724d0,-389707470d0,   !   0850  1.00000
     . -8871285d0,-143935990d0, 335592997d0, 166634288d0,-389258437d0,   !   0851  1.00000
     . -8777403d0,-144036568d0, 335040392d0, 166864884d0,-388806829d0,   !   0852  1.00000
     . -8683349d0,-144136564d0, 334485727d0, 167094511d0,-388352651d0,   !   0853  1.00000
     . -8589122d0,-144235979d0, 333929008d0, 167323170d0,-387895913d0,   !   0854  1.00000
     . -8494724d0,-144334814d0, 333370243d0, 167550863d0,-387436622d0,   !   0855  1.00000
     . -8400155d0,-144433070d0, 332809437d0, 167777590d0,-386974786d0,   !   0856  1.00000
     . -8305415d0,-144530746d0, 332246597d0, 168003352d0,-386510411d0,   !   0857  1.00000
     . -8210507d0,-144627843d0, 331681730d0, 168228150d0,-386043506d0,   !   0858  1.00000
     . -8115429d0,-144724362d0, 331114842d0, 168451985d0,-385574079d0,   !   0859  1.00000
     . -8020183d0,-144820303d0, 330545940d0, 168674857d0,-385102136d0,   !   0860  1.00000
     . -7924769d0,-144915667d0, 329975030d0, 168896768d0,-384627686d0,   !   0861  1.00000
     . -7829189d0,-145010454d0, 329402119d0, 169117718d0,-384150737d0,   !   0862  1.00000
     . -7733442d0,-145104665d0, 328827213d0, 169337708d0,-383671295d0,   !   0863  1.00000
     . -7637530d0,-145198301d0, 328250319d0, 169556739d0,-383189368d0,   !   0864  1.00000
     . -7541452d0,-145291361d0, 327671444d0, 169774811d0,-382704965d0,   !   0865  1.00000
     . -7445211d0,-145383846d0, 327090593d0, 169991926d0,-382218092d0,   !   0866  1.00000
     . -7348805d0,-145475757d0, 326507774d0, 170208085d0,-381728757d0,   !   0867  1.00000
     . -7252237d0,-145567094d0, 325922992d0, 170423288d0,-381236968d0,   !   0868  1.00000
     . -7155506d0,-145657858d0, 325336255d0, 170637536d0,-380742732d0,   !   0869  1.00000
     . -7058613d0,-145748050d0, 324747569d0, 170850830d0,-380246057d0,   !   0870  1.00000
     . -6961560d0,-145837669d0, 324156940d0, 171063170d0,-379746950d0,   !   0871  1.00000
     . -6864346d0,-145926716d0, 323564374d0, 171274558d0,-379245419d0,   !   0872  1.00000
     . -6766972d0,-146015192d0, 322969879d0, 171484995d0,-378741472d0,   !   0873  1.00000
     . -6669439d0,-146103097d0, 322373460d0, 171694480d0,-378235115d0,   !   0874  1.00000
     . -6571747d0,-146190432d0, 321775125d0, 171903016d0,-377726357d0,   !   0875  1.00000
     . -6473898d0,-146277198d0, 321174878d0, 172110603d0,-377215205d0,   !   0876  1.00000
     . -6375891d0,-146363394d0, 320572728d0, 172317242d0,-376701666d0,   !   0877  1.00000
     . -6277728d0,-146449021d0, 319968681d0, 172522933d0,-376185749d0,   !   0878  1.00000
     . -6179409d0,-146534080d0, 319362742d0, 172727677d0,-375667459d0,   !   0879  1.00000
     . -6080935d0,-146618572d0, 318754918d0, 172931476d0,-375146806d0,   !   0880  1.00000
     . -5982306d0,-146702496d0, 318145216d0, 173134330d0,-374623796d0,   !   0881  1.00000
     . -5883523d0,-146785853d0, 317533642d0, 173336240d0,-374098437d0,   !   0882  1.00000
     . -5784586d0,-146868645d0, 316920202d0, 173537207d0,-373570736d0,   !   0883  1.00000
     . -5685497d0,-146950870d0, 316304903d0, 173737231d0,-373040701d0,   !   0884  1.00000
     . -5586256d0,-147032531d0, 315687751d0, 173936314d0,-372508340d0,   !   0885  1.00000
     . -5486863d0,-147113626d0, 315068754d0, 174134455d0,-371973659d0,   !   0886  1.00000
     . -5387320d0,-147194158d0, 314447916d0, 174331658d0,-371436666d0,   !   0887  1.00000
     . -5287626d0,-147274125d0, 313825244d0, 174527921d0,-370897369d0,   !   0888  1.00000
     . -5187783d0,-147353530d0, 313200745d0, 174723245d0,-370355775d0,   !   0889  1.00000
     . -5087791d0,-147432372d0, 312574425d0, 174917633d0,-369811891d0,   !   0890  1.00000
     . -4987650d0,-147510651d0, 311946291d0, 175111084d0,-369265725d0,   !   0891  1.00000
     . -4887362d0,-147588369d0, 311316348d0, 175303599d0,-368717284d0,   !   0892  1.00000
     . -4786927d0,-147665526d0, 310684604d0, 175495179d0,-368166576d0,   !   0893  1.00000
     . -4686346d0,-147742121d0, 310051064d0, 175685825d0,-367613608d0,   !   0894  1.00000
     . -4585619d0,-147818157d0, 309415735d0, 175875538d0,-367058387d0,   !   0895  1.00000
     . -4484747d0,-147893633d0, 308778623d0, 176064319d0,-366500921d0,   !   0896  1.00000
     . -4383730d0,-147968550d0, 308139734d0, 176252168d0,-365941218d0,   !   0897  1.00000
     . -4282570d0,-148042908d0, 307499076d0, 176439087d0,-365379283d0,   !   0898  1.00000
     . -4181266d0,-148116708d0, 306856653d0, 176625075d0,-364815126d0,   !   0899  1.00000
     . -4079820d0,-148189950d0, 306212473d0, 176810135d0,-364248753d0,   !   0900  1.00000
     . -3978232d0,-148262635d0, 305566541d0, 176994266d0,-363680172d0,   !   0901  1.00000
     . -3876503d0,-148334764d0, 304918865d0, 177177470d0,-363109389d0,   !   0902  1.00000
     . -3774633d0,-148406336d0, 304269450d0, 177359748d0,-362536413d0,   !   0903  1.00000
     . -3672623d0,-148477352d0, 303618302d0, 177541099d0,-361961250d0,   !   0904  1.00000
     . -3570473d0,-148547813d0, 302965428d0, 177721526d0,-361383908d0,   !   0905  1.00000
     . -3468185d0,-148617720d0, 302310834d0, 177901029d0,-360804395d0,   !   0906  1.00000
     . -3365758d0,-148687072d0, 301654527d0, 178079609d0,-360222717d0,   !   0907  1.00000
     . -3263194d0,-148755871d0, 300996512d0, 178257266d0,-359638882d0,   !   0908  1.00000
     . -3160493d0,-148824117d0, 300336796d0, 178434002d0,-359052897d0,   !   0909  1.00000
     . -3057656d0,-148891809d0, 299675386d0, 178609818d0,-358464769d0,   !   0910  1.00000
     . -2954683d0,-148958950d0, 299012287d0, 178784713d0,-357874506d0,   !   0911  1.00000
     . -2851575d0,-149025539d0, 298347505d0, 178958689d0,-357282115d0,   !   0912  1.00000
     . -2748333d0,-149091577d0, 297681047d0, 179131748d0,-356687603d0,   !   0913  1.00000
     . -2644957d0,-149157064d0, 297012919d0, 179303888d0,-356090978d0,   !   0914  1.00000
     . -2541447d0,-149222001d0, 296343128d0, 179475113d0,-355492246d0,   !   0915  1.00000
     . -2437805d0,-149286388d0, 295671679d0, 179645422d0,-354891415d0,   !   0916  1.00000
     . -2334032d0,-149350226d0, 294998579d0, 179814815d0,-354288493d0,   !   0917  1.00000
     . -2230127d0,-149413516d0, 294323834d0, 179983295d0,-353683486d0,   !   0918  1.00000
     . -2126091d0,-149476257d0, 293647449d0, 180150862d0,-353076402d0,   !   0919  1.00000
     . -2021925d0,-149538451d0, 292969433d0, 180317517d0,-352467247d0,   !   0920  1.00000
     . -1917630d0,-149600097d0, 292289789d0, 180483260d0,-351856030d0,   !   0921  1.00000
     . -1813206d0,-149661197d0, 291608526d0, 180648092d0,-351242757d0,   !   0922  1.00000
     . -1708653d0,-149721751d0, 290925648d0, 180812014d0,-350627436d0,   !   0923  1.00000
     . -1603973d0,-149781759d0, 290241162d0, 180975028d0,-350010073d0,   !   0924  1.00000
     . -1499166d0,-149841221d0, 289555075d0, 181137133d0,-349390676d0,   !   0925  1.00000
     . -1394233d0,-149900140d0, 288867391d0, 181298331d0,-348769253d0,   !   0926  1.00000
     . -1289174d0,-149958514d0, 288178119d0, 181458623d0,-348145809d0,   !   0927  1.00000
     . -1183990d0,-150016344d0, 287487263d0, 181618009d0,-347520353d0,   !   0928  1.00000
     . -1078681d0,-150073631d0, 286794829d0, 181776490d0,-346892892d0,   !   0929  1.00000
     .  -973249d0,-150130376d0, 286100825d0, 181934067d0,-346263432d0,   !   0930  1.00000
     .  -867693d0,-150186578d0, 285405256d0, 182090742d0,-345631981d0,   !   0931  1.00000
     .  -762014d0,-150242239d0, 284708128d0, 182246513d0,-344998546d0,   !   0932  1.00000
     .  -656214d0,-150297358d0, 284009447d0, 182401384d0,-344363135d0,   !   0933  1.00000
     .  -550292d0,-150351937d0, 283309219d0, 182555353d0,-343725754d0,   !   0934  1.00000
     .  -444249d0,-150405976d0, 282607451d0, 182708423d0,-343086410d0,   !   0935  1.00000
     .  -338086d0,-150459475d0, 281904149d0, 182860594d0,-342445110d0,   !   0936  1.00000
     .  -231803d0,-150512435d0, 281199319d0, 183011867d0,-341801863d0,   !   0937  1.00000
     .  -125401d0,-150564857d0, 280492966d0, 183162243d0,-341156674d0,   !   0938  1.00000
     .   -18881d0,-150616740d0, 279785097d0, 183311722d0,-340509550d0,   !   0939  1.00000
     .    87756d0,-150668086d0, 279075718d0, 183460305d0,-339860500d0,   !   0940  1.00000
     .   194511d0,-150718894d0, 278364835d0, 183607994d0,-339209530d0,   !   0941  1.00000
     .   301382d0,-150769166d0, 277652454d0, 183754789d0,-338556647d0,   !   0942  1.00000
     .   408370d0,-150818902d0, 276938582d0, 183900691d0,-337901858d0,   !   0943  1.00000
     .   515473d0,-150868102d0, 276223223d0, 184045700d0,-337245170d0,   !   0944  1.00000
     .   622691d0,-150916767d0, 275506385d0, 184189818d0,-336586590d0,   !   0945  1.00000
     .   730023d0,-150964897d0, 274788073d0, 184333046d0,-335926126d0,   !   0946  1.00000
     .   837468d0,-151012493d0, 274068294d0, 184475383d0,-335263784d0,   !   0947  1.00000
     .   945027d0,-151059556d0, 273347053d0, 184616832d0,-334599572d0,   !   0948  1.00000
     .  1052698d0,-151106085d0, 272624356d0, 184757393d0,-333933496d0,   !   0949  1.00000
     .  1160480d0,-151152082d0, 271900210d0, 184897066d0,-333265563d0,   !   0950  1.00000
     .  1268374d0,-151197547d0, 271174620d0, 185035853d0,-332595781d0,   !   0951  1.00000
     .  1376378d0,-151242480d0, 270447593d0, 185173754d0,-331924156d0,   !   0952  1.00000
     .  1484492d0,-151286882d0, 269719134d0, 185310771d0,-331250696d0,   !   0953  1.00000
     .  1592716d0,-151330753d0, 268989250d0, 185446904d0,-330575407d0,   !   0954  1.00000
     .  1701048d0,-151374095d0, 268257946d0, 185582153d0,-329898297d0,   !   0955  1.00000
     .  1809489d0,-151416906d0, 267525228d0, 185716521d0,-329219372d0,   !   0956  1.00000
     .  1918037d0,-151459189d0, 266791103d0, 185850007d0,-328538639d0,   !   0957  1.00000
     .  2026691d0,-151500943d0, 266055577d0, 185982612d0,-327856106d0,   !   0958  1.00000
     .  2135453d0,-151542169d0, 265318655d0, 186114338d0,-327171779d0,   !   0959  1.00000
     .  2244320d0,-151582867d0, 264580343d0, 186245185d0,-326485666d0,   !   0960  1.00000
     .  2353292d0,-151623038d0, 263840647d0, 186375154d0,-325797773d0,   !   0961  1.00000
     .  2462369d0,-151662682d0, 263099574d0, 186504245d0,-325108107d0,   !   0962  1.00000
     .  2571549d0,-151701801d0, 262357129d0, 186632460d0,-324416675d0,   !   0963  1.00000
     .  2680833d0,-151740394d0, 261613318d0, 186759800d0,-323723484d0,   !   0964  1.00000
     .  2790220d0,-151778461d0, 260868147d0, 186886265d0,-323028541d0,   !   0965  1.00000
     .  2899709d0,-151816004d0, 260121622d0, 187011856d0,-322331853d0,   !   0966  1.00000
     .  3009300d0,-151853024d0, 259373749d0, 187136573d0,-321633427d0,   !   0967  1.00000
     .  3118991d0,-151889519d0, 258624534d0, 187260419d0,-320933269d0,   !   0968  1.00000
     .  3228783d0,-151925491d0, 257873983d0, 187383393d0,-320231387d0,   !   0969  1.00000
     .  3338675d0,-151960941d0, 257122102d0, 187505497d0,-319527788d0,   !   0970  1.00000
     .  3448666d0,-151995869d0, 256368896d0, 187626731d0,-318822478d0,   !   0971  1.00000
     .  3558756d0,-152030275d0, 255614371d0, 187747096d0,-318115464d0,   !   0972  1.00000
     .  3668943d0,-152064160d0, 254858534d0, 187866593d0,-317406753d0,   !   0973  1.00000
     .  3779228d0,-152097524d0, 254101390d0, 187985223d0,-316696352d0,   !   0974  1.00000
     .  3889610d0,-152130369d0, 253342945d0, 188102987d0,-315984269d0,   !   0975  1.00000
     .  4000089d0,-152162693d0, 252583206d0, 188219885d0,-315270508d0,   !   0976  1.00000
     .  4110663d0,-152194499d0, 251822177d0, 188335918d0,-314555079d0,   !   0977  1.00000
     .  4221332d0,-152225786d0, 251059865d0, 188451088d0,-313837987d0,   !   0978  1.00000
     .  4332095d0,-152256555d0, 250296275d0, 188565394d0,-313119238d0,   !   0979  1.00000
     .  4442953d0,-152286807d0, 249531414d0, 188678838d0,-312398841d0,   !   0980  1.00000
     .  4553904d0,-152316541d0, 248765288d0, 188791421d0,-311676802d0,   !   0981  1.00000
     .  4664948d0,-152345759d0, 247997901d0, 188903144d0,-310953128d0,   !   0982  1.00000
     .  4776083d0,-152374461d0, 247229261d0, 189014006d0,-310227825d0,   !   0983  1.00000
     .  4887311d0,-152402647d0, 246459372d0, 189124010d0,-309500900d0,   !   0984  1.00000
     .  4998629d0,-152430318d0, 245688241d0, 189233156d0,-308772360d0,   !   0985  1.00000
     .  5110038d0,-152457474d0, 244915873d0, 189341445d0,-308042212d0,   !   0986  1.00000
     .  5221537d0,-152484117d0, 244142275d0, 189448877d0,-307310463d0,   !   0987  1.00000
     .  5333125d0,-152510246d0, 243367452d0, 189555454d0,-306577119d0,   !   0988  1.00000
     .  5444802d0,-152535862d0, 242591409d0, 189661176d0,-305842187d0,   !   0989  1.00000
     .  5556566d0,-152560965d0, 241814153d0, 189766044d0,-305105674d0,   !   0990  1.00000
     .  5668419d0,-152585556d0, 241035690d0, 189870060d0,-304367587d0,   !   0991  1.00000
     .  5780358d0,-152609636d0, 240256025d0, 189973223d0,-303627932d0,   !   0992  1.00000
     .  5892384d0,-152633204d0, 239475164d0, 190075534d0,-302886716d0,   !   0993  1.00000
     .  6004495d0,-152656262d0, 238693113d0, 190176996d0,-302143946d0,   !   0994  1.00000
     .  6116692d0,-152678810d0, 237909877d0, 190277607d0,-301399629d0,   !   0995  1.00000
     .  6228973d0,-152700849d0, 237125463d0, 190377370d0,-300653771d0,   !   0996  1.00000
     .  6341338d0,-152722378d0, 236339875d0, 190476285d0,-299906379d0,   !   0997  1.00000
     .  6453787d0,-152743399d0, 235553121d0, 190574353d0,-299157460d0,   !   0998  1.00000
     .  6566319d0,-152763912d0, 234765205d0, 190671574d0,-298407020d0,   !   0999  1.00000
     .  6678933d0,-152783917d0, 233976133d0, 190767950d0,-297655066d0,   !   1000  1.00000
     .  6791629d0,-152803416d0, 233185912d0, 190863481d0,-296901605d0,   !   1001  1.00000
     .  6904406d0,-152822407d0, 232394546d0, 190958169d0,-296146644d0,   !   1002  1.00000
     .  7017264d0,-152840893d0, 231602042d0, 191052013d0,-295390189d0,   !   1003  1.00000
     .  7130202d0,-152858873d0, 230808404d0, 191145015d0,-294632246d0,   !   1004  1.00000
     .  7243219d0,-152876348d0, 230013640d0, 191237176d0,-293872823d0,   !   1005  1.00000
     .  7356315d0,-152893319d0, 229217755d0, 191328497d0,-293111926d0,   !   1006  1.00000
     .  7469490d0,-152909785d0, 228420753d0, 191418978d0,-292349562d0,   !   1007  1.00000
     .  7582742d0,-152925748d0, 227622642d0, 191508620d0,-291585738d0,   !   1008  1.00000
     .  7696071d0,-152941208d0, 226823427d0, 191597424d0,-290820459d0,   !   1009  1.00000
     .  7809477d0,-152956166d0, 226023112d0, 191685391d0,-290053733d0,   !   1010  1.00000
     .  7922960d0,-152970621d0, 225221705d0, 191772522d0,-289285566d0,   !   1011  1.00000
     .  8036517d0,-152984575d0, 224419211d0, 191858817d0,-288515965d0,   !   1012  1.00000
     .  8150150d0,-152998028d0, 223615635d0, 191944278d0,-287744937d0,   !   1013  1.00000
     .  8263857d0,-153010980d0, 222810983d0, 192028905d0,-286972488d0,   !   1014  1.00000
     .  8377638d0,-153023432d0, 222005261d0, 192112699d0,-286198624d0,   !   1015  1.00000
     .  8491492d0,-153035385d0, 221198474d0, 192195660d0,-285423353d0,   !   1016  1.00000
     .  8605418d0,-153046839d0, 220390628d0, 192277791d0,-284646681d0,   !   1017  1.00000
     .  8719417d0,-153057794d0, 219581729d0, 192359091d0,-283868614d0,   !   1018  1.00000
     .  8833487d0,-153068251d0, 218771782d0, 192439561d0,-283089159d0,   !   1019  1.00000
     .  8947629d0,-153078211d0, 217960793d0, 192519203d0,-282308323d0,   !   1020  1.00000
     .  9061840d0,-153087674d0, 217148767d0, 192598017d0,-281526111d0,   !   1021  1.00000
     .  9176122d0,-153096640d0, 216335711d0, 192676003d0,-280742532d0,   !   1022  1.00000
     .  9290473d0,-153105110d0, 215521629d0, 192753164d0,-279957590d0,   !   1023  1.00000
     .  9404892d0,-153113085d0, 214706527d0, 192829499d0,-279171293d0,   !   1024  1.00000
     .  9519380d0,-153120565d0, 213890412d0, 192905009d0,-278383648d0,   !   1025  1.00000
     .  9633936d0,-153127550d0, 213073288d0, 192979696d0,-277594660d0,   !   1026  1.00000
     .  9748558d0,-153134041d0, 212255161d0, 193053559d0,-276804337d0,   !   1027  1.00000
     .  9863247d0,-153140039d0, 211436037d0, 193126601d0,-276012684d0,   !   1028  1.00000
     .  9978002d0,-153145543d0, 210615921d0, 193198821d0,-275219709d0,   !   1029  1.00000
     . 10092822d0,-153150556d0, 209794819d0, 193270221d0,-274425418d0,   !   1030  1.00000
     . 10207707d0,-153155076d0, 208972736d0, 193340802d0,-273629816d0,   !   1031  1.00000
     . 10322656d0,-153159104d0, 208149679d0, 193410563d0,-272832912d0,   !   1032  1.00000
     . 10437669d0,-153162642d0, 207325651d0, 193479507d0,-272034711d0,   !   1033  1.00000
     . 10552746d0,-153165689d0, 206500660d0, 193547634d0,-271235219d0,   !   1034  1.00000
     . 10667885d0,-153168247d0, 205674710d0, 193614944d0,-270434444d0,   !   1035  1.00000
     . 10783085d0,-153170314d0, 204847808d0, 193681439d0,-269632391d0,   !   1036  1.00000
     . 10898348d0,-153171893d0, 204019958d0, 193747119d0,-268829067d0,   !   1037  1.00000
     . 11013671d0,-153172983d0, 203191166d0, 193811986d0,-268024479d0,   !   1038  1.00000
     . 11129055d0,-153173586d0, 202361438d0, 193876040d0,-267218633d0,   !   1039  1.00000
     . 11244498d0,-153173700d0, 201530780d0, 193939281d0,-266411535d0,   !   1040  1.00000
     . 11360001d0,-153173328d0, 200699195d0, 194001712d0,-265603191d0,   !   1041  1.00000
     . 11475563d0,-153172469d0, 199866692d0, 194063332d0,-264793609d0,   !   1042  1.00000
     . 11591183d0,-153171125d0, 199033273d0, 194124142d0,-263982795d0,   !   1043  1.00000
     . 11706860d0,-153169295d0, 198198946d0, 194184144d0,-263170755d0,   !   1044  1.00000
     . 11822594d0,-153166980d0, 197363716d0, 194243338d0,-262357495d0,   !   1045  1.00000
     . 11938386d0,-153164180d0, 196527588d0, 194301725d0,-261543021d0,   !   1046  1.00000
     . 12054233d0,-153160896d0, 195690567d0, 194359306d0,-260727341d0,   !   1047  1.00000
     . 12170135d0,-153157129d0, 194852660d0, 194416081d0,-259910461d0,   !   1048  1.00000
     . 12286093d0,-153152879d0, 194013871d0, 194472052d0,-259092386d0,   !   1049  1.00000
     . 12402104d0,-153148147d0, 193174206d0, 194527219d0,-258273124d0,   !   1050  1.00000
     . 12518170d0,-153142932d0, 192333671d0, 194581583d0,-257452680d0,   !   1051  1.00000
     . 12634289d0,-153137237d0, 191492270d0, 194635145d0,-256631061d0,   !   1052  1.00000
     . 12750460d0,-153131060d0, 190650010d0, 194687906d0,-255808274d0,   !   1053  1.00000
     . 12866684d0,-153124402d0, 189806896d0, 194739866d0,-254984324d0,   !   1054  1.00000
     . 12982959d0,-153117265d0, 188962932d0, 194791027d0,-254159218d0,   !   1055  1.00000
     . 13099286d0,-153109648d0, 188118126d0, 194841389d0,-253332963d0,   !   1056  1.00000
     . 13215663d0,-153101553d0, 187272481d0, 194890954d0,-252505564d0,   !   1057  1.00000
     . 13332090d0,-153092979d0, 186426004d0, 194939721d0,-251677028d0,   !   1058  1.00000
     . 13448567d0,-153083926d0, 185578700d0, 194987693d0,-250847361d0,   !   1059  1.00000
     . 13565092d0,-153074397d0, 184730574d0, 195034868d0,-250016570d0,   !   1060  1.00000
     . 13681666d0,-153064390d0, 183881631d0, 195081250d0,-249184660d0,   !   1061  1.00000
     . 13798288d0,-153053907d0, 183031878d0, 195126837d0,-248351639d0,   !   1062  1.00000
     . 13914957d0,-153042949d0, 182181319d0, 195171632d0,-247517512d0,   !   1063  1.00000
     . 14031672d0,-153031514d0, 181329960d0, 195215635d0,-246682285d0,   !   1064  1.00000
     . 14148435d0,-153019605d0, 180477806d0, 195258847d0,-245845966d0,   !   1065  1.00000
     . 14265242d0,-153007221d0, 179624862d0, 195301268d0,-245008559d0,   !   1066  1.00000
     . 14382095d0,-152994363d0, 178771135d0, 195342900d0,-244170072d0,   !   1067  1.00000
     . 14498993d0,-152981032d0, 177916628d0, 195383743d0,-243330511d0,   !   1068  1.00000
     . 14615934d0,-152967228d0, 177061348d0, 195423798d0,-242489882d0,   !   1069  1.00000
     . 14732920d0,-152952952d0, 176205300d0, 195463066d0,-241648190d0,   !   1070  1.00000
     . 14849948d0,-152938203d0, 175348490d0, 195501549d0,-240805443d0,   !   1071  1.00000
     . 14967019d0,-152922984d0, 174490922d0, 195539245d0,-239961647d0,   !   1072  1.00000
     . 15084131d0,-152907293d0, 173632602d0, 195576158d0,-239116807d0,   !   1073  1.00000
     . 15201285d0,-152891132d0, 172773535d0, 195612286d0,-238270930d0,   !   1074  1.00000
     . 15318480d0,-152874500d0, 171913726d0, 195647632d0,-237424022d0,   !   1075  1.00000
     . 15435716d0,-152857400d0, 171053182d0, 195682196d0,-236576090d0,   !   1076  1.00000
     . 15552991d0,-152839830d0, 170191906d0, 195715979d0,-235727139d0,   !   1077  1.00000
     . 15670305d0,-152821792d0, 169329905d0, 195748981d0,-234877176d0,   !   1078  1.00000
     . 15787659d0,-152803286d0, 168467184d0, 195781204d0,-234026207d0,   !   1079  1.00000
     . 15905050d0,-152784313d0, 167603748d0, 195812648d0,-233174237d0,   !   1080  1.00000
     . 16022480d0,-152764873d0, 166739602d0, 195843314d0,-232321274d0,   !   1081  1.00000
     . 16139946d0,-152744966d0, 165874752d0, 195873204d0,-231467324d0,   !   1082  1.00000
     . 16257449d0,-152724594d0, 165009203d0, 195902317d0,-230612391d0,   !   1083  1.00000
     . 16374988d0,-152703756d0, 164142959d0, 195930655d0,-229756484d0,   !   1084  1.00000
     . 16492563d0,-152682453d0, 163276027d0, 195958219d0,-228899607d0,   !   1085  1.00000
     . 16610173d0,-152660686d0, 162408412d0, 195985009d0,-228041767d0,   !   1086  1.00000
     . 16727818d0,-152638455d0, 161540119d0, 196011026d0,-227182970d0,   !   1087  1.00000
     . 16845496d0,-152615760d0, 160671152d0, 196036271d0,-226323222d0,   !   1088  1.00000
     . 16963209d0,-152592602d0, 159801518d0, 196060745d0,-225462529d0,   !   1089  1.00000
     . 17080954d0,-152568983d0, 158931222d0, 196084449d0,-224600898d0,   !   1090  1.00000
     . 17198731d0,-152544901d0, 158060268d0, 196107383d0,-223738333d0,   !   1091  1.00000
     . 17316541d0,-152520357d0, 157188662d0, 196129549d0,-222874843d0,   !   1092  1.00000
     . 17434382d0,-152495353d0, 156316410d0, 196150946d0,-222010432d0,   !   1093  1.00000
     . 17552253d0,-152469889d0, 155443516d0, 196171577d0,-221145106d0,   !   1094  1.00000
     . 17670156d0,-152443964d0, 154569985d0, 196191442d0,-220278872d0,   !   1095  1.00000
     . 17788088d0,-152417581d0, 153695824d0, 196210542d0,-219411736d0,   !   1096  1.00000
     . 17906049d0,-152390738d0, 152821036d0, 196228877d0,-218543704d0,   !   1097  1.00000
     . 18024039d0,-152363437d0, 151945628d0, 196246448d0,-217674782d0,   !   1098  1.00000
     . 18142058d0,-152335678d0, 151069604d0, 196263257d0,-216804975d0,   !   1099  1.00000
     . 18260104d0,-152307461d0, 150192969d0, 196279304d0,-215934290d0,   !   1100  1.00000
     . 18378178d0,-152278788d0, 149315729d0, 196294589d0,-215062734d0,   !   1101  1.00000
     . 18496278d0,-152249658d0, 148437889d0, 196309115d0,-214190311d0,   !   1102  1.00000
     . 18614405d0,-152220073d0, 147559454d0, 196322881d0,-213317028d0,   !   1103  1.00000
     . 18732557d0,-152190032d0, 146680430d0, 196335888d0,-212442891d0,   !   1104  1.00000
     . 18850735d0,-152159536d0, 145800820d0, 196348138d0,-211567906d0,   !   1105  1.00000
     . 18968937d0,-152128586d0, 144920631d0, 196359630d0,-210692079d0,   !   1106  1.00000
     . 19087164d0,-152097182d0, 144039868d0, 196370367d0,-209815416d0,   !   1107  1.00000
     . 19205414d0,-152065324d0, 143158535d0, 196380348d0,-208937922d0,   !   1108  1.00000
     . 19323687d0,-152033014d0, 142276638d0, 196389575d0,-208059605d0,   !   1109  1.00000
     . 19441983d0,-152000251d0, 141394182d0, 196398049d0,-207180470d0,   !   1110  1.00000
     . 19560301d0,-151967037d0, 140511172d0, 196405769d0,-206300522d0,   !   1111  1.00000
     . 19678641d0,-151933371d0, 139627614d0, 196412738d0,-205419768d0,   !   1112  1.00000
     . 19797002d0,-151899254d0, 138743511d0, 196418956d0,-204538214d0,   !   1113  1.00000
     . 19915383d0,-151864687d0, 137858870d0, 196424423d0,-203655865d0,   !   1114  1.00000
     . 20033785d0,-151829669d0, 136973695d0, 196429141d0,-202772728d0,   !   1115  1.00000
     . 20152206d0,-151794203d0, 136087992d0, 196433111d0,-201888808d0,   !   1116  1.00000
     . 20270646d0,-151758288d0, 135201766d0, 196436332d0,-201004112d0,   !   1117  1.00000
     . 20389105d0,-151721924d0, 134315021d0, 196438807d0,-200118645d0,   !   1118  1.00000
     . 20507582d0,-151685113d0, 133427763d0, 196440536d0,-199232414d0,   !   1119  1.00000
     . 20626076d0,-151647854d0, 132539996d0, 196441519d0,-198345423d0,   !   1120  1.00000
     . 20744588d0,-151610148d0, 131651727d0, 196441758d0,-197457679d0,   !   1121  1.00000
     . 20863116d0,-151571996d0, 130762959d0, 196441254d0,-196569189d0,   !   1122  1.00000
     . 20981660d0,-151533398d0, 129873698d0, 196440007d0,-195679957d0,   !   1123  1.00000
     . 21100219d0,-151494354d0, 128983949d0, 196438018d0,-194789990d0,   !   1124  1.00000
     . 21218794d0,-151454866d0, 128093717d0, 196435288d0,-193899293d0,   !   1125  1.00000
     . 21337383d0,-151414934d0, 127203008d0, 196431817d0,-193007873d0,   !   1126  1.00000
     . 21455986d0,-151374557d0, 126311825d0, 196427608d0,-192115735d0,   !   1127  1.00000
     . 21574603d0,-151333738d0, 125420174d0, 196422659d0,-191222885d0,   !   1128  1.00000
     . 21693233d0,-151292475d0, 124528060d0, 196416973d0,-190329329d0,   !   1129  1.00000
     . 21811875d0,-151250770d0, 123635489d0, 196410551d0,-189435072d0,   !   1130  1.00000
     . 21930530d0,-151208623d0, 122742464d0, 196403392d0,-188540121d0,   !   1131  1.00000
     . 22049196d0,-151166035d0, 121848992d0, 196395498d0,-187644482d0,   !   1132  1.00000
     . 22167873d0,-151123006d0, 120955076d0, 196386869d0,-186748160d0,   !   1133  1.00000
     . 22286560d0,-151079537d0, 120060723d0, 196377507d0,-185851161d0,   !   1134  1.00000
     . 22405258d0,-151035628d0, 119165936d0, 196367413d0,-184953490d0,   !   1135  1.00000
     . 22523965d0,-150991279d0, 118270721d0, 196356586d0,-184055155d0,   !   1136  1.00000
     . 22642681d0,-150946492d0, 117375084d0, 196345029d0,-183156159d0,   !   1137  1.00000
     . 22761405d0,-150901266d0, 116479028d0, 196332741d0,-182256510d0,   !   1138  1.00000
     . 22880138d0,-150855603d0, 115582559d0, 196319724d0,-181356213d0,   !   1139  1.00000
     . 22998878d0,-150809502d0, 114685682d0, 196305978d0,-180455273d0,   !   1140  1.00000
     . 23117626d0,-150762965d0, 113788401d0, 196291505d0,-179553697d0,   !   1141  1.00000
     . 23236380d0,-150715991d0, 112890722d0, 196276305d0,-178651490d0,   !   1142  1.00000
     . 23355140d0,-150668581d0, 111992650d0, 196260379d0,-177748658d0,   !   1143  1.00000
     . 23473905d0,-150620736d0, 111094189d0, 196243727d0,-176845207d0,   !   1144  1.00000
     . 23592676d0,-150572456d0, 110195344d0, 196226352d0,-175941142d0,   !   1145  1.00000
     . 23711451d0,-150523742d0, 109296121d0, 196208253d0,-175036469d0,   !   1146  1.00000
     . 23830231d0,-150474593d0, 108396524d0, 196189431d0,-174131194d0,   !   1147  1.00000
     . 23949014d0,-150425012d0, 107496559d0, 196169887d0,-173225323d0,   !   1148  1.00000
     . 24067800d0,-150374998d0, 106596229d0, 196149622d0,-172318861d0,   !   1149  1.00000
     . 24186588d0,-150324551d0, 105695540d0, 196128638d0,-171411814d0,   !   1150  1.00000
     . 24305379d0,-150273673d0, 104794497d0, 196106933d0,-170504188d0,   !   1151  1.00000
     . 24424172d0,-150222363d0, 103893104d0, 196084511d0,-169595988d0,   !   1152  1.00000
     . 24542966d0,-150170623d0, 102991367d0, 196061370d0,-168687221d0,   !   1153  1.00000
     . 24661760d0,-150118452d0, 102089291d0, 196037513d0,-167777891d0,   !   1154  1.00000
     . 24780555d0,-150065852d0, 101186879d0, 196012940d0,-166868004d0,   !   1155  1.00000
     . 24899349d0,-150012822d0, 100284138d0, 195987652d0,-165957567d0,   !   1156  1.00000
     . 25018143d0,-149959363d0,  99381072d0, 195961649d0,-165046584d0,   !   1157  1.00000
     . 25136935d0,-149905476d0,  98477686d0, 195934933d0,-164135061d0,   !   1158  1.00000
     . 25255726d0,-149851161d0,  97573984d0, 195907505d0,-163223005d0,   !   1159  1.00000
     . 25374514d0,-149796419d0,  96669971d0, 195879364d0,-162310420d0,   !   1160  1.00000
     . 25493300d0,-149741250d0,  95765653d0, 195850513d0,-161397313d0,   !   1161  1.00000
     . 25612083d0,-149685655d0,  94861035d0, 195820951d0,-160483688d0,   !   1162  1.00000
     . 25730862d0,-149629634d0,  93956120d0, 195790680d0,-159569552d0,   !   1163  1.00000
     . 25849636d0,-149573188d0,  93050913d0, 195759701d0,-158654910d0,   !   1164  1.00000
     . 25968407d0,-149516317d0,  92145421d0, 195728014d0,-157739768d0,   !   1165  1.00000
     . 26087172d0,-149459022d0,  91239646d0, 195695620d0,-156824131d0,   !   1166  1.00000
     . 26205931d0,-149401302d0,  90333595d0, 195662520d0,-155908004d0,   !   1167  1.00000
     . 26324684d0,-149343160d0,  89427272d0, 195628715d0,-154991395d0,   !   1168  1.00000
     . 26443431d0,-149284595d0,  88520681d0, 195594206d0,-154074307d0,   !   1169  1.00000
     . 26562171d0,-149225607d0,  87613828d0, 195558993d0,-153156747d0,   !   1170  1.00000
     . 26680903d0,-149166198d0,  86706717d0, 195523078d0,-152238721d0,   !   1171  1.00000
     . 26799628d0,-149106367d0,  85799353d0, 195486461d0,-151320233d0,   !   1172  1.00000
     . 26918344d0,-149046116d0,  84891741d0, 195449143d0,-150401289d0,   !   1173  1.00000
     . 27037051d0,-148985444d0,  83983886d0, 195411125d0,-149481895d0,   !   1174  1.00000
     . 27155748d0,-148924353d0,  83075791d0, 195372408d0,-148562056d0,   !   1175  1.00000
     . 27274436d0,-148862842d0,  82167463d0, 195332992d0,-147641779d0,   !   1176  1.00000
     . 27393113d0,-148800913d0,  81258905d0, 195292878d0,-146721067d0,   !   1177  1.00000
     . 27511780d0,-148738565d0,  80350123d0, 195252068d0,-145799928d0,   !   1178  1.00000
     . 27630435d0,-148675800d0,  79441121d0, 195210562d0,-144878366d0,   !   1179  1.00000
     . 27749079d0,-148612617d0,  78531904d0, 195168361d0,-143956387d0,   !   1180  1.00000
     . 27867710d0,-148549017d0,  77622477d0, 195125465d0,-143033997d0,   !   1181  1.00000
     . 27986328d0,-148485002d0,  76712844d0, 195081877d0,-142111200d0,   !   1182  1.00000
     . 28104934d0,-148420570d0,  75803010d0, 195037595d0,-141188003d0,   !   1183  1.00000
     . 28223525d0,-148355723d0,  74892980d0, 194992622d0,-140264410d0,   !   1184  1.00000
     . 28342103d0,-148290462d0,  73982758d0, 194946958d0,-139340428d0,   !   1185  1.00000
     . 28460666d0,-148224786d0,  73072349d0, 194900603d0,-138416062d0,   !   1186  1.00000
     . 28579214d0,-148158696d0,  72161758d0, 194853560d0,-137491317d0,   !   1187  1.00000
     . 28697746d0,-148092193d0,  71250990d0, 194805828d0,-136566199d0,   !   1188  1.00000
     . 28816263d0,-148025278d0,  70340048d0, 194757408d0,-135640713d0,   !   1189  1.00000
     . 28934763d0,-147957950d0,  69428939d0, 194708302d0,-134714865d0,   !   1190  1.00000
     . 29053246d0,-147890210d0,  68517666d0, 194658510d0,-133788659d0,   !   1191  1.00000
     . 29171712d0,-147822060d0,  67606234d0, 194608033d0,-132862102d0,   !   1192  1.00000
     . 29290159d0,-147753498d0,  66694647d0, 194556871d0,-131935199d0,   !   1193  1.00000
     . 29408589d0,-147684526d0,  65782911d0, 194505026d0,-131007955d0,   !   1194  1.00000
     . 29527000d0,-147615144d0,  64871031d0, 194452499d0,-130080376d0,   !   1195  1.00000
     . 29645391d0,-147545354d0,  63959009d0, 194399289d0,-129152467d0,   !   1196  1.00000
     . 29763763d0,-147475154d0,  63046852d0, 194345399d0,-128224233d0,   !   1197  1.00000
     . 29882115d0,-147404546d0,  62134564d0, 194290829d0,-127295681d0,   !   1198  1.00000
     . 30000446d0,-147333531d0,  61222149d0, 194235580d0,-126366814d0,   !   1199  1.00000
     . 30118755d0,-147262108d0,  60309613d0, 194179652d0,-125437639d0,   !   1200  1.00000
     . 30237044d0,-147190278d0,  59396959d0, 194123047d0,-124508161d0,   !   1201  1.00000
     . 30355310d0,-147118042d0,  58484192d0, 194065765d0,-123578385d0,   !   1202  1.00000
     . 30473554d0,-147045401d0,  57571318d0, 194007807d0,-122648317d0,   !   1203  1.00000
     . 30591775d0,-146972354d0,  56658339d0, 193949174d0,-121717962d0,   !   1204  1.00000
     . 30709973d0,-146898903d0,  55745262d0, 193889867d0,-120787325d0,   !   1205  1.00000
     . 30828146d0,-146825047d0,  54832090d0, 193829886d0,-119856412d0,   !   1206  1.00000
     . 30946296d0,-146750787d0,  53918829d0, 193769233d0,-118925228d0,   !   1207  1.00000
     . 31064421d0,-146676125d0,  53005482d0, 193707908d0,-117993778d0,   !   1208  1.00000
     . 31182520d0,-146601059d0,  52092055d0, 193645912d0,-117062068d0,   !   1209  1.00000
     . 31300594d0,-146525592d0,  51178551d0, 193583247d0,-116130103d0,   !   1210  1.00000
     . 31418642d0,-146449722d0,  50264976d0, 193519912d0,-115197888d0,   !   1211  1.00000
     . 31536663d0,-146373452d0,  49351333d0, 193455908d0,-114265429d0,   !   1212  1.00000
     . 31654657d0,-146296781d0,  48437628d0, 193391238d0,-113332730d0,   !   1213  1.00000
     . 31772623d0,-146219709d0,  47523865d0, 193325900d0,-112399797d0,   !   1214  1.00000
     . 31890562d0,-146142238d0,  46610049d0, 193259897d0,-111466636d0,   !   1215  1.00000
     . 32008472d0,-146064367d0,  45696183d0, 193193228d0,-110533251d0,   !   1216  1.00000
     . 32126354d0,-145986098d0,  44782274d0, 193125896d0,-109599648d0,   !   1217  1.00000
     . 32244206d0,-145907431d0,  43868324d0, 193057900d0,-108665833d0,   !   1218  1.00000
     . 32362029d0,-145828366d0,  42954338d0, 192989241d0,-107731809d0,   !   1219  1.00000
     . 32479821d0,-145748904d0,  42040322d0, 192919921d0,-106797584d0,   !   1220  1.00000
     . 32597583d0,-145669045d0,  41126279d0, 192849940d0,-105863161d0,   !   1221  1.00000
     . 32715313d0,-145588789d0,  40212214d0, 192779299d0,-104928546d0,   !   1222  1.00000
     . 32833012d0,-145508138d0,  39298131d0, 192707999d0,-103993744d0,   !   1223  1.00000
     . 32950679d0,-145427092d0,  38384035d0, 192636041d0,-103058761d0,   !   1224  1.00000
     . 33068314d0,-145345652d0,  37469931d0, 192563425d0,-102123601d0,   !   1225  1.00000
     . 33185916d0,-145263817d0,  36555823d0, 192490152d0,-101188271d0,   !   1226  1.00000
     . 33303484d0,-145181588d0,  35641715d0, 192416223d0,-100252774d0,   !   1227  1.00000
     . 33421019d0,-145098966d0,  34727611d0, 192341640d0, -99317117d0,   !   1228  1.00000
     . 33538519d0,-145015951d0,  33813517d0, 192266402d0, -98381304d0,   !   1229  1.00000
     . 33655985d0,-144932545d0,  32899436d0, 192190511d0, -97445341d0,   !   1230  1.00000
     . 33773416d0,-144848746d0,  31985374d0, 192113968d0, -96509233d0,   !   1231  1.00000
     . 33890811d0,-144764557d0,  31071333d0, 192036773d0, -95572984d0,   !   1232  1.00000
     . 34008170d0,-144679976d0,  30157320d0, 191958927d0, -94636601d0,   !   1233  1.00000
     . 34125493d0,-144595006d0,  29243338d0, 191880431d0, -93700088d0,   !   1234  1.00000
     . 34242779d0,-144509645d0,  28329392d0, 191801286d0, -92763450d0,   !   1235  1.00000
     . 34360028d0,-144423896d0,  27415486d0, 191721492d0, -91826693d0,   !   1236  1.00000
     . 34477239d0,-144337758d0,  26501624d0, 191641051d0, -90889821d0,   !   1237  1.00000
     . 34594412d0,-144251232d0,  25587811d0, 191559964d0, -89952840d0,   !   1238  1.00000
     . 34711546d0,-144164318d0,  24674051d0, 191478230d0, -89015755d0,   !   1239  1.00000
     . 34828641d0,-144077017d0,  23760350d0, 191395852d0, -88078571d0,   !   1240  1.00000
     . 34945697d0,-143989329d0,  22846710d0, 191312830d0, -87141293d0,   !   1241  1.00000
     . 35062712d0,-143901255d0,  21933136d0, 191229164d0, -86203926d0,   !   1242  1.00000
     . 35179688d0,-143812795d0,  21019633d0, 191144856d0, -85266475d0,   !   1243  1.00000
     . 35296622d0,-143723950d0,  20106206d0, 191059906d0, -84328946d0,   !   1244  1.00000
     . 35413516d0,-143634721d0,  19192857d0, 190974315d0, -83391343d0,   !   1245  1.00000
     . 35530367d0,-143545107d0,  18279593d0, 190888085d0, -82453671d0,   !   1246  1.00000
     . 35647177d0,-143455110d0,  17366417d0, 190801215d0, -81515936d0,   !   1247  1.00000
     . 35763944d0,-143364729d0,  16453333d0, 190713707d0, -80578143d0,   !   1248  1.00000
     . 35880668d0,-143273966d0,  15540346d0, 190625561d0, -79640296d0,   !   1249  1.00000
     . 35997349d0,-143182821d0,  14627460d0, 190536779d0, -78702401d0,   !   1250  1.00000
     . 36113985d0,-143091294d0,  13714679d0, 190447361d0, -77764462d0,   !   1251  1.00000
     . 36230578d0,-142999386d0,  12802008d0, 190357309d0, -76826485d0,   !   1252  1.00000
     . 36347126d0,-142907097d0,  11889452d0, 190266622d0, -75888475d0,   !   1253  1.00000
     . 36463628d0,-142814428d0,  10977013d0, 190175301d0, -74950437d0,   !   1254  1.00000
     . 36580086d0,-142721380d0,  10064697d0, 190083349d0, -74012376d0,   !   1255  1.00000
     . 36696497d0,-142627952d0,   9152509d0, 189990765d0, -73074296d0,   !   1256  1.00000
     . 36812862d0,-142534145d0,   8240451d0, 189897550d0, -72136203d0,   !   1257  1.00000
     . 36929180d0,-142439961d0,   7328529d0, 189803705d0, -71198101d0,   !   1258  1.00000
     . 37045450d0,-142345399d0,   6416746d0, 189709231d0, -70259997d0,   !   1259  1.00000
     . 37161673d0,-142250459d0,   5505108d0, 189614128d0, -69321894d0,   !   1260  1.00000
     . 37277848d0,-142155143d0,   4593618d0, 189518399d0, -68383797d0,   !   1261  1.00000
     . 37393974d0,-142059451d0,   3682280d0, 189422042d0, -67445713d0,   !   1262  1.00000
     . 37510052d0,-141963383d0,   2771100d0, 189325060d0, -66507644d0,   !   1263  1.00000
     . 37626080d0,-141866940d0,   1860080d0, 189227453d0, -65569598d0,   !   1264  1.00000
     . 37742058d0,-141770123d0,    949225d0, 189129222d0, -64631578d0,   !   1265  1.00000
     . 37857986d0,-141672931d0,     38540d0, 189030368d0, -63693589d0,   !   1266  1.00000
     . 37973863d0,-141575366d0,   -871970d0, 188930891d0, -62755636d0,   !   1267  1.00000
     . 38089689d0,-141477427d0,  -1782303d0, 188830793d0, -61817725d0,   !   1268  1.00000
     . 38205464d0,-141379116d0,  -2692455d0, 188730073d0, -60879860d0,   !   1269  1.00000
     . 38321186d0,-141280433d0,  -3602420d0, 188628734d0, -59942046d0,   !   1270  1.00000
     . 38436857d0,-141181378d0,  -4512194d0, 188526776d0, -59004288d0,   !   1271  1.00000
     . 38552474d0,-141081952d0,  -5421773d0, 188424200d0, -58066590d0,   !   1272  1.00000
     . 38668039d0,-140982156d0,  -6331154d0, 188321006d0, -57128959d0,   !   1273  1.00000
     . 38783549d0,-140881989d0,  -7240332d0, 188217196d0, -56191398d0,   !   1274  1.00000
     . 38899006d0,-140781452d0,  -8149302d0, 188112770d0, -55253912d0,   !   1275  1.00000
     . 39014408d0,-140680547d0,  -9058060d0, 188007729d0, -54316507d0,   !   1276  1.00000
     . 39129756d0,-140579273d0,  -9966603d0, 187902074d0, -53379187d0,   !   1277  1.00000
     . 39245048d0,-140477631d0, -10874927d0, 187795806d0, -52441958d0,   !   1278  1.00000
     . 39360285d0,-140375621d0, -11783026d0, 187688926d0, -51504823d0,   !   1279  1.00000
     . 39475465d0,-140273244d0, -12690897d0, 187581433d0, -50567788d0,   !   1280  1.00000
     . 39590589d0,-140170500d0, -13598536d0, 187473331d0, -49630857d0,   !   1281  1.00000
     . 39705656d0,-140067390d0, -14505938d0, 187364618d0, -48694036d0,   !   1282  1.00000
     . 39820665d0,-139963915d0, -15413100d0, 187255296d0, -47757329d0,   !   1283  1.00000
     . 39935617d0,-139860075d0, -16320017d0, 187145366d0, -46820741d0,   !   1284  1.00000
     . 40050511d0,-139755870d0, -17226686d0, 187034829d0, -45884277d0,   !   1285  1.00000
     . 40165346d0,-139651300d0, -18133101d0, 186923685d0, -44947942d0,   !   1286  1.00000
     . 40280122d0,-139546368d0, -19039260d0, 186811935d0, -44011740d0,   !   1287  1.00000
     . 40394838d0,-139441072d0, -19945157d0, 186699581d0, -43075676d0,   !   1288  1.00000
     . 40509495d0,-139335413d0, -20850789d0, 186586622d0, -42139756d0,   !   1289  1.00000
     . 40624092d0,-139229393d0, -21756152d0, 186473060d0, -41203983d0,   !   1290  1.00000
     . 40738628d0,-139123011d0, -22661241d0, 186358896d0, -40268363d0,   !   1291  1.00000
     . 40853102d0,-139016267d0, -23566053d0, 186244131d0, -39332900d0,   !   1292  1.00000
     . 40967516d0,-138909164d0, -24470584d0, 186128765d0, -38397599d0,   !   1293  1.00000
     . 41081867d0,-138801700d0, -25374828d0, 186012798d0, -37462465d0,   !   1294  1.00000
     . 41196157d0,-138693876d0, -26278783d0, 185896233d0, -36527503d0,   !   1295  1.00000
     . 41310383d0,-138585693d0, -27182444d0, 185779070d0, -35592717d0,   !   1296  1.00000
     . 41424547d0,-138477152d0, -28085807d0, 185661309d0, -34658112d0,   !   1297  1.00000
     . 41538647d0,-138368253d0, -28988869d0, 185542952d0, -33723693d0,   !   1298  1.00000
     . 41652683d0,-138258996d0, -29891624d0, 185423999d0, -32789464d0,   !   1299  1.00000
     . 41766655d0,-138149382d0, -30794070d0, 185304451d0, -31855430d0,   !   1300  1.00000
     . 41880562d0,-138039411d0, -31696201d0, 185184309d0, -30921597d0,   !   1301  1.00000
     . 41994404d0,-137929085d0, -32598015d0, 185063574d0, -29987968d0,   !   1302  1.00000
     . 42108181d0,-137818403d0, -33499506d0, 184942247d0, -29054548d0,   !   1303  1.00000
     . 42221891d0,-137707365d0, -34400672d0, 184820328d0, -28121342d0,   !   1304  1.00000
     . 42335535d0,-137595973d0, -35301507d0, 184697819d0, -27188355d0,   !   1305  1.00000
     . 42449113d0,-137484227d0, -36202008d0, 184574720d0, -26255591d0,   !   1306  1.00000
     . 42562624d0,-137372128d0, -37102171d0, 184451031d0, -25323055d0,   !   1307  1.00000
     . 42676067d0,-137259675d0, -38001992d0, 184326755d0, -24390752d0,   !   1308  1.00000
     . 42789442d0,-137146870d0, -38901467d0, 184201891d0, -23458686d0,   !   1309  1.00000
     . 42902748d0,-137033713d0, -39800592d0, 184076441d0, -22526861d0,   !   1310  1.00000
     . 43015987d0,-136920204d0, -40699363d0, 183950405d0, -21595284d0,   !   1311  1.00000
     . 43129156d0,-136806344d0, -41597776d0, 183823784d0, -20663958d0,   !   1312  1.00000
     . 43242255d0,-136692133d0, -42495827d0, 183696580d0, -19732887d0,   !   1313  1.00000
     . 43355285d0,-136577573d0, -43393512d0, 183568792d0, -18802077d0,   !   1314  1.00000
     . 43468244d0,-136462663d0, -44290828d0, 183440422d0, -17871532d0,   !   1315  1.00000
     . 43581132d0,-136347404d0, -45187769d0, 183311470d0, -16941256d0,   !   1316  1.00000
     . 43693950d0,-136231796d0, -46084333d0, 183181938d0, -16011255d0,   !   1317  1.00000
     . 43806696d0,-136115840d0, -46980515d0, 183051826d0, -15081533d0,   !   1318  1.00000
     . 43919370d0,-135999537d0, -47876311d0, 182921134d0, -14152094d0,   !   1319  1.00000
     . 44031972d0,-135882886d0, -48771718d0, 182789865d0, -13222942d0,   !   1320  1.00000
     . 44144501d0,-135765889d0, -49666732d0, 182658019d0, -12294084d0,   !   1321  1.00000
     . 44256958d0,-135648546d0, -50561348d0, 182525596d0, -11365523d0,   !   1322  1.00000
     . 44369341d0,-135530858d0, -51455562d0, 182392597d0, -10437263d0,   !   1323  1.00000
     . 44481649d0,-135412824d0, -52349372d0, 182259024d0,  -9509309d0,   !   1324  1.00000
     . 44593884d0,-135294446d0, -53242772d0, 182124876d0,  -8581667d0,   !   1325  1.00000
     . 44706044d0,-135175724d0, -54135760d0, 181990156d0,  -7654339d0,   !   1326  1.00000
     . 44818130d0,-135056658d0, -55028330d0, 181854863d0,  -6727332d0,   !   1327  1.00000
     . 44930139d0,-134937249d0, -55920480d0, 181718999d0,  -5800648d0,   !   1328  1.00000
     . 45042073d0,-134817498d0, -56812206d0, 181582564d0,  -4874294d0,   !   1329  1.00000
     . 45153931d0,-134697405d0, -57703503d0, 181445560d0,  -3948273d0,   !   1330  1.00000
     . 45265713d0,-134576970d0, -58594367d0, 181307986d0,  -3022590d0,   !   1331  1.00000
     . 45377417d0,-134456194d0, -59484795d0, 181169844d0,  -2097249d0,   !   1332  1.00000
     . 45489044d0,-134335077d0, -60374784d0, 181031135d0,  -1172256d0,   !   1333  1.00000
     . 45600594d0,-134213621d0, -61264328d0, 180891860d0,   -247613d0,   !   1334  1.00000
     . 45712065d0,-134091825d0, -62153425d0, 180752019d0,    676673d0,   !   1335  1.00000
     . 45823458d0,-133969690d0, -63042070d0, 180611613d0,   1600599d0,   !   1336  1.00000
     . 45934772d0,-133847216d0, -63930259d0, 180470644d0,   2524160d0,   !   1337  1.00000
     . 46046007d0,-133724405d0, -64817990d0, 180329111d0,   3447353d0,   !   1338  1.00000
     . 46157162d0,-133601255d0, -65705257d0, 180187016d0,   4370171d0,   !   1339  1.00000
     . 46268237d0,-133477769d0, -66592058d0, 180044359d0,   5292612d0,   !   1340  1.00000
     . 46379232d0,-133353947d0, -67478387d0, 179901142d0,   6214670d0,   !   1341  1.00000
     . 46490146d0,-133229788d0, -68364243d0, 179757365d0,   7136341d0,   !   1342  1.00000
     . 46600979d0,-133105294d0, -69249620d0, 179613030d0,   8057620d0,   !   1343  1.00000
     . 46711730d0,-132980464d0, -70134515d0, 179468136d0,   8978504d0,   !   1344  1.00000
     . 46822400d0,-132855300d0, -71018925d0, 179322685d0,   9898987d0,   !   1345  1.00000
     . 46932987d0,-132729803d0, -71902844d0, 179176677d0,  10819066d0,   !   1346  1.00000
     . 47043491d0,-132603971d0, -72786271d0, 179030114d0,  11738736d0,   !   1347  1.00000
     . 47153913d0,-132477807d0, -73669200d0, 178882996d0,  12657992d0,   !   1348  1.00000
     . 47264251d0,-132351310d0, -74551628d0, 178735325d0,  13576831d0,   !   1349  1.00000
     . 47374506d0,-132224481d0, -75433552d0, 178587100d0,  14495247d0,   !   1350  1.00000
     . 47484676d0,-132097321d0, -76314967d0, 178438323d0,  15413236d0,   !   1351  1.00000
     . 47594762d0,-131969830d0, -77195870d0, 178288994d0,  16330795d0,   !   1352  1.00000
     . 47704762d0,-131842008d0, -78076258d0, 178139115d0,  17247919d0,   !   1353  1.00000
     . 47814678d0,-131713856d0, -78956125d0, 177988687d0,  18164602d0,   !   1354  1.00000
     . 47924508d0,-131585374d0, -79835469d0, 177837709d0,  19080842d0,   !   1355  1.00000
     . 48034252d0,-131456564d0, -80714287d0, 177686183d0,  19996633d0,   !   1356  1.00000
     . 48143909d0,-131327425d0, -81592573d0, 177534110d0,  20911972d0,   !   1357  1.00000
     . 48253480d0,-131197958d0, -82470325d0, 177381491d0,  21826854d0,   !   1358  1.00000
     . 48362964d0,-131068164d0, -83347539d0, 177228326d0,  22741274d0,   !   1359  1.00000
     . 48472360d0,-130938043d0, -84224211d0, 177074617d0,  23655229d0,   !   1360  1.00000
     . 48581668d0,-130807595d0, -85100337d0, 176920363d0,  24568714d0,   !   1361  1.00000
     . 48690888d0,-130676821d0, -85975914d0, 176765567d0,  25481725d0,   !   1362  1.00000
     . 48800019d0,-130545722d0, -86850938d0, 176610228d0,  26394257d0,   !   1363  1.00000
     . 48909062d0,-130414297d0, -87725406d0, 176454348d0,  27306306d0,   !   1364  1.00000
     . 49018015d0,-130282548d0, -88599313d0, 176297927d0,  28217868d0,   !   1365  1.00000
     . 49126878d0,-130150476d0, -89472656d0, 176140967d0,  29128939d0,   !   1366  1.00000
     . 49235652d0,-130018079d0, -90345432d0, 175983468d0,  30039514d0,   !   1367  1.00000
     . 49344334d0,-129885360d0, -91217636d0, 175825431d0,  30949589d0,   !   1368  1.00000
     . 49452926d0,-129752318d0, -92089265d0, 175666856d0,  31859160d0,   !   1369  1.00000
     . 49561427d0,-129618955d0, -92960316d0, 175507746d0,  32768223d0,   !   1370  1.00000
     . 49669837d0,-129485269d0, -93830785d0, 175348099d0,  33676773d0,   !   1371  1.00000
     . 49778154d0,-129351263d0, -94700668d0, 175187918d0,  34584807d0,   !   1372  1.00000
     . 49886379d0,-129216936d0, -95569961d0, 175027203d0,  35492319d0,   !   1373  1.00000
     . 49994512d0,-129082290d0, -96438662d0, 174865955d0,  36399306d0,   !   1374  1.00000
     . 50102552d0,-128947324d0, -97306765d0, 174704175d0,  37305763d0,   !   1375  1.00000
     . 50210498d0,-128812039d0, -98174269d0, 174541864d0,  38211687d0,   !   1376  1.00000
     . 50318350d0,-128676435d0, -99041168d0, 174379022d0,  39117072d0,   !   1377  1.00000
     . 50426109d0,-128540514d0, -99907460d0, 174215650d0,  40021916d0,   !   1378  1.00000
     . 50533773d0,-128404275d0,-100773141d0, 174051750d0,  40926213d0,   !   1379  1.00000
     . 50641342d0,-128267719d0,-101638207d0, 173887322d0,  41829960d0,   !   1380  1.00000
     . 50748816d0,-128130846d0,-102502655d0, 173722366d0,  42733152d0,   !   1381  1.00000
     . 50856195d0,-127993658d0,-103366481d0, 173556884d0,  43635785d0,   !   1382  1.00000
     . 50963478d0,-127856154d0,-104229682d0, 173390876d0,  44537856d0,   !   1383  1.00000
     . 51070664d0,-127718335d0,-105092254d0, 173224344d0,  45439359d0,   !   1384  1.00000
     . 51177754d0,-127580202d0,-105954193d0, 173057288d0,  46340291d0,   !   1385  1.00000
     . 51284747d0,-127441755d0,-106815497d0, 172889709d0,  47240647d0,   !   1386  1.00000
     . 51391643d0,-127302994d0,-107676160d0, 172721608d0,  48140424d0,   !   1387  1.00000
     . 51498441d0,-127163921d0,-108536181d0, 172552986d0,  49039617d0,   !   1388  1.00000
     . 51605141d0,-127024535d0,-109395555d0, 172383843d0,  49938222d0,   !   1389  1.00000
     . 51711743d0,-126884837d0,-110254279d0, 172214180d0,  50836236d0,   !   1390  1.00000
     . 51818247d0,-126744827d0,-111112349d0, 172043999d0,  51733653d0,   !   1391  1.00000
     . 51924651d0,-126604507d0,-111969761d0, 171873299d0,  52630470d0,   !   1392  1.00000
     . 52030955d0,-126463876d0,-112826513d0, 171702083d0,  53526683d0,   !   1393  1.00000
     . 52137160d0,-126322935d0,-113682601d0, 171530350d0,  54422288d0,   !   1394  1.00000
     . 52243265d0,-126181685d0,-114538021d0, 171358102d0,  55317281d0,   !   1395  1.00000
     . 52349270d0,-126040126d0,-115392770d0, 171185339d0,  56211656d0,   !   1396  1.00000
     . 52455173d0,-125898258d0,-116246844d0, 171012062d0,  57105412d0,   !   1397  1.00000
     . 52560976d0,-125756083d0,-117100240d0, 170838272d0,  57998542d0,   !   1398  1.00000
     . 52666677d0,-125613600d0,-117952954d0, 170663970d0,  58891044d0,   !   1399  1.00000
     . 52772276d0,-125470810d0,-118804983d0, 170489157d0,  59782914d0,   !   1400  1.00000
     . 52877773d0,-125327713d0,-119656323d0, 170313834d0,  60674146d0,   !   1401  1.00000
     . 52983168d0,-125184311d0,-120506971d0, 170138000d0,  61564738d0,   !   1402  1.00000
     . 53088459d0,-125040603d0,-121356924d0, 169961658d0,  62454685d0,   !   1403  1.00000
     . 53193648d0,-124896591d0,-122206178d0, 169784808d0,  63343983d0,   !   1404  1.00000
     . 53298732d0,-124752274d0,-123054729d0, 169607451d0,  64232629d0,   !   1405  1.00000
     . 53403713d0,-124607653d0,-123902574d0, 169429588d0,  65120617d0,   !   1406  1.00000
     . 53508590d0,-124462728d0,-124749710d0, 169251219d0,  66007945d0,   !   1407  1.00000
     . 53613362d0,-124317501d0,-125596134d0, 169072345d0,  66894607d0,   !   1408  1.00000
     . 53718029d0,-124171971d0,-126441841d0, 168892968d0,  67780601d0,   !   1409  1.00000
     . 53822591d0,-124026140d0,-127286828d0, 168713088d0,  68665922d0,   !   1410  1.00000
     . 53927047d0,-123880007d0,-128131093d0, 168532705d0,  69550567d0,   !   1411  1.00000
     . 54031398d0,-123733573d0,-128974631d0, 168351822d0,  70434530d0,   !   1412  1.00000
     . 54135642d0,-123586838d0,-129817440d0, 168170438d0,  71317809d0,   !   1413  1.00000
     . 54239779d0,-123439804d0,-130659515d0, 167988554d0,  72200399d0,   !   1414  1.00000
     . 54343809d0,-123292470d0,-131500854d0, 167806172d0,  73082297d0,   !   1415  1.00000
     . 54447732d0,-123144838d0,-132341452d0, 167623291d0,  73963498d0,   !   1416  1.00000
     . 54551548d0,-122996906d0,-133181308d0, 167439914d0,  74843998d0,   !   1417  1.00000
     . 54655255d0,-122848678d0,-134020417d0, 167256040d0,  75723794d0,   !   1418  1.00000
     . 54758854d0,-122700151d0,-134858775d0, 167071671d0,  76602882d0,   !   1419  1.00000
     . 54862345d0,-122551328d0,-135696381d0, 166886807d0,  77481257d0,   !   1420  1.00000
     . 54965726d0,-122402208d0,-136533229d0, 166701450d0,  78358916d0,   !   1421  1.00000
     . 55068998d0,-122252793d0,-137369318d0, 166515600d0,  79235855d0,   !   1422  1.00000
     . 55172161d0,-122103082d0,-138204643d0, 166329257d0,  80112070d0,   !   1423  1.00000
     . 55275213d0,-121953076d0,-139039201d0, 166142423d0,  80987558d0,   !   1424  1.00000
     . 55378155d0,-121802776d0,-139872989d0, 165955099d0,  81862313d0,   !   1425  1.00000
     . 55480987d0,-121652182d0,-140706004d0, 165767285d0,  82736333d0,   !   1426  1.00000
     . 55583707d0,-121501294d0,-141538242d0, 165578983d0,  83609613d0,   !   1427  1.00000
     . 55686316d0,-121350114d0,-142369699d0, 165390193d0,  84482150d0,   !   1428  1.00000
     . 55788814d0,-121198641d0,-143200374d0, 165200915d0,  85353940d0,   !   1429  1.00000
     . 55891199d0,-121046876d0,-144030262d0, 165011152d0,  86224979d0,   !   1430  1.00000
     . 55993472d0,-120894820d0,-144859360d0, 164820902d0,  87095262d0,   !   1431  1.00000
     . 56095633d0,-120742473d0,-145687664d0, 164630169d0,  87964787d0,   !   1432  1.00000
     . 56197680d0,-120589836d0,-146515172d0, 164438951d0,  88833550d0,   !   1433  1.00000
     . 56299614d0,-120436908d0,-147341880d0, 164247251d0,  89701546d0,   !   1434  1.00000
     . 56401435d0,-120283692d0,-148167786d0, 164055069d0,  90568772d0,   !   1435  1.00000
     . 56503141d0,-120130186d0,-148992884d0, 163862405d0,  91435223d0,   !   1436  1.00000
     . 56604734d0,-119976392d0,-149817173d0, 163669261d0,  92300897d0,   !   1437  1.00000
     . 56706211d0,-119822311d0,-150640650d0, 163475638d0,  93165790d0,   !   1438  1.00000
     . 56807574d0,-119667942d0,-151463310d0, 163281536d0,  94029896d0,   !   1439  1.00000
     . 56908822d0,-119513286d0,-152285151d0, 163086956d0,  94893214d0,   !   1440  1.00000
     . 57009953d0,-119358343d0,-153106169d0, 162891899d0,  95755739d0,   !   1441  1.00000
     . 57110969d0,-119203115d0,-153926361d0, 162696365d0,  96617467d0,   !   1442  1.00000
     . 57211869d0,-119047601d0,-154745724d0, 162500357d0,  97478394d0,   !   1443  1.00000
     . 57312652d0,-118891803d0,-155564255d0, 162303873d0,  98338517d0,   !   1444  1.00000
     . 57413319d0,-118735720d0,-156381950d0, 162106916d0,  99197832d0,   !   1445  1.00000
     . 57513868d0,-118579353d0,-157198807d0, 161909487d0, 100056336d0,   !   1446  1.00000
     . 57614299d0,-118422703d0,-158014822d0, 161711585d0, 100914024d0,   !   1447  1.00000
     . 57714613d0,-118265770d0,-158829991d0, 161513212d0, 101770892d0,   !   1448  1.00000
     . 57814809d0,-118108555d0,-159644313d0, 161314368d0, 102626938d0,   !   1449  1.00000
     . 57914886d0,-117951058d0,-160457782d0, 161115055d0, 103482157d0,   !   1450  1.00000
     . 58014844d0,-117793279d0,-161270397d0, 160915274d0, 104336546d0,   !   1451  1.00000
     . 58114683d0,-117635220d0,-162082155d0, 160715024d0, 105190100d0,   !   1452  1.00000
     . 58214403d0,-117476880d0,-162893051d0, 160514308d0, 106042817d0,   !   1453  1.00000
     . 58314003d0,-117318261d0,-163703083d0, 160313125d0, 106894692d0,   !   1454  1.00000
     . 58413483d0,-117159362d0,-164512247d0, 160111477d0, 107745722d0,   !   1455  1.00000
     . 58512843d0,-117000184d0,-165320541d0, 159909364d0, 108595903d0,   !   1456  1.00000
     . 58612082d0,-116840728d0,-166127962d0, 159706788d0, 109445232d0,   !   1457  1.00000
     . 58711199d0,-116680994d0,-166934505d0, 159503749d0, 110293705d0,   !   1458  1.00000
     . 58810196d0,-116520983d0,-167740169d0, 159300248d0, 111141317d0,   !   1459  1.00000
     . 58909071d0,-116360695d0,-168544949d0, 159096286d0, 111988066d0,   !   1460  1.00000
     . 59007824d0,-116200130d0,-169348844d0, 158891863d0, 112833948d0,   !   1461  1.00000
     . 59106455d0,-116039290d0,-170151849d0, 158686981d0, 113678960d0,   !   1462  1.00000
     . 59204963d0,-115878174d0,-170953961d0, 158481640d0, 114523097d0,   !   1463  1.00000
     . 59303348d0,-115716784d0,-171755178d0, 158275842d0, 115366356d0,   !   1464  1.00000
     . 59401610d0,-115555119d0,-172555497d0, 158069586d0, 116208733d0,   !   1465  1.00000
     . 59499748d0,-115393180d0,-173354913d0, 157862875d0, 117050225d0,   !   1466  1.00000
     . 59597763d0,-115230968d0,-174153425d0, 157655708d0, 117890828d0,   !   1467  1.00000
     . 59695654d0,-115068483d0,-174951029d0, 157448086d0, 118730538d0,   !   1468  1.00000
     . 59793420d0,-114905726d0,-175747722d0, 157240011d0, 119569353d0,   !   1469  1.00000
     . 59891061d0,-114742697d0,-176543500d0, 157031483d0, 120407268d0,   !   1470  1.00000
     . 59988578d0,-114579397d0,-177338362d0, 156822503d0, 121244280d0,   !   1471  1.00000
     . 60085969d0,-114415826d0,-178132303d0, 156613072d0, 122080385d0,   !   1472  1.00000
     . 60183234d0,-114251985d0,-178925321d0, 156403190d0, 122915579d0,   !   1473  1.00000
     . 60280374d0,-114087873d0,-179717413d0, 156192859d0, 123749860d0,   !   1474  1.00000
     . 60377387d0,-113923493d0,-180508575d0, 155982080d0, 124583223d0,   !   1475  1.00000
     . 60474273d0,-113758843d0,-181298805d0, 155770852d0, 125415665d0,   !   1476  1.00000
     . 60571033d0,-113593926d0,-182088100d0, 155559178d0, 126247183d0,   !   1477  1.00000
     . 60667666d0,-113428740d0,-182876456d0, 155347058d0, 127077772d0,   !   1478  1.00000
     . 60764171d0,-113263287d0,-183663870d0, 155134492d0, 127907430d0,   !   1479  1.00000
     . 60860548d0,-113097568d0,-184450340d0, 154921481d0, 128736153d0,   !   1480  1.00000
     . 60956797d0,-112931582d0,-185235863d0, 154708028d0, 129563937d0,   !   1481  1.00000
     . 61052918d0,-112765330d0,-186020435d0, 154494131d0, 130390779d0,   !   1482  1.00000
     . 61148910d0,-112598813d0,-186804053d0, 154279792d0, 131216675d0,   !   1483  1.00000
     . 61244773d0,-112432031d0,-187586715d0, 154065013d0, 132041623d0,   !   1484  1.00000
     . 61340507d0,-112264985d0,-188368417d0, 153849793d0, 132865617d0,   !   1485  1.00000
     . 61436111d0,-112097675d0,-189149157d0, 153634133d0, 133688656d0,   !   1486  1.00000
     . 61531585d0,-111930101d0,-189928932d0, 153418035d0, 134510734d0,   !   1487  1.00000
     . 61626930d0,-111762265d0,-190707738d0, 153201499d0, 135331850d0,   !   1488  1.00000
     . 61722143d0,-111594167d0,-191485572d0, 152984526d0, 136152000d0,   !   1489  1.00000
     . 61817226d0,-111425807d0,-192262432d0, 152767118d0, 136971179d0,   !   1490  1.00000
     . 61912178d0,-111257185d0,-193038315d0, 152549273d0, 137789385d0,   !   1491  1.00000
     . 62006998d0,-111088303d0,-193813217d0, 152330995d0, 138606614d0,   !   1492  1.00000
     . 62101687d0,-110919160d0,-194587136d0, 152112282d0, 139422863d0,   !   1493  1.00000
     . 62196244d0,-110749758d0,-195360070d0, 151893137d0, 140238128d0,   !   1494  1.00000
     . 62290669d0,-110580096d0,-196132013d0, 151673561d0, 141052406d0,   !   1495  1.00000
     . 62384961d0,-110410176d0,-196902965d0, 151453553d0, 141865694d0,   !   1496  1.00000
     . 62479120d0,-110239997d0,-197672922d0, 151233114d0, 142677987d0,   !   1497  1.00000
     . 62573146d0,-110069561d0,-198441881d0, 151012247d0, 143489283d0,   !   1498  1.00000
     . 62667039d0,-109898867d0,-199209840d0, 150790950d0, 144299579d0,   !   1499  1.00000
     . 62760798d0,-109727916d0,-199976794d0, 150569226d0, 145108870d0,   !   1500  1.00000
     . 62854423d0,-109556710d0,-200742742d0, 150347075d0, 145917154d0,   !   1501  1.00000
     . 62947914d0,-109385247d0,-201507681d0, 150124498d0, 146724427d0,   !   1502  1.00000
     . 63041270d0,-109213529d0,-202271607d0, 149901496d0, 147530685d0,   !   1503  1.00000
     . 63134491d0,-109041556d0,-203034518d0, 149678069d0, 148335926d0,   !   1504  1.00000
     . 63227578d0,-108869329d0,-203796410d0, 149454219d0, 149140146d0,   !   1505  1.00000
     . 63320528d0,-108696849d0,-204557282d0, 149229946d0, 149943342d0,   !   1506  1.00000
     . 63413343d0,-108524115d0,-205317130d0, 149005251d0, 150745510d0,   !   1507  1.00000
     . 63506022d0,-108351128d0,-206075951d0, 148780135d0, 151546647d0,   !   1508  1.00000
     . 63598565d0,-108177889d0,-206833743d0, 148554599d0, 152346750d0,   !   1509  1.00000
     . 63690971d0,-108004398d0,-207590502d0, 148328643d0, 153145815d0,   !   1510  1.00000
     . 63783240d0,-107830656d0,-208346226d0, 148102269d0, 153943840d0,   !   1511  1.00000
     . 63875373d0,-107656663d0,-209100911d0, 147875477d0, 154740820d0,   !   1512  1.00000
     . 63967367d0,-107482420d0,-209854556d0, 147648268d0, 155536752d0,   !   1513  1.00000
     . 64059224d0,-107307927d0,-210607157d0, 147420643d0, 156331634d0,   !   1514  1.00000
     . 64150943d0,-107133185d0,-211358712d0, 147192603d0, 157125461d0,   !   1515  1.00000
     . 64242524d0,-106958194d0,-212109217d0, 146964149d0, 157918231d0,   !   1516  1.00000
     . 64333966d0,-106782955d0,-212858670d0, 146735281d0, 158709940d0,   !   1517  1.00000
     . 64425269d0,-106607468d0,-213607068d0, 146506000d0, 159500586d0,   !   1518  1.00000
     . 64516434d0,-106431734d0,-214354408d0, 146276307d0, 160290164d0,   !   1519  1.00000
     . 64607458d0,-106255754d0,-215100687d0, 146046204d0, 161078671d0,   !   1520  1.00000
     . 64698344d0,-106079526d0,-215845904d0, 145815690d0, 161866105d0,   !   1521  1.00000
     . 64789089d0,-105903054d0,-216590054d0, 145584767d0, 162652462d0,   !   1522  1.00000
     . 64879694d0,-105726335d0,-217333135d0, 145353435d0, 163437739d0,   !   1523  1.00000
     . 64970158d0,-105549373d0,-218075144d0, 145121695d0, 164221932d0,   !   1524  1.00000
     . 65060482d0,-105372165d0,-218816079d0, 144889549d0, 165005039d0,   !   1525  1.00000
     . 65150664d0,-105194715d0,-219555936d0, 144656997d0, 165787056d0,   !   1526  1.00000
     . 65240706d0,-105017021d0,-220294714d0, 144424039d0, 166567979d0,   !   1527  1.00000
     . 65330605d0,-104839084d0,-221032409d0, 144190677d0, 167347807d0,   !   1528  1.00000
     . 65420363d0,-104660905d0,-221769019d0, 143956912d0, 168126535d0,   !   1529  1.00000
     . 65509978d0,-104482484d0,-222504540d0, 143722743d0, 168904160d0,   !   1530  1.00000
     . 65599452d0,-104303822d0,-223238970d0, 143488173d0, 169680679d0,   !   1531  1.00000
     . 65688782d0,-104124919d0,-223972307d0, 143253202d0, 170456090d0,   !   1532  1.00000
     . 65777970d0,-103945776d0,-224704547d0, 143017831d0, 171230388d0,   !   1533  1.00000
     . 65867014d0,-103766393d0,-225435689d0, 142782061d0, 172003571d0,   !   1534  1.00000
     . 65955915d0,-103586772d0,-226165728d0, 142545892d0, 172775635d0,   !   1535  1.00000
     . 66044671d0,-103406911d0,-226894663d0, 142309325d0, 173546578d0,   !   1536  1.00000
     . 66133284d0,-103226812d0,-227622491d0, 142072361d0, 174316395d0,   !   1537  1.00000
     . 66221753d0,-103046475d0,-228349209d0, 141835002d0, 175085085d0,   !   1538  1.00000
     . 66310077d0,-102865902d0,-229074814d0, 141597247d0, 175852644d0,   !   1539  1.00000
     . 66398256d0,-102685091d0,-229799305d0, 141359098d0, 176619069d0,   !   1540  1.00000
     . 66486290d0,-102504044d0,-230522677d0, 141120555d0, 177384356d0,   !   1541  1.00000
     . 66574178d0,-102322762d0,-231244929d0, 140881620d0, 178148503d0,   !   1542  1.00000
     . 66661921d0,-102141244d0,-231966057d0, 140642293d0, 178911506d0,   !   1543  1.00000
     . 66749518d0,-101959492d0,-232686060d0, 140402575d0, 179673363d0,   !   1544  1.00000
     . 66836969d0,-101777505d0,-233404934d0, 140162466d0, 180434070d0,   !   1545  1.00000
     . 66924273d0,-101595285d0,-234122677d0, 139921969d0, 181193624d0,   !   1546  1.00000
     . 67011431d0,-101412831d0,-234839286d0, 139681083d0, 181952023d0,   !   1547  1.00000
     . 67098441d0,-101230145d0,-235554758d0, 139439809d0, 182709262d0,   !   1548  1.00000
     . 67185305d0,-101047226d0,-236269092d0, 139198148d0, 183465340d0,   !   1549  1.00000
     . 67272021d0,-100864076d0,-236982284d0, 138956102d0, 184220252d0,   !   1550  1.00000
     . 67358589d0,-100680695d0,-237694331d0, 138713670d0, 184973996d0,   !   1551  1.00000
     . 67445009d0,-100497083d0,-238405232d0, 138470854d0, 185726570d0,   !   1552  1.00000
     . 67531280d0,-100313240d0,-239114983d0, 138227655d0, 186477969d0,   !   1553  1.00000
     . 67617404d0,-100129168d0,-239823582d0, 137984073d0, 187228190d0,   !   1554  1.00000
     . 67703378d0, -99944867d0,-240531026d0, 137740109d0, 187977232d0,   !   1555  1.00000
     . 67789203d0, -99760337d0,-241237313d0, 137495764d0, 188725090d0,   !   1556  1.00000
     . 67874879d0, -99575579d0,-241942440d0, 137251039d0, 189471762d0,   !   1557  1.00000
     . 67960406d0, -99390594d0,-242646405d0, 137005935d0, 190217245d0,   !   1558  1.00000
     . 68045782d0, -99205381d0,-243349204d0, 136760452d0, 190961535d0,   !   1559  1.00000
     . 68131008d0, -99019941d0,-244050836d0, 136514592d0, 191704630d0,   !   1560  1.00000
     . 68216085d0, -98834275d0,-244751298d0, 136268355d0, 192446527d0,   !   1561  1.00000
     . 68301010d0, -98648384d0,-245450587d0, 136021741d0, 193187223d0,   !   1562  1.00000
     . 68385785d0, -98462267d0,-246148701d0, 135774753d0, 193926715d0,   !   1563  1.00000
     . 68470408d0, -98275926d0,-246845637d0, 135527391d0, 194664999d0,   !   1564  1.00000
     . 68554880d0, -98089360d0,-247541393d0, 135279655d0, 195402073d0,   !   1565  1.00000
     . 68639200d0, -97902571d0,-248235965d0, 135031546d0, 196137934d0,   !   1566  1.00000
     . 68723369d0, -97715558d0,-248929353d0, 134783066d0, 196872580d0,   !   1567  1.00000
     . 68807385d0, -97528323d0,-249621552d0, 134534214d0, 197606006d0,   !   1568  1.00000
     . 68891249d0, -97340866d0,-250312561d0, 134284993d0, 198338210d0,   !   1569  1.00000
     . 68974961d0, -97153187d0,-251002378d0, 134035402d0, 199069190d0,   !   1570  1.00000
     . 69058519d0, -96965287d0,-251690998d0, 133785443d0, 199798942d0,   !   1571  1.00000
     . 69141925d0, -96777166d0,-252378421d0, 133535116d0, 200527463d0,   !   1572  1.00000
     . 69225177d0, -96588825d0,-253064643d0, 133284423d0, 201254750d0,   !   1573  1.00000
     . 69308275d0, -96400265d0,-253749663d0, 133033364d0, 201980802d0,   !   1574  1.00000
     . 69391220d0, -96211485d0,-254433477d0, 132781939d0, 202705614d0,   !   1575  1.00000
     . 69474010d0, -96022486d0,-255116083d0, 132530150d0, 203429183d0,   !   1576  1.00000
     . 69556646d0, -95833270d0,-255797479d0, 132277998d0, 204151508d0,   !   1577  1.00000
     . 69639127d0, -95643835d0,-256477662d0, 132025483d0, 204872585d0,   !   1578  1.00000
     . 69721454d0, -95454184d0,-257156630d0, 131772607d0, 205592410d0,   !   1579  1.00000
     . 69803625d0, -95264316d0,-257834380d0, 131519369d0, 206310983d0,   !   1580  1.00000
     . 69885641d0, -95074232d0,-258510910d0, 131265772d0, 207028298d0,   !   1581  1.00000
     . 69967501d0, -94883932d0,-259186217d0, 131011815d0, 207744354d0,   !   1582  1.00000
     . 70049206d0, -94693417d0,-259860300d0, 130757500d0, 208459148d0,   !   1583  1.00000
     . 70130754d0, -94502687d0,-260533155d0, 130502827d0, 209172677d0,   !   1584  1.00000
     . 70212146d0, -94311743d0,-261204780d0, 130247798d0, 209884939d0,   !   1585  1.00000
     . 70293382d0, -94120585d0,-261875172d0, 129992412d0, 210595929d0,   !   1586  1.00000
     . 70374461d0, -93929215d0,-262544330d0, 129736672d0, 211305646d0,   !   1587  1.00000
     . 70455382d0, -93737631d0,-263212251d0, 129480577d0, 212014086d0,   !   1588  1.00000
     . 70536147d0, -93545836d0,-263878933d0, 129224129d0, 212721248d0,   !   1589  1.00000
     . 70616753d0, -93353829d0,-264544372d0, 128967328d0, 213427127d0,   !   1590  1.00000
     . 70697202d0, -93161611d0,-265208568d0, 128710176d0, 214131722d0,   !   1591  1.00000
     . 70777494d0, -92969182d0,-265871516d0, 128452673d0, 214835030d0,   !   1592  1.00000
     . 70857626d0, -92776543d0,-266533216d0, 128194819d0, 215537047d0,   !   1593  1.00000
     . 70937601d0, -92583694d0,-267193664d0, 127936617d0, 216237771d0,   !   1594  1.00000
     . 71017416d0, -92390636d0,-267852858d0, 127678066d0, 216937199d0,   !   1595  1.00000
     . 71097073d0, -92197370d0,-268510796d0, 127419167d0, 217635329d0,   !   1596  1.00000
     . 71176570d0, -92003895d0,-269167476d0, 127159922d0, 218332158d0,   !   1597  1.00000
     . 71255908d0, -91810213d0,-269822895d0, 126900331d0, 219027683d0,   !   1598  1.00000
     . 71335086d0, -91616324d0,-270477050d0, 126640395d0, 219721900d0,   !   1599  1.00000
     . 71414105d0, -91422228d0,-271129940d0, 126380114d0, 220414809d0,   !   1600  1.00000
     . 71492963d0, -91227926d0,-271781562d0, 126119491d0, 221106405d0,   !   1601  1.00000
     . 71571661d0, -91033418d0,-272431914d0, 125858524d0, 221796686d0,   !   1602  1.00000
     . 71650198d0, -90838705d0,-273080994d0, 125597217d0, 222485650d0,   !   1603  1.00000
     . 71728575d0, -90643788d0,-273728799d0, 125335568d0, 223173294d0,   !   1604  1.00000
     . 71806790d0, -90448666d0,-274375326d0, 125073579d0, 223859614d0,   !   1605  1.00000
     . 71884844d0, -90253341d0,-275020574d0, 124811251d0, 224544609d0,   !   1606  1.00000
     . 71962737d0, -90057812d0,-275664540d0, 124548584d0, 225228275d0,   !   1607  1.00000
     . 72040468d0, -89862082d0,-276307223d0, 124285580d0, 225910611d0,   !   1608  1.00000
     . 72118036d0, -89666148d0,-276948619d0, 124022240d0, 226591612d0,   !   1609  1.00000
     . 72195443d0, -89470014d0,-277588726d0, 123758563d0, 227271278d0,   !   1610  1.00000
     . 72272687d0, -89273678d0,-278227542d0, 123494552d0, 227949604d0,   !   1611  1.00000
     . 72349769d0, -89077142d0,-278865066d0, 123230206d0, 228626589d0,   !   1612  1.00000
     . 72426687d0, -88880405d0,-279501294d0, 122965527d0, 229302229d0,   !   1613  1.00000
     . 72503443d0, -88683469d0,-280136224d0, 122700516d0, 229976522d0,   !   1614  1.00000
     . 72580035d0, -88486334d0,-280769854d0, 122435172d0, 230649466d0,   !   1615  1.00000
     . 72656464d0, -88289000d0,-281402182d0, 122169498d0, 231321058d0,   !   1616  1.00000
     . 72732729d0, -88091468d0,-282033206d0, 121903494d0, 231991295d0,   !   1617  1.00000
     . 72808829d0, -87893739d0,-282662923d0, 121637161d0, 232660174d0,   !   1618  1.00000
     . 72884766d0, -87695812d0,-283291331d0, 121370500d0, 233327694d0,   !   1619  1.00000
     . 72960538d0, -87497689d0,-283918428d0, 121103511d0, 233993850d0,   !   1620  1.00000
     . 73036146d0, -87299369d0,-284544212d0, 120836195d0, 234658642d0,   !   1621  1.00000
     . 73111588d0, -87100855d0,-285168680d0, 120568553d0, 235322066d0,   !   1622  1.00000
     . 73186865d0, -86902145d0,-285791831d0, 120300587d0, 235984120d0,   !   1623  1.00000
     . 73261978d0, -86703240d0,-286413662d0, 120032296d0, 236644800d0,   !   1624  1.00000
     . 73336924d0, -86504141d0,-287034170d0, 119763682d0, 237304105d0,   !   1625  1.00000
     . 73411705d0, -86304849d0,-287653354d0, 119494745d0, 237962033d0,   !   1626  1.00000
     . 73486319d0, -86105364d0,-288271212d0, 119225487d0, 238618579d0,   !   1627  1.00000
     . 73560768d0, -85905686d0,-288887740d0, 118955908d0, 239273743d0,   !   1628  1.00000
     . 73635050d0, -85705816d0,-289502938d0, 118686009d0, 239927520d0,   !   1629  1.00000
     . 73709165d0, -85505754d0,-290116803d0, 118415791d0, 240579910d0,   !   1630  1.00000
     . 73783114d0, -85305501d0,-290729333d0, 118145255d0, 241230909d0,   !   1631  1.00000
     . 73856895d0, -85105058d0,-291340525d0, 117874401d0, 241880515d0,   !   1632  1.00000
     . 73930509d0, -84904425d0,-291950378d0, 117603231d0, 242528725d0,   !   1633  1.00000
     . 74003956d0, -84703602d0,-292558890d0, 117331744d0, 243175537d0,   !   1634  1.00000
     . 74077235d0, -84502590d0,-293166057d0, 117059943d0, 243820948d0,   !   1635  1.00000
     . 74150346d0, -84301389d0,-293771878d0, 116787828d0, 244464956d0,   !   1636  1.00000
     . 74223289d0, -84100001d0,-294376352d0, 116515399d0, 245107558d0,   !   1637  1.00000
     . 74296063d0, -83898424d0,-294979475d0, 116242658d0, 245748752d0,   !   1638  1.00000
     . 74368669d0, -83696661d0,-295581246d0, 115969605d0, 246388535d0,   !   1639  1.00000
     . 74441106d0, -83494711d0,-296181663d0, 115696241d0, 247026906d0,   !   1640  1.00000
     . 74513374d0, -83292575d0,-296780723d0, 115422568d0, 247663860d0,   !   1641  1.00000
     . 74585473d0, -83090254d0,-297378424d0, 115148585d0, 248299397d0,   !   1642  1.00000
     . 74657402d0, -82887747d0,-297974765d0, 114874294d0, 248933513d0,   !   1643  1.00000
     . 74729162d0, -82685056d0,-298569743d0, 114599696d0, 249566207d0,   !   1644  1.00000
     . 74800752d0, -82482181d0,-299163356d0, 114324791d0, 250197475d0,   !   1645  1.00000
     . 74872172d0, -82279122d0,-299755602d0, 114049580d0, 250827315d0,   !   1646  1.00000
     . 74943421d0, -82075880d0,-300346479d0, 113774064d0, 251455726d0,   !   1647  1.00000
     . 75014500d0, -81872456d0,-300935985d0, 113498244d0, 252082704d0,   !   1648  1.00000
     . 75085409d0, -81668850d0,-301524118d0, 113222121d0, 252708246d0,   !   1649  1.00000
     . 75156146d0, -81465062d0,-302110876d0, 112945696d0, 253332352d0,   !   1650  1.00000
     . 75226713d0, -81261093d0,-302696256d0, 112668969d0, 253955018d0,   !   1651  1.00000
     . 75297108d0, -81056943d0,-303280257d0, 112391941d0, 254576241d0,   !   1652  1.00000
     . 75367332d0, -80852614d0,-303862877d0, 112114613d0, 255196020d0,   !   1653  1.00000
     . 75437384d0, -80648104d0,-304444113d0, 111836986d0, 255814353d0,   !   1654  1.00000
     . 75507264d0, -80443416d0,-305023964d0, 111559061d0, 256431236d0,   !   1655  1.00000
     . 75576972d0, -80238550d0,-305602428d0, 111280838d0, 257046667d0,   !   1656  1.00000
     . 75646508d0, -80033505d0,-306179502d0, 111002319d0, 257660645d0,   !   1657  1.00000
     . 75715871d0, -79828282d0,-306755185d0, 110723504d0, 258273166d0,   !   1658  1.00000
     . 75785062d0, -79622883d0,-307329474d0, 110444394d0, 258884228d0,   !   1659  1.00000
     . 75854080d0, -79417307d0,-307902368d0, 110164990d0, 259493830d0,   !   1660  1.00000
     . 75922924d0, -79211555d0,-308473864d0, 109885293d0, 260101969d0,   !   1661  1.00000
     . 75991595d0, -79005628d0,-309043961d0, 109605304d0, 260708641d0,   !   1662  1.00000
     . 76060093d0, -78799525d0,-309612657d0, 109325022d0, 261313846d0,   !   1663  1.00000
     . 76128417d0, -78593248d0,-310179949d0, 109044450d0, 261917581d0,   !   1664  1.00000
     . 76196568d0, -78386797d0,-310745836d0, 108763589d0, 262519843d0,   !   1665  1.00000
     . 76264544d0, -78180173d0,-311310315d0, 108482438d0, 263120631d0,   !   1666  1.00000
     . 76332345d0, -77973375d0,-311873386d0, 108200999d0, 263719941d0,   !   1667  1.00000
     . 76399973d0, -77766405d0,-312435045d0, 107919272d0, 264317772d0,   !   1668  1.00000
     . 76467426d0, -77559263d0,-312995291d0, 107637259d0, 264914122d0,   !   1669  1.00000
     . 76534703d0, -77351949d0,-313554121d0, 107354960d0, 265508988d0,   !   1670  1.00000
     . 76601806d0, -77144464d0,-314111535d0, 107072377d0, 266102367d0,   !   1671  1.00000
     . 76668734d0, -76936809d0,-314667530d0, 106789509d0, 266694258d0,   !   1672  1.00000
     . 76735486d0, -76728984d0,-315222104d0, 106506358d0, 267284659d0,   !   1673  1.00000
     . 76802062d0, -76520989d0,-315775255d0, 106222925d0, 267873566d0,   !   1674  1.00000
     . 76868463d0, -76312826d0,-316326981d0, 105939210d0, 268460979d0,   !   1675  1.00000
     . 76934687d0, -76104493d0,-316877281d0, 105655214d0, 269046894d0,   !   1676  1.00000
     . 77000736d0, -75895993d0,-317426152d0, 105370939d0, 269631310d0,   !   1677  1.00000
     . 77066608d0, -75687325d0,-317973593d0, 105086384d0, 270214224d0,   !   1678  1.00000
     . 77132303d0, -75478490d0,-318519601d0, 104801551d0, 270795634d0,   !   1679  1.00000
     . 77197822d0, -75269489d0,-319064176d0, 104516441d0, 271375537d0,   !   1680  1.00000
     . 77263163d0, -75060322d0,-319607314d0, 104231055d0, 271953933d0,   !   1681  1.00000
     . 77328328d0, -74850989d0,-320149014d0, 103945392d0, 272530818d0,   !   1682  1.00000
     . 77393315d0, -74641491d0,-320689274d0, 103659455d0, 273106190d0,   !   1683  1.00000
     . 77458125d0, -74431829d0,-321228093d0, 103373244d0, 273680048d0,   !   1684  1.00000
     . 77522756d0, -74222003d0,-321765468d0, 103086759d0, 274252388d0,   !   1685  1.00000
     . 77587210d0, -74012013d0,-322301398d0, 102800003d0, 274823209d0,   !   1686  1.00000
     . 77651486d0, -73801860d0,-322835880d0, 102512974d0, 275392509d0,   !   1687  1.00000
     . 77715584d0, -73591545d0,-323368914d0, 102225675d0, 275960285d0,   !   1688  1.00000
     . 77779503d0, -73381067d0,-323900496d0, 101938107d0, 276526536d0,   !   1689  1.00000
     . 77843244d0, -73170429d0,-324430626d0, 101650269d0, 277091259d0,   !   1690  1.00000
     . 77906806d0, -72959629d0,-324959301d0, 101362163d0, 277654452d0,   !   1691  1.00000
     . 77970188d0, -72748668d0,-325486520d0, 101073789d0, 278216113d0,   !   1692  1.00000
     . 78033392d0, -72537548d0,-326012281d0, 100785149d0, 278776240d0,   !   1693  1.00000
     . 78096416d0, -72326269d0,-326536581d0, 100496244d0, 279334831d0,   !   1694  1.00000
     . 78159261d0, -72114830d0,-327059420d0, 100207073d0, 279891884d0,   !   1695  1.00000
     . 78221926d0, -71903233d0,-327580795d0,  99917639d0, 280447396d0,   !   1696  1.00000
     . 78284411d0, -71691477d0,-328100705d0,  99627941d0, 281001365d0,   !   1697  1.00000
     . 78346716d0, -71479565d0,-328619148d0,  99337981d0, 281553790d0,   !   1698  1.00000
     . 78408841d0, -71267495d0,-329136122d0,  99047760d0, 282104669d0,   !   1699  1.00000
     . 78470785d0, -71055269d0,-329651625d0,  98757278d0, 282653998d0,   !   1700  1.00000
     . 78532549d0, -70842887d0,-330165655d0,  98466536d0, 283201777d0,   !   1701  1.00000
     . 78594132d0, -70630349d0,-330678212d0,  98175535d0, 283748003d0,   !   1702  1.00000
     . 78655534d0, -70417657d0,-331189292d0,  97884276d0, 284292674d0,   !   1703  1.00000
     . 78716755d0, -70204810d0,-331698895d0,  97592760d0, 284835789d0,   !   1704  1.00000
     . 78777794d0, -69991809d0,-332207018d0,  97300987d0, 285377344d0,   !   1705  1.00000
     . 78838652d0, -69778654d0,-332713660d0,  97008958d0, 285917339d0,   !   1706  1.00000
     . 78899329d0, -69565347d0,-333218819d0,  96716675d0, 286455770d0,   !   1707  1.00000
     . 78959823d0, -69351887d0,-333722493d0,  96424138d0, 286992636d0,   !   1708  1.00000
     . 79020136d0, -69138275d0,-334224681d0,  96131347d0, 287527936d0,   !   1709  1.00000
     . 79080266d0, -68924512d0,-334725381d0,  95838304d0, 288061666d0,   !   1710  1.00000
     . 79140215d0, -68710598d0,-335224591d0,  95545010d0, 288593826d0,   !   1711  1.00000
     . 79199980d0, -68496533d0,-335722310d0,  95251465d0, 289124412d0,   !   1712  1.00000
     . 79259563d0, -68282318d0,-336218535d0,  94957671d0, 289653424d0,   !   1713  1.00000
     . 79318963d0, -68067954d0,-336713265d0,  94663627d0, 290180858d0,   !   1714  1.00000
     . 79378180d0, -67853441d0,-337206499d0,  94369335d0, 290706714d0,   !   1715  1.00000
     . 79437214d0, -67638780d0,-337698235d0,  94074796d0, 291230989d0,   !   1716  1.00000
     . 79496064d0, -67423970d0,-338188470d0,  93780011d0, 291753681d0,   !   1717  1.00000
     . 79554731d0, -67209013d0,-338677204d0,  93484980d0, 292274788d0,   !   1718  1.00000
     . 79613215d0, -66993909d0,-339164435d0,  93189704d0, 292794309d0,   !   1719  1.00000
     . 79671514d0, -66778659d0,-339650161d0,  92894184d0, 293312241d0,   !   1720  1.00000
     . 79729629d0, -66563263d0,-340134380d0,  92598421d0, 293828582d0,   !   1721  1.00000
     . 79787561d0, -66347721d0,-340617091d0,  92302416d0, 294343331d0,   !   1722  1.00000
     . 79845308d0, -66132035d0,-341098292d0,  92006169d0, 294856486d0,   !   1723  1.00000
     . 79902870d0, -65916203d0,-341577981d0,  91709681d0, 295368044d0,   !   1724  1.00000
     . 79960248d0, -65700228d0,-342056158d0,  91412954d0, 295878004d0,   !   1725  1.00000
     . 80017440d0, -65484110d0,-342532819d0,  91115988d0, 296386364d0,   !   1726  1.00000
     . 80074448d0, -65267849d0,-343007964d0,  90818784d0, 296893122d0,   !   1727  1.00000
     . 80131271d0, -65051445d0,-343481592d0,  90521343d0, 297398276d0,   !   1728  1.00000
     . 80187908d0, -64834899d0,-343953699d0,  90223665d0, 297901825d0,   !   1729  1.00000
     . 80244360d0, -64618212d0,-344424286d0,  89925751d0, 298403766d0,   !   1730  1.00000
     . 80300626d0, -64401384d0,-344893349d0,  89627603d0, 298904097d0,   !   1731  1.00000
     . 80356707d0, -64184415d0,-345360888d0,  89329221d0, 299402817d0,   !   1732  1.00000
     . 80412601d0, -63967307d0,-345826902d0,  89030606d0, 299899924d0,   !   1733  1.00000
     . 80468310d0, -63750059d0,-346291387d0,  88731758d0, 300395416d0,   !   1734  1.00000
     . 80523832d0, -63532672d0,-346754344d0,  88432679d0, 300889291d0,   !   1735  1.00000
     . 80579167d0, -63315147d0,-347215770d0,  88133370d0, 301381548d0,   !   1736  1.00000
     . 80634317d0, -63097484d0,-347675663d0,  87833831d0, 301872184d0,   !   1737  1.00000
     . 80689279d0, -62879683d0,-348134023d0,  87534062d0, 302361198d0,   !   1738  1.00000
     . 80744054d0, -62661746d0,-348590848d0,  87234066d0, 302848587d0,   !   1739  1.00000
     . 80798643d0, -62443672d0,-349046135d0,  86933843d0, 303334351d0,   !   1740  1.00000
     . 80853044d0, -62225462d0,-349499884d0,  86633393d0, 303818487d0,   !   1741  1.00000
     . 80907258d0, -62007117d0,-349952094d0,  86332717d0, 304300993d0,   !   1742  1.00000
     . 80961284d0, -61788636d0,-350402761d0,  86031817d0, 304781869d0,   !   1743  1.00000
     . 81015123d0, -61570022d0,-350851886d0,  85730693d0, 305261111d0,   !   1744  1.00000
     . 81068774d0, -61351273d0,-351299466d0,  85429346d0, 305738718d0,   !   1745  1.00000
     . 81122237d0, -61132392d0,-351745500d0,  85127776d0, 306214689d0,   !   1746  1.00000
     . 81175512d0, -60913377d0,-352189986d0,  84825985d0, 306689022d0,   !   1747  1.00000
     . 81228598d0, -60694230d0,-352632924d0,  84523974d0, 307161714d0,   !   1748  1.00000
     . 81281497d0, -60474951d0,-353074311d0,  84221743d0, 307632765d0,   !   1749  1.00000
     . 81334206d0, -60255540d0,-353514145d0,  83919293d0, 308102172d0,   !   1750  1.00000
     . 81386727d0, -60035999d0,-353952427d0,  83616625d0, 308569934d0,   !   1751  1.00000
     . 81439059d0, -59816327d0,-354389153d0,  83313740d0, 309036049d0,   !   1752  1.00000
     . 81491202d0, -59596525d0,-354824323d0,  83010638d0, 309500515d0,   !   1753  1.00000
     . 81543155d0, -59376594d0,-355257935d0,  82707321d0, 309963331d0,   !   1754  1.00000
     . 81594920d0, -59156534d0,-355689987d0,  82403789d0, 310424494d0,   !   1755  1.00000
     . 81646495d0, -58936346d0,-356120479d0,  82100043d0, 310884004d0,   !   1756  1.00000
     . 81697880d0, -58716030d0,-356549408d0,  81796084d0, 311341858d0,   !   1757  1.00000
     . 81749076d0, -58495586d0,-356976774d0,  81491913d0, 311798055d0,   !   1758  1.00000
     . 81800082d0, -58275016d0,-357402574d0,  81187531d0, 312252593d0,   !   1759  1.00000
     . 81850897d0, -58054319d0,-357826808d0,  80882938d0, 312705471d0,   !   1760  1.00000
     . 81901523d0, -57833496d0,-358249474d0,  80578135d0, 313156687d0,   !   1761  1.00000
     . 81951958d0, -57612548d0,-358670570d0,  80273124d0, 313606238d0,   !   1762  1.00000
     . 82002203d0, -57391475d0,-359090096d0,  79967904d0, 314054124d0,   !   1763  1.00000
     . 82052257d0, -57170277d0,-359508049d0,  79662478d0, 314500344d0,   !   1764  1.00000
     . 82102121d0, -56948956d0,-359924429d0,  79356844d0, 314944894d0,   !   1765  1.00000
     . 82151793d0, -56727512d0,-360339233d0,  79051006d0, 315387774d0,   !   1766  1.00000
     . 82201275d0, -56505944d0,-360752461d0,  78744963d0, 315828981d0,   !   1767  1.00000
     . 82250566d0, -56284254d0,-361164111d0,  78438715d0, 316268516d0,   !   1768  1.00000
     . 82299665d0, -56062442d0,-361574182d0,  78132265d0, 316706375d0,   !   1769  1.00000
     . 82348573d0, -55840509d0,-361982673d0,  77825613d0, 317142557d0,   !   1770  1.00000
     . 82397289d0, -55618455d0,-362389581d0,  77518759d0, 317577060d0,   !   1771  1.00000
     . 82445813d0, -55396281d0,-362794906d0,  77211705d0, 318009884d0,   !   1772  1.00000
     . 82494146d0, -55173987d0,-363198646d0,  76904451d0, 318441026d0,   !   1773  1.00000
     . 82542287d0, -54951573d0,-363600801d0,  76596998d0, 318870484d0,   !   1774  1.00000
     . 82590235d0, -54729041d0,-364001368d0,  76289348d0, 319298258d0,   !   1775  1.00000
     . 82637992d0, -54506390d0,-364400346d0,  75981500d0, 319724346d0,   !   1776  1.00000
     . 82685556d0, -54283621d0,-364797734d0,  75673456d0, 320148746d0,   !   1777  1.00000
     . 82732927d0, -54060735d0,-365193531d0,  75365216d0, 320571456d0,   !   1778  1.00000
     . 82780106d0, -53837732d0,-365587735d0,  75056782d0, 320992476d0,   !   1779  1.00000
     . 82827092d0, -53614613d0,-365980345d0,  74748153d0, 321411802d0,   !   1780  1.00000
     . 82873885d0, -53391378d0,-366371359d0,  74439332d0, 321829435d0,   !   1781  1.00000
     . 82920486d0, -53168027d0,-366760777d0,  74130319d0, 322245372d0,   !   1782  1.00000
     . 82966893d0, -52944562d0,-367148597d0,  73821114d0, 322659612d0,   !   1783  1.00000
     . 83013106d0, -52720982d0,-367534818d0,  73511719d0, 323072154d0,   !   1784  1.00000
     . 83059127d0, -52497289d0,-367919438d0,  73202134d0, 323482995d0,   !   1785  1.00000
     . 83104953d0, -52273482d0,-368302457d0,  72892361d0, 323892135d0,   !   1786  1.00000
     . 83150586d0, -52049562d0,-368683872d0,  72582399d0, 324299572d0,   !   1787  1.00000
     . 83196026d0, -51825530d0,-369063683d0,  72272250d0, 324705304d0,   !   1788  1.00000
     . 83241271d0, -51601386d0,-369441889d0,  71961915d0, 325109330d0,   !   1789  1.00000
     . 83286322d0, -51377131d0,-369818487d0,  71651395d0, 325511649d0,   !   1790  1.00000
     . 83331179d0, -51152765d0,-370193477d0,  71340690d0, 325912258d0,   !   1791  1.00000
     . 83375842d0, -50928289d0,-370566858d0,  71029801d0, 326311158d0,   !   1792  1.00000
     . 83420311d0, -50703703d0,-370938628d0,  70718729d0, 326708345d0,   !   1793  1.00000
     . 83464584d0, -50479007d0,-371308787d0,  70407475d0, 327103819d0,   !   1794  1.00000
     . 83508663d0, -50254203d0,-371677332d0,  70096040d0, 327497578d0,   !   1795  1.00000
     . 83552548d0, -50029290d0,-372044263d0,  69784424d0, 327889621d0,   !   1796  1.00000
     . 83596237d0, -49804269d0,-372409578d0,  69472628d0, 328279946d0,   !   1797  1.00000
     . 83639732d0, -49579142d0,-372773276d0,  69160654d0, 328668553d0,   !   1798  1.00000
     . 83683031d0, -49353907d0,-373135357d0,  68848502d0, 329055438d0,   !   1799  1.00000
     . 83726135d0, -49128566d0,-373495818d0,  68536173d0, 329440602d0,   !   1800  1.00000
     . 83769043d0, -48903119d0,-373854659d0,  68223667d0, 329824043d0,   !   1801  1.00000
     . 83811756d0, -48677567d0,-374211878d0,  67910986d0, 330205759d0,   !   1802  1.00000
     . 83854274d0, -48451910d0,-374567475d0,  67598131d0, 330585749d0,   !   1803  1.00000
     . 83896595d0, -48226149d0,-374921447d0,  67285101d0, 330964012d0,   !   1804  1.00000
     . 83938721d0, -48000284d0,-375273795d0,  66971899d0, 331340546d0,   !   1805  1.00000
     . 83980651d0, -47774315d0,-375624516d0,  66658525d0, 331715349d0,   !   1806  1.00000
     . 84022384d0, -47548244d0,-375973610d0,  66344979d0, 332088421d0,   !   1807  1.00000
     . 84063922d0, -47322070d0,-376321075d0,  66031263d0, 332459761d0,   !   1808  1.00000
     . 84105263d0, -47095795d0,-376666910d0,  65717377d0, 332829366d0,   !   1809  1.00000
     . 84146407d0, -46869418d0,-377011114d0,  65403323d0, 333197235d0,   !   1810  1.00000
     . 84187355d0, -46642941d0,-377353687d0,  65089100d0, 333563368d0,   !   1811  1.00000
     . 84228106d0, -46416363d0,-377694626d0,  64774711d0, 333927763d0,   !   1812  1.00000
     . 84268661d0, -46189686d0,-378033930d0,  64460155d0, 334290417d0,   !   1813  1.00000
     . 84309018d0, -45962909d0,-378371599d0,  64145434d0, 334651332d0,   !   1814  1.00000
     . 84349179d0, -45736033d0,-378707632d0,  63830548d0, 335010504d0,   !   1815  1.00000
     . 84389142d0, -45509059d0,-379042027d0,  63515499d0, 335367932d0,   !   1816  1.00000
     . 84428909d0, -45281987d0,-379374783d0,  63200286d0, 335723616d0,   !   1817  1.00000
     . 84468477d0, -45054818d0,-379705899d0,  62884912d0, 336077554d0,   !   1818  1.00000
     . 84507849d0, -44827552d0,-380035375d0,  62569376d0, 336429744d0,   !   1819  1.00000
     . 84547023d0, -44600190d0,-380363208d0,  62253680d0, 336780186d0,   !   1820  1.00000
     . 84585999d0, -44372732d0,-380689398d0,  61937824d0, 337128879d0,   !   1821  1.00000
     . 84624777d0, -44145179d0,-381013943d0,  61621810d0, 337475820d0,   !   1822  1.00000
     . 84663357d0, -43917531d0,-381336844d0,  61305637d0, 337821008d0,   !   1823  1.00000
     . 84701740d0, -43689789d0,-381658098d0,  60989308d0, 338164444d0,   !   1824  1.00000
     . 84739924d0, -43461953d0,-381977704d0,  60672822d0, 338506124d0,   !   1825  1.00000
     . 84777910d0, -43234023d0,-382295662d0,  60356181d0, 338846048d0,   !   1826  1.00000
     . 84815698d0, -43006001d0,-382611970d0,  60039385d0, 339184215d0,   !   1827  1.00000
     . 84853287d0, -42777887d0,-382926628d0,  59722436d0, 339520624d0,   !   1828  1.00000
     . 84890678d0, -42549680d0,-383239634d0,  59405334d0, 339855272d0,   !   1829  1.00000
     . 84927870d0, -42321383d0,-383550988d0,  59088080d0, 340188160d0,   !   1830  1.00000
     . 84964863d0, -42092995d0,-383860687d0,  58770674d0, 340519286d0,   !   1831  1.00000
     . 85001658d0, -41864516d0,-384168732d0,  58453118d0, 340848648d0,   !   1832  1.00000
     . 85038253d0, -41635948d0,-384475121d0,  58135413d0, 341176246d0,   !   1833  1.00000
     . 85074650d0, -41407290d0,-384779854d0,  57817559d0, 341502079d0,   !   1834  1.00000
     . 85110847d0, -41178544d0,-385082928d0,  57499557d0, 341826144d0,   !   1835  1.00000
     . 85146846d0, -40949709d0,-385384344d0,  57181408d0, 342148441d0,   !   1836  1.00000
     . 85182644d0, -40720787d0,-385684100d0,  56863112d0, 342468969d0,   !   1837  1.00000
     . 85218244d0, -40491777d0,-385982196d0,  56544672d0, 342787727d0,   !   1838  1.00000
     . 85253644d0, -40262681d0,-386278629d0,  56226087d0, 343104714d0,   !   1839  1.00000
     . 85288844d0, -40033498d0,-386573400d0,  55907358d0, 343419927d0,   !   1840  1.00000
     . 85323844d0, -39804229d0,-386866507d0,  55588486d0, 343733367d0,   !   1841  1.00000
     . 85358645d0, -39574876d0,-387157949d0,  55269472d0, 344045033d0,   !   1842  1.00000
     . 85393246d0, -39345437d0,-387447726d0,  54950317d0, 344354922d0,   !   1843  1.00000
     . 85427646d0, -39115915d0,-387735836d0,  54631022d0, 344663034d0,   !   1844  1.00000
     . 85461847d0, -38886308d0,-388022279d0,  54311587d0, 344969368d0,   !   1845  1.00000
     . 85495847d0, -38656618d0,-388307053d0,  53992014d0, 345273922d0,   !   1846  1.00000
     . 85529647d0, -38426846d0,-388590158d0,  53672302d0, 345576696d0,   !   1847  1.00000
     . 85563247d0, -38196991d0,-388871593d0,  53352454d0, 345877689d0,   !   1848  1.00000
     . 85596646d0, -37967055d0,-389151356d0,  53032470d0, 346176899d0,   !   1849  1.00000
     . 85629844d0, -37737037d0,-389429447d0,  52712350d0, 346474325d0,   !   1850  1.00000
     . 85662842d0, -37506938d0,-389705865d0,  52392095d0, 346769967d0,   !   1851  1.00000
     . 85695639d0, -37276760d0,-389980609d0,  52071707d0, 347063822d0,   !   1852  1.00000
     . 85728235d0, -37046501d0,-390253678d0,  51751186d0, 347355891d0,   !   1853  1.00000
     . 85760631d0, -36816164d0,-390525071d0,  51430533d0, 347646172d0,   !   1854  1.00000
     . 85792825d0, -36585747d0,-390794787d0,  51109749d0, 347934665d0,   !   1855  1.00000
     . 85824818d0, -36355252d0,-391062826d0,  50788835d0, 348221367d0,   !   1856  1.00000
     . 85856610d0, -36124680d0,-391329186d0,  50467791d0, 348506278d0,   !   1857  1.00000
     . 85888200d0, -35894030d0,-391593867d0,  50146618d0, 348789397d0,   !   1858  1.00000
     . 85919590d0, -35663304d0,-391856868d0,  49825318d0, 349070723d0,   !   1859  1.00000
     . 85950777d0, -35432501d0,-392118188d0,  49503891d0, 349350256d0,   !   1860  1.00000
     . 85981763d0, -35201623d0,-392377826d0,  49182337d0, 349627993d0,   !   1861  1.00000
     . 86012548d0, -34970669d0,-392635781d0,  48860658d0, 349903934d0,   !   1862  1.00000
     . 86043131d0, -34739641d0,-392892052d0,  48538855d0, 350178078d0,   !   1863  1.00000
     . 86073512d0, -34508539d0,-393146639d0,  48216928d0, 350450424d0,   !   1864  1.00000
     . 86103691d0, -34277362d0,-393399540d0,  47894878d0, 350720971d0,   !   1865  1.00000
     . 86133668d0, -34046113d0,-393650755d0,  47572706d0, 350989718d0,   !   1866  1.00000
     . 86163443d0, -33814791d0,-393900284d0,  47250413d0, 351256665d0,   !   1867  1.00000
     . 86193016d0, -33583396d0,-394148124d0,  46928000d0, 351521809d0,   !   1868  1.00000
     . 86222386d0, -33351930d0,-394394276d0,  46605468d0, 351785151d0,   !   1869  1.00000
     . 86251554d0, -33120393d0,-394638739d0,  46282816d0, 352046688d0,   !   1870  1.00000
     . 86280520d0, -32888785d0,-394881511d0,  45960047d0, 352306422d0,   !   1871  1.00000
     . 86309284d0, -32657106d0,-395122592d0,  45637161d0, 352564349d0,   !   1872  1.00000
     . 86337845d0, -32425358d0,-395361982d0,  45314159d0, 352820470d0,   !   1873  1.00000
     . 86366203d0, -32193541d0,-395599678d0,  44991041d0, 353074784d0,   !   1874  1.00000
     . 86394358d0, -31961655d0,-395835682d0,  44667809d0, 353327289d0,   !   1875  1.00000
     . 86422311d0, -31729701d0,-396069991d0,  44344464d0, 353577985d0,   !   1876  1.00000
     . 86450061d0, -31497679d0,-396302606d0,  44021005d0, 353826871d0,   !   1877  1.00000
     . 86477608d0, -31265590d0,-396533524d0,  43697435d0, 354073946d0,   !   1878  1.00000
     . 86504952d0, -31033435d0,-396762747d0,  43373754d0, 354319209d0,   !   1879  1.00000
     . 86532093d0, -30801213d0,-396990272d0,  43049962d0, 354562659d0,   !   1880  1.00000
     . 86559031d0, -30568925d0,-397216099d0,  42726061d0, 354804295d0,   !   1881  1.00000
     . 86585765d0, -30336572d0,-397440227d0,  42402051d0, 355044117d0,   !   1882  1.00000
     . 86612296d0, -30104155d0,-397662656d0,  42077933d0, 355282123d0,   !   1883  1.00000
     . 86638624d0, -29871673d0,-397883385d0,  41753709d0, 355518313d0,   !   1884  1.00000
     . 86664748d0, -29639127d0,-398102414d0,  41429378d0, 355752687d0,   !   1885  1.00000
     . 86690669d0, -29406519d0,-398319740d0,  41104943d0, 355985242d0,   !   1886  1.00000
     . 86716386d0, -29173848d0,-398535364d0,  40780402d0, 356215978d0,   !   1887  1.00000
     . 86741900d0, -28941114d0,-398749286d0,  40455758d0, 356444895d0,   !   1888  1.00000
     . 86767210d0, -28708319d0,-398961503d0,  40131012d0, 356671992d0,   !   1889  1.00000
     . 86792315d0, -28475462d0,-399172016d0,  39806163d0, 356897267d0,   !   1890  1.00000
     . 86817218d0, -28242545d0,-399380824d0,  39481213d0, 357120720d0,   !   1891  1.00000
     . 86841916d0, -28009568d0,-399587927d0,  39156164d0, 357342351d0,   !   1892  1.00000
     . 86866410d0, -27776531d0,-399793323d0,  38831014d0, 357562158d0,   !   1893  1.00000
     . 86890700d0, -27543435d0,-399997011d0,  38505766d0, 357780141d0,   !   1894  1.00000
     . 86914786d0, -27310280d0,-400198992d0,  38180421d0, 357996298d0,   !   1895  1.00000
     . 86938667d0, -27077067d0,-400399265d0,  37854978d0, 358210630d0,   !   1896  1.00000
     . 86962345d0, -26843797d0,-400597829d0,  37529439d0, 358423134d0,   !   1897  1.00000
     . 86985818d0, -26610469d0,-400794682d0,  37203806d0, 358633812d0,   !   1898  1.00000
     . 87009087d0, -26377085d0,-400989826d0,  36878077d0, 358842661d0,   !   1899  1.00000
     . 87032151d0, -26143644d0,-401183258d0,  36552255d0, 359049681d0,   !   1900  1.00000
     . 87055010d0, -25910148d0,-401374979d0,  36226341d0, 359254872d0,   !   1901  1.00000
     . 87077665d0, -25676596d0,-401564987d0,  35900334d0, 359458232d0,   !   1902  1.00000
     . 87100116d0, -25442990d0,-401753283d0,  35574237d0, 359659761d0,   !   1903  1.00000
     . 87122362d0, -25209330d0,-401939865d0,  35248049d0, 359859459d0,   !   1904  1.00000
     . 87144402d0, -24975616d0,-402124733d0,  34921772d0, 360057324d0,   !   1905  1.00000
     . 87166238d0, -24741849d0,-402307886d0,  34595406d0, 360253355d0,   !   1906  1.00000
     . 87187870d0, -24508029d0,-402489324d0,  34268952d0, 360447553d0,   !   1907  1.00000
     . 87209296d0, -24274157d0,-402669046d0,  33942412d0, 360639916d0,   !   1908  1.00000
     . 87230517d0, -24040233d0,-402847052d0,  33615786d0, 360830443d0,   !   1909  1.00000
     . 87251533d0, -23806258d0,-403023340d0,  33289074d0, 361019135d0,   !   1910  1.00000
     . 87272344d0, -23572233d0,-403197911d0,  32962278d0, 361205990d0,   !   1911  1.00000
     . 87292950d0, -23338157d0,-403370763d0,  32635398d0, 361391007d0,   !   1912  1.00000
     . 87313350d0, -23104032d0,-403541897d0,  32308436d0, 361574187d0,   !   1913  1.00000
     . 87333545d0, -22869858d0,-403711312d0,  31981392d0, 361755528d0,   !   1914  1.00000
     . 87353535d0, -22635635d0,-403879006d0,  31654267d0, 361935030d0,   !   1915  1.00000
     . 87373320d0, -22401364d0,-404044980d0,  31327061d0, 362112691d0,   !   1916  1.00000
     . 87392899d0, -22167045d0,-404209233d0,  30999777d0, 362288513d0,   !   1917  1.00000
     . 87412272d0, -21932679d0,-404371765d0,  30672413d0, 362462493d0,   !   1918  1.00000
     . 87431440d0, -21698266d0,-404532574d0,  30344972d0, 362634631d0,   !   1919  1.00000
     . 87450402d0, -21463807d0,-404691661d0,  30017454d0, 362804927d0,   !   1920  1.00000
     . 87469159d0, -21229303d0,-404849025d0,  29689861d0, 362973379d0,   !   1921  1.00000
     . 87487709d0, -20994754d0,-405004665d0,  29362192d0, 363139989d0,   !   1922  1.00000
     . 87506054d0, -20760159d0,-405158581d0,  29034448d0, 363304754d0,   !   1923  1.00000
     . 87524194d0, -20525521d0,-405310773d0,  28706631d0, 363467674d0,   !   1924  1.00000
     . 87542127d0, -20290839d0,-405461239d0,  28378742d0, 363628749d0,   !   1925  1.00000
     . 87559854d0, -20056114d0,-405609980d0,  28050780d0, 363787978d0,   !   1926  1.00000
     . 87577376d0, -19821346d0,-405756995d0,  27722748d0, 363945360d0,   !   1927  1.00000
     . 87594691d0, -19586537d0,-405902283d0,  27394645d0, 364100895d0,   !   1928  1.00000
     . 87611800d0, -19351685d0,-406045844d0,  27066473d0, 364254583d0,   !   1929  1.00000
     . 87628704d0, -19116793d0,-406187677d0,  26738232d0, 364406423d0,   !   1930  1.00000
     . 87645401d0, -18881859d0,-406327783d0,  26409924d0, 364556413d0,   !   1931  1.00000
     . 87661891d0, -18646886d0,-406466160d0,  26081549d0, 364704555d0,   !   1932  1.00000
     . 87678176d0, -18411873d0,-406602808d0,  25753108d0, 364850847d0,   !   1933  1.00000
     . 87694254d0, -18176821d0,-406737727d0,  25424602d0, 364995289d0,   !   1934  1.00000
     . 87710126d0, -17941731d0,-406870916d0,  25096031d0, 365137879d0,   !   1935  1.00000
     . 87725791d0, -17706602d0,-407002375d0,  24767397d0, 365278619d0,   !   1936  1.00000
     . 87741251d0, -17471435d0,-407132104d0,  24438701d0, 365417507d0,   !   1937  1.00000
     . 87756503d0, -17236232d0,-407260101d0,  24109942d0, 365554542d0,   !   1938  1.00000
     . 87771549d0, -17000992d0,-407386367d0,  23781123d0, 365689725d0,   !   1939  1.00000
     . 87786389d0, -16765716d0,-407510901d0,  23452243d0, 365823054d0,   !   1940  1.00000
     . 87801022d0, -16530404d0,-407633703d0,  23123304d0, 365954530d0,   !   1941  1.00000
     . 87815448d0, -16295057d0,-407754772d0,  22794307d0, 366084152d0,   !   1942  1.00000
     . 87829667d0, -16059675d0,-407874108d0,  22465252d0, 366211918d0,   !   1943  1.00000
     . 87843680d0, -15824260d0,-407991710d0,  22136140d0, 366337830d0,   !   1944  1.00000
     . 87857486d0, -15588810d0,-408107579d0,  21806972d0, 366461886d0,   !   1945  1.00000
     . 87871086d0, -15353328d0,-408221713d0,  21477749d0, 366584087d0,   !   1946  1.00000
     . 87884478d0, -15117813d0,-408334113d0,  21148471d0, 366704431d0,   !   1947  1.00000
     . 87897664d0, -14882266d0,-408444778d0,  20819141d0, 366822918d0,   !   1948  1.00000
     . 87910643d0, -14646687d0,-408553708d0,  20489757d0, 366939547d0,   !   1949  1.00000
     . 87923414d0, -14411078d0,-408660902d0,  20160322d0, 367054319d0,   !   1950  1.00000
     . 87935979d0, -14175437d0,-408766360d0,  19830836d0, 367167233d0,   !   1951  1.00000
     . 87948337d0, -13939767d0,-408870081d0,  19501299d0, 367278289d0,   !   1952  1.00000
     . 87960488d0, -13704067d0,-408972066d0,  19171714d0, 367387485d0,   !   1953  1.00000
     . 87972431d0, -13468338d0,-409072314d0,  18842080d0, 367494823d0,   !   1954  1.00000
     . 87984168d0, -13232580d0,-409170824d0,  18512398d0, 367600301d0,   !   1955  1.00000
     . 87995697d0, -12996794d0,-409267597d0,  18182670d0, 367703918d0,   !   1956  1.00000
     . 88007020d0, -12760981d0,-409362632d0,  17852895d0, 367805676d0,   !   1957  1.00000
     . 88018135d0, -12525140d0,-409455928d0,  17523076d0, 367905572d0,   !   1958  1.00000
     . 88029042d0, -12289273d0,-409547486d0,  17193212d0, 368003608d0,   !   1959  1.00000
     . 88039743d0, -12053380d0,-409637305d0,  16863305d0, 368099782d0,   !   1960  1.00000
     . 88050236d0, -11817462d0,-409725385d0,  16533355d0, 368194094d0,   !   1961  1.00000
     . 88060522d0, -11581518d0,-409811725d0,  16203363d0, 368286544d0,   !   1962  1.00000
     . 88070601d0, -11345549d0,-409896326d0,  15873331d0, 368377132d0,   !   1963  1.00000
     . 88080472d0, -11109557d0,-409979186d0,  15543258d0, 368465857d0,   !   1964  1.00000
     . 88090136d0, -10873541d0,-410060307d0,  15213146d0, 368552719d0,   !   1965  1.00000
     . 88099592d0, -10637502d0,-410139686d0,  14882996d0, 368637717d0,   !   1966  1.00000
     . 88108841d0, -10401440d0,-410217325d0,  14552808d0, 368720852d0,   !   1967  1.00000
     . 88117883d0, -10165357d0,-410293223d0,  14222583d0, 368802123d0,   !   1968  1.00000
     . 88126717d0,  -9929252d0,-410367379d0,  13892322d0, 368881529d0,   !   1969  1.00000
     . 88135343d0,  -9693125d0,-410439794d0,  13562026d0, 368959071d0,   !   1970  1.00000
     . 88143762d0,  -9456978d0,-410510467d0,  13231696d0, 369034748d0,   !   1971  1.00000
     . 88151973d0,  -9220811d0,-410579398d0,  12901333d0, 369108560d0,   !   1972  1.00000
     . 88159977d0,  -8984625d0,-410646587d0,  12570937d0, 369180507d0,   !   1973  1.00000
     . 88167773d0,  -8748419d0,-410712033d0,  12240509d0, 369250588d0,   !   1974  1.00000
     . 88175362d0,  -8512195d0,-410775737d0,  11910051d0, 369318803d0,   !   1975  1.00000
     . 88182742d0,  -8275953d0,-410837698d0,  11579562d0, 369385153d0,   !   1976  1.00000
     . 88189915d0,  -8039693d0,-410897915d0,  11249044d0, 369449635d0,   !   1977  1.00000
     . 88196881d0,  -7803417d0,-410956390d0,  10918498d0, 369512252d0,   !   1978  1.00000
     . 88203639d0,  -7567123d0,-411013121d0,  10587924d0, 369573001d0,   !   1979  1.00000
     . 88210189d0,  -7330814d0,-411068108d0,  10257323d0, 369631884d0,   !   1980  1.00000
     . 88216531d0,  -7094489d0,-411121351d0,   9926696d0, 369688899d0,   !   1981  1.00000
     . 88222665d0,  -6858149d0,-411172851d0,   9596045d0, 369744047d0,   !   1982  1.00000
     . 88228592d0,  -6621794d0,-411222606d0,   9265369d0, 369797327d0,   !   1983  1.00000
     . 88234311d0,  -6385426d0,-411270617d0,   8934669d0, 369848740d0,   !   1984  1.00000
     . 88239822d0,  -6149044d0,-411316883d0,   8603948d0, 369898284d0,   !   1985  1.00000
     . 88245125d0,  -5912649d0,-411361405d0,   8273204d0, 369945961d0,   !   1986  1.00000
     . 88250220d0,  -5676241d0,-411404182d0,   7942440d0, 369991769d0,   !   1987  1.00000
     . 88255108d0,  -5439821d0,-411445214d0,   7611655d0, 370035708d0,   !   1988  1.00000
     . 88259787d0,  -5203390d0,-411484500d0,   7280851d0, 370077779d0,   !   1989  1.00000
     . 88264259d0,  -4966948d0,-411522042d0,   6950029d0, 370117981d0,   !   1990  1.00000
     . 88268523d0,  -4730495d0,-411557839d0,   6619190d0, 370156315d0,   !   1991  1.00000
     . 88272578d0,  -4494032d0,-411591890d0,   6288333d0, 370192779d0,   !   1992  1.00000
     . 88276426d0,  -4257560d0,-411624195d0,   5957461d0, 370227374d0,   !   1993  1.00000
     . 88280066d0,  -4021078d0,-411654755d0,   5626574d0, 370260100d0,   !   1994  1.00000
     . 88283498d0,  -3784588d0,-411683569d0,   5295673d0, 370290956d0,   !   1995  1.00000
     . 88286722d0,  -3548091d0,-411710637d0,   4964758d0, 370319943d0,   !   1996  1.00000
     . 88289739d0,  -3311585d0,-411735959d0,   4633830d0, 370347060d0,   !   1997  1.00000
     . 88292547d0,  -3075073d0,-411759536d0,   4302891d0, 370372307d0,   !   1998  1.00000
     . 88295147d0,  -2838554d0,-411781366d0,   3971941d0, 370395685d0,   !   1999  1.00000
     . 88297539d0,  -2602029d0,-411801450d0,   3640981d0, 370417192d0,   !   2000  1.00000
     . 88299723d0,  -2365498d0,-411819788d0,   3310012d0, 370436830d0,   !   2001  1.00000
     . 88301699d0,  -2128962d0,-411836379d0,   2979035d0, 370454598d0,   !   2002  1.00000
     . 88303468d0,  -1892422d0,-411851224d0,   2648050d0, 370470495d0,   !   2003  1.00000
     . 88305028d0,  -1655878d0,-411864323d0,   2317058d0, 370484523d0,   !   2004  1.00000
     . 88306380d0,  -1419331d0,-411875676d0,   1986061d0, 370496680d0,   !   2005  1.00000
     . 88307524d0,  -1182780d0,-411885282d0,   1655058d0, 370506967d0,   !   2006  1.00000
     . 88308460d0,   -946227d0,-411893141d0,   1324051d0, 370515383d0,   !   2007  1.00000
     . 88309188d0,   -709672d0,-411899254d0,    993041d0, 370521930d0,   !   2008  1.00000
     . 88309708d0,   -473115d0,-411903620d0,    662029d0, 370526606d0,   !   2009  1.00000
     . 88310020d0,   -236558d0,-411906240d0,    331015d0, 370529411d0,   !   2010  1.00000
     . 88310124d0,         0d0,-411907114d0,         0d0, 370530346d0,   !   2011  1.00000
     . 88310020d0,    236558d0,-411906240d0,   -331015d0, 370529411d0,   !   2012  1.00000
     . 88309708d0,    473115d0,-411903620d0,   -662029d0, 370526606d0,   !   2013  1.00000
     . 88309188d0,    709672d0,-411899254d0,   -993041d0, 370521930d0,   !   2014  1.00000
     . 88308460d0,    946227d0,-411893141d0,  -1324051d0, 370515383d0,   !   2015  1.00000
     . 88307524d0,   1182780d0,-411885282d0,  -1655058d0, 370506967d0,   !   2016  1.00000
     . 88306380d0,   1419331d0,-411875676d0,  -1986061d0, 370496680d0,   !   2017  1.00000
     . 88305028d0,   1655878d0,-411864323d0,  -2317058d0, 370484523d0,   !   2018  1.00000
     . 88303468d0,   1892422d0,-411851224d0,  -2648050d0, 370470495d0,   !   2019  1.00000
     . 88301699d0,   2128962d0,-411836379d0,  -2979035d0, 370454598d0,   !   2020  1.00000
     . 88299723d0,   2365498d0,-411819788d0,  -3310012d0, 370436830d0,   !   2021  1.00000
     . 88297539d0,   2602029d0,-411801450d0,  -3640981d0, 370417192d0,   !   2022  1.00000
     . 88295147d0,   2838554d0,-411781366d0,  -3971941d0, 370395685d0,   !   2023  1.00000
     . 88292547d0,   3075073d0,-411759536d0,  -4302891d0, 370372307d0,   !   2024  1.00000
     . 88289739d0,   3311585d0,-411735959d0,  -4633830d0, 370347060d0,   !   2025  1.00000
     . 88286722d0,   3548091d0,-411710637d0,  -4964758d0, 370319943d0,   !   2026  1.00000
     . 88283498d0,   3784588d0,-411683569d0,  -5295673d0, 370290956d0,   !   2027  1.00000
     . 88280066d0,   4021078d0,-411654755d0,  -5626574d0, 370260100d0,   !   2028  1.00000
     . 88276426d0,   4257560d0,-411624195d0,  -5957461d0, 370227374d0,   !   2029  1.00000
     . 88272578d0,   4494032d0,-411591890d0,  -6288333d0, 370192779d0,   !   2030  1.00000
     . 88268523d0,   4730495d0,-411557839d0,  -6619190d0, 370156315d0,   !   2031  1.00000
     . 88264259d0,   4966948d0,-411522042d0,  -6950029d0, 370117981d0,   !   2032  1.00000
     . 88259787d0,   5203390d0,-411484500d0,  -7280851d0, 370077779d0,   !   2033  1.00000
     . 88255108d0,   5439821d0,-411445214d0,  -7611655d0, 370035708d0,   !   2034  1.00000
     . 88250220d0,   5676241d0,-411404182d0,  -7942440d0, 369991769d0,   !   2035  1.00000
     . 88245125d0,   5912649d0,-411361405d0,  -8273204d0, 369945961d0,   !   2036  1.00000
     . 88239822d0,   6149044d0,-411316883d0,  -8603948d0, 369898284d0,   !   2037  1.00000
     . 88234311d0,   6385426d0,-411270617d0,  -8934669d0, 369848740d0,   !   2038  1.00000
     . 88228592d0,   6621794d0,-411222606d0,  -9265369d0, 369797327d0,   !   2039  1.00000
     . 88222665d0,   6858149d0,-411172851d0,  -9596045d0, 369744047d0,   !   2040  1.00000
     . 88216531d0,   7094489d0,-411121351d0,  -9926696d0, 369688899d0,   !   2041  1.00000
     . 88210189d0,   7330814d0,-411068108d0, -10257323d0, 369631884d0,   !   2042  1.00000
     . 88203639d0,   7567123d0,-411013121d0, -10587924d0, 369573001d0,   !   2043  1.00000
     . 88196881d0,   7803417d0,-410956390d0, -10918498d0, 369512252d0,   !   2044  1.00000
     . 88189915d0,   8039693d0,-410897915d0, -11249044d0, 369449635d0,   !   2045  1.00000
     . 88182742d0,   8275953d0,-410837698d0, -11579562d0, 369385153d0,   !   2046  1.00000
     . 88175362d0,   8512195d0,-410775737d0, -11910051d0, 369318803d0,   !   2047  1.00000
     . 88167773d0,   8748419d0,-410712033d0, -12240509d0, 369250588d0,   !   2048  1.00000
     . 88159977d0,   8984625d0,-410646587d0, -12570937d0, 369180507d0,   !   2049  1.00000
     . 88151973d0,   9220811d0,-410579398d0, -12901333d0, 369108560d0,   !   2050  1.00000
     . 88143762d0,   9456978d0,-410510467d0, -13231696d0, 369034748d0,   !   2051  1.00000
     . 88135343d0,   9693125d0,-410439794d0, -13562026d0, 368959071d0,   !   2052  1.00000
     . 88126717d0,   9929252d0,-410367379d0, -13892322d0, 368881529d0,   !   2053  1.00000
     . 88117883d0,  10165357d0,-410293223d0, -14222583d0, 368802123d0,   !   2054  1.00000
     . 88108841d0,  10401440d0,-410217325d0, -14552808d0, 368720852d0,   !   2055  1.00000
     . 88099592d0,  10637502d0,-410139686d0, -14882996d0, 368637717d0,   !   2056  1.00000
     . 88090136d0,  10873541d0,-410060307d0, -15213146d0, 368552719d0,   !   2057  1.00000
     . 88080472d0,  11109557d0,-409979186d0, -15543258d0, 368465857d0,   !   2058  1.00000
     . 88070601d0,  11345549d0,-409896326d0, -15873331d0, 368377132d0,   !   2059  1.00000
     . 88060522d0,  11581518d0,-409811725d0, -16203363d0, 368286544d0,   !   2060  1.00000
     . 88050236d0,  11817462d0,-409725385d0, -16533355d0, 368194094d0,   !   2061  1.00000
     . 88039743d0,  12053380d0,-409637305d0, -16863305d0, 368099782d0,   !   2062  1.00000
     . 88029042d0,  12289273d0,-409547486d0, -17193212d0, 368003608d0,   !   2063  1.00000
     . 88018135d0,  12525140d0,-409455928d0, -17523076d0, 367905572d0,   !   2064  1.00000
     . 88007020d0,  12760981d0,-409362632d0, -17852895d0, 367805676d0,   !   2065  1.00000
     . 87995697d0,  12996794d0,-409267597d0, -18182670d0, 367703918d0,   !   2066  1.00000
     . 87984168d0,  13232580d0,-409170824d0, -18512398d0, 367600301d0,   !   2067  1.00000
     . 87972431d0,  13468338d0,-409072314d0, -18842080d0, 367494823d0,   !   2068  1.00000
     . 87960488d0,  13704067d0,-408972066d0, -19171714d0, 367387485d0,   !   2069  1.00000
     . 87948337d0,  13939767d0,-408870081d0, -19501299d0, 367278289d0,   !   2070  1.00000
     . 87935979d0,  14175437d0,-408766360d0, -19830836d0, 367167233d0,   !   2071  1.00000
     . 87923414d0,  14411078d0,-408660902d0, -20160322d0, 367054319d0,   !   2072  1.00000
     . 87910643d0,  14646687d0,-408553708d0, -20489757d0, 366939547d0,   !   2073  1.00000
     . 87897664d0,  14882266d0,-408444778d0, -20819141d0, 366822918d0,   !   2074  1.00000
     . 87884478d0,  15117813d0,-408334113d0, -21148471d0, 366704431d0,   !   2075  1.00000
     . 87871086d0,  15353328d0,-408221713d0, -21477749d0, 366584087d0,   !   2076  1.00000
     . 87857486d0,  15588810d0,-408107579d0, -21806972d0, 366461886d0,   !   2077  1.00000
     . 87843680d0,  15824260d0,-407991710d0, -22136140d0, 366337830d0,   !   2078  1.00000
     . 87829667d0,  16059675d0,-407874108d0, -22465252d0, 366211918d0,   !   2079  1.00000
     . 87815448d0,  16295057d0,-407754772d0, -22794307d0, 366084152d0,   !   2080  1.00000
     . 87801022d0,  16530404d0,-407633703d0, -23123304d0, 365954530d0,   !   2081  1.00000
     . 87786389d0,  16765716d0,-407510901d0, -23452243d0, 365823054d0,   !   2082  1.00000
     . 87771549d0,  17000992d0,-407386367d0, -23781123d0, 365689725d0,   !   2083  1.00000
     . 87756503d0,  17236232d0,-407260101d0, -24109942d0, 365554542d0,   !   2084  1.00000
     . 87741251d0,  17471435d0,-407132104d0, -24438701d0, 365417507d0,   !   2085  1.00000
     . 87725791d0,  17706602d0,-407002375d0, -24767397d0, 365278619d0,   !   2086  1.00000
     . 87710126d0,  17941731d0,-406870916d0, -25096031d0, 365137879d0,   !   2087  1.00000
     . 87694254d0,  18176821d0,-406737727d0, -25424602d0, 364995289d0,   !   2088  1.00000
     . 87678176d0,  18411873d0,-406602808d0, -25753108d0, 364850847d0,   !   2089  1.00000
     . 87661891d0,  18646886d0,-406466160d0, -26081549d0, 364704555d0,   !   2090  1.00000
     . 87645401d0,  18881859d0,-406327783d0, -26409924d0, 364556413d0,   !   2091  1.00000
     . 87628704d0,  19116793d0,-406187677d0, -26738232d0, 364406423d0,   !   2092  1.00000
     . 87611800d0,  19351685d0,-406045844d0, -27066473d0, 364254583d0,   !   2093  1.00000
     . 87594691d0,  19586537d0,-405902283d0, -27394645d0, 364100895d0,   !   2094  1.00000
     . 87577376d0,  19821346d0,-405756995d0, -27722748d0, 363945360d0,   !   2095  1.00000
     . 87559854d0,  20056114d0,-405609980d0, -28050780d0, 363787978d0,   !   2096  1.00000
     . 87542127d0,  20290839d0,-405461239d0, -28378742d0, 363628749d0,   !   2097  1.00000
     . 87524194d0,  20525521d0,-405310773d0, -28706631d0, 363467674d0,   !   2098  1.00000
     . 87506054d0,  20760159d0,-405158581d0, -29034448d0, 363304754d0,   !   2099  1.00000
     . 87487709d0,  20994754d0,-405004665d0, -29362192d0, 363139989d0,   !   2100  1.00000
     . 87469159d0,  21229303d0,-404849025d0, -29689861d0, 362973379d0,   !   2101  1.00000
     . 87450402d0,  21463807d0,-404691661d0, -30017454d0, 362804927d0,   !   2102  1.00000
     . 87431440d0,  21698266d0,-404532574d0, -30344972d0, 362634631d0,   !   2103  1.00000
     . 87412272d0,  21932679d0,-404371765d0, -30672413d0, 362462493d0,   !   2104  1.00000
     . 87392899d0,  22167045d0,-404209233d0, -30999777d0, 362288513d0,   !   2105  1.00000
     . 87373320d0,  22401364d0,-404044980d0, -31327061d0, 362112691d0,   !   2106  1.00000
     . 87353535d0,  22635635d0,-403879006d0, -31654267d0, 361935030d0,   !   2107  1.00000
     . 87333545d0,  22869858d0,-403711312d0, -31981392d0, 361755528d0,   !   2108  1.00000
     . 87313350d0,  23104032d0,-403541897d0, -32308436d0, 361574187d0,   !   2109  1.00000
     . 87292950d0,  23338157d0,-403370763d0, -32635398d0, 361391007d0,   !   2110  1.00000
     . 87272344d0,  23572233d0,-403197911d0, -32962278d0, 361205990d0,   !   2111  1.00000
     . 87251533d0,  23806258d0,-403023340d0, -33289074d0, 361019135d0,   !   2112  1.00000
     . 87230517d0,  24040233d0,-402847052d0, -33615786d0, 360830443d0,   !   2113  1.00000
     . 87209296d0,  24274157d0,-402669046d0, -33942412d0, 360639916d0,   !   2114  1.00000
     . 87187870d0,  24508029d0,-402489324d0, -34268952d0, 360447553d0,   !   2115  1.00000
     . 87166238d0,  24741849d0,-402307886d0, -34595406d0, 360253355d0,   !   2116  1.00000
     . 87144402d0,  24975616d0,-402124733d0, -34921772d0, 360057324d0,   !   2117  1.00000
     . 87122362d0,  25209330d0,-401939865d0, -35248049d0, 359859459d0,   !   2118  1.00000
     . 87100116d0,  25442990d0,-401753283d0, -35574237d0, 359659761d0,   !   2119  1.00000
     . 87077665d0,  25676596d0,-401564987d0, -35900334d0, 359458232d0,   !   2120  1.00000
     . 87055010d0,  25910148d0,-401374979d0, -36226341d0, 359254872d0,   !   2121  1.00000
     . 87032151d0,  26143644d0,-401183258d0, -36552255d0, 359049681d0,   !   2122  1.00000
     . 87009087d0,  26377085d0,-400989826d0, -36878077d0, 358842661d0,   !   2123  1.00000
     . 86985818d0,  26610469d0,-400794682d0, -37203806d0, 358633812d0,   !   2124  1.00000
     . 86962345d0,  26843797d0,-400597829d0, -37529439d0, 358423134d0,   !   2125  1.00000
     . 86938667d0,  27077067d0,-400399265d0, -37854978d0, 358210630d0,   !   2126  1.00000
     . 86914786d0,  27310280d0,-400198992d0, -38180421d0, 357996298d0,   !   2127  1.00000
     . 86890700d0,  27543435d0,-399997011d0, -38505766d0, 357780141d0,   !   2128  1.00000
     . 86866410d0,  27776531d0,-399793323d0, -38831014d0, 357562158d0,   !   2129  1.00000
     . 86841916d0,  28009568d0,-399587927d0, -39156164d0, 357342351d0,   !   2130  1.00000
     . 86817218d0,  28242545d0,-399380824d0, -39481213d0, 357120720d0,   !   2131  1.00000
     . 86792315d0,  28475462d0,-399172016d0, -39806163d0, 356897267d0,   !   2132  1.00000
     . 86767210d0,  28708319d0,-398961503d0, -40131012d0, 356671992d0,   !   2133  1.00000
     . 86741900d0,  28941114d0,-398749286d0, -40455758d0, 356444895d0,   !   2134  1.00000
     . 86716386d0,  29173848d0,-398535364d0, -40780402d0, 356215978d0,   !   2135  1.00000
     . 86690669d0,  29406519d0,-398319740d0, -41104943d0, 355985242d0,   !   2136  1.00000
     . 86664748d0,  29639127d0,-398102414d0, -41429378d0, 355752687d0,   !   2137  1.00000
     . 86638624d0,  29871673d0,-397883385d0, -41753709d0, 355518313d0,   !   2138  1.00000
     . 86612296d0,  30104155d0,-397662656d0, -42077933d0, 355282123d0,   !   2139  1.00000
     . 86585765d0,  30336572d0,-397440227d0, -42402051d0, 355044117d0,   !   2140  1.00000
     . 86559031d0,  30568925d0,-397216099d0, -42726061d0, 354804295d0,   !   2141  1.00000
     . 86532093d0,  30801213d0,-396990272d0, -43049962d0, 354562659d0,   !   2142  1.00000
     . 86504952d0,  31033435d0,-396762747d0, -43373754d0, 354319209d0,   !   2143  1.00000
     . 86477608d0,  31265590d0,-396533524d0, -43697435d0, 354073946d0,   !   2144  1.00000
     . 86450061d0,  31497679d0,-396302606d0, -44021005d0, 353826871d0,   !   2145  1.00000
     . 86422311d0,  31729701d0,-396069991d0, -44344464d0, 353577985d0,   !   2146  1.00000
     . 86394358d0,  31961655d0,-395835682d0, -44667809d0, 353327289d0,   !   2147  1.00000
     . 86366203d0,  32193541d0,-395599678d0, -44991041d0, 353074784d0,   !   2148  1.00000
     . 86337845d0,  32425358d0,-395361982d0, -45314159d0, 352820470d0,   !   2149  1.00000
     . 86309284d0,  32657106d0,-395122592d0, -45637161d0, 352564349d0,   !   2150  1.00000
     . 86280520d0,  32888785d0,-394881511d0, -45960047d0, 352306422d0,   !   2151  1.00000
     . 86251554d0,  33120393d0,-394638739d0, -46282816d0, 352046688d0,   !   2152  1.00000
     . 86222386d0,  33351930d0,-394394276d0, -46605468d0, 351785151d0,   !   2153  1.00000
     . 86193016d0,  33583396d0,-394148124d0, -46928000d0, 351521809d0,   !   2154  1.00000
     . 86163443d0,  33814791d0,-393900284d0, -47250413d0, 351256665d0,   !   2155  1.00000
     . 86133668d0,  34046113d0,-393650755d0, -47572706d0, 350989718d0,   !   2156  1.00000
     . 86103691d0,  34277362d0,-393399540d0, -47894878d0, 350720971d0,   !   2157  1.00000
     . 86073512d0,  34508539d0,-393146639d0, -48216928d0, 350450424d0,   !   2158  1.00000
     . 86043131d0,  34739641d0,-392892052d0, -48538855d0, 350178078d0,   !   2159  1.00000
     . 86012548d0,  34970669d0,-392635781d0, -48860658d0, 349903934d0,   !   2160  1.00000
     . 85981763d0,  35201623d0,-392377826d0, -49182337d0, 349627993d0,   !   2161  1.00000
     . 85950777d0,  35432501d0,-392118188d0, -49503891d0, 349350256d0,   !   2162  1.00000
     . 85919590d0,  35663304d0,-391856868d0, -49825318d0, 349070723d0,   !   2163  1.00000
     . 85888200d0,  35894030d0,-391593867d0, -50146618d0, 348789397d0,   !   2164  1.00000
     . 85856610d0,  36124680d0,-391329186d0, -50467791d0, 348506278d0,   !   2165  1.00000
     . 85824818d0,  36355252d0,-391062826d0, -50788835d0, 348221367d0,   !   2166  1.00000
     . 85792825d0,  36585747d0,-390794787d0, -51109749d0, 347934665d0,   !   2167  1.00000
     . 85760631d0,  36816164d0,-390525071d0, -51430533d0, 347646172d0,   !   2168  1.00000
     . 85728235d0,  37046501d0,-390253678d0, -51751186d0, 347355891d0,   !   2169  1.00000
     . 85695639d0,  37276760d0,-389980609d0, -52071707d0, 347063822d0,   !   2170  1.00000
     . 85662842d0,  37506938d0,-389705865d0, -52392095d0, 346769967d0,   !   2171  1.00000
     . 85629844d0,  37737037d0,-389429447d0, -52712350d0, 346474325d0,   !   2172  1.00000
     . 85596646d0,  37967055d0,-389151356d0, -53032470d0, 346176899d0,   !   2173  1.00000
     . 85563247d0,  38196991d0,-388871593d0, -53352454d0, 345877689d0,   !   2174  1.00000
     . 85529647d0,  38426846d0,-388590158d0, -53672302d0, 345576696d0,   !   2175  1.00000
     . 85495847d0,  38656618d0,-388307053d0, -53992014d0, 345273922d0,   !   2176  1.00000
     . 85461847d0,  38886308d0,-388022279d0, -54311587d0, 344969368d0,   !   2177  1.00000
     . 85427646d0,  39115915d0,-387735836d0, -54631022d0, 344663034d0,   !   2178  1.00000
     . 85393246d0,  39345437d0,-387447726d0, -54950317d0, 344354922d0,   !   2179  1.00000
     . 85358645d0,  39574876d0,-387157949d0, -55269472d0, 344045033d0,   !   2180  1.00000
     . 85323844d0,  39804229d0,-386866507d0, -55588486d0, 343733367d0,   !   2181  1.00000
     . 85288844d0,  40033498d0,-386573400d0, -55907358d0, 343419927d0,   !   2182  1.00000
     . 85253644d0,  40262681d0,-386278629d0, -56226087d0, 343104714d0,   !   2183  1.00000
     . 85218244d0,  40491777d0,-385982196d0, -56544672d0, 342787727d0,   !   2184  1.00000
     . 85182644d0,  40720787d0,-385684100d0, -56863112d0, 342468969d0,   !   2185  1.00000
     . 85146846d0,  40949709d0,-385384344d0, -57181408d0, 342148441d0,   !   2186  1.00000
     . 85110847d0,  41178544d0,-385082928d0, -57499557d0, 341826144d0,   !   2187  1.00000
     . 85074650d0,  41407290d0,-384779854d0, -57817559d0, 341502079d0,   !   2188  1.00000
     . 85038253d0,  41635948d0,-384475121d0, -58135413d0, 341176246d0,   !   2189  1.00000
     . 85001658d0,  41864516d0,-384168732d0, -58453118d0, 340848648d0,   !   2190  1.00000
     . 84964863d0,  42092995d0,-383860687d0, -58770674d0, 340519286d0,   !   2191  1.00000
     . 84927870d0,  42321383d0,-383550988d0, -59088080d0, 340188160d0,   !   2192  1.00000
     . 84890678d0,  42549680d0,-383239634d0, -59405334d0, 339855272d0,   !   2193  1.00000
     . 84853287d0,  42777887d0,-382926628d0, -59722436d0, 339520624d0,   !   2194  1.00000
     . 84815698d0,  43006001d0,-382611970d0, -60039385d0, 339184215d0,   !   2195  1.00000
     . 84777910d0,  43234023d0,-382295662d0, -60356181d0, 338846048d0,   !   2196  1.00000
     . 84739924d0,  43461953d0,-381977704d0, -60672822d0, 338506124d0,   !   2197  1.00000
     . 84701740d0,  43689789d0,-381658098d0, -60989308d0, 338164444d0,   !   2198  1.00000
     . 84663357d0,  43917531d0,-381336844d0, -61305637d0, 337821008d0,   !   2199  1.00000
     . 84624777d0,  44145179d0,-381013943d0, -61621810d0, 337475820d0,   !   2200  1.00000
     . 84585999d0,  44372732d0,-380689398d0, -61937824d0, 337128879d0,   !   2201  1.00000
     . 84547023d0,  44600190d0,-380363208d0, -62253680d0, 336780186d0,   !   2202  1.00000
     . 84507849d0,  44827552d0,-380035375d0, -62569376d0, 336429744d0,   !   2203  1.00000
     . 84468477d0,  45054818d0,-379705899d0, -62884912d0, 336077554d0,   !   2204  1.00000
     . 84428909d0,  45281987d0,-379374783d0, -63200286d0, 335723616d0,   !   2205  1.00000
     . 84389142d0,  45509059d0,-379042027d0, -63515499d0, 335367932d0,   !   2206  1.00000
     . 84349179d0,  45736033d0,-378707632d0, -63830548d0, 335010504d0,   !   2207  1.00000
     . 84309018d0,  45962909d0,-378371599d0, -64145434d0, 334651332d0,   !   2208  1.00000
     . 84268661d0,  46189686d0,-378033930d0, -64460155d0, 334290417d0,   !   2209  1.00000
     . 84228106d0,  46416363d0,-377694626d0, -64774711d0, 333927763d0,   !   2210  1.00000
     . 84187355d0,  46642941d0,-377353687d0, -65089100d0, 333563368d0,   !   2211  1.00000
     . 84146407d0,  46869418d0,-377011114d0, -65403323d0, 333197235d0,   !   2212  1.00000
     . 84105263d0,  47095795d0,-376666910d0, -65717377d0, 332829366d0,   !   2213  1.00000
     . 84063922d0,  47322070d0,-376321075d0, -66031263d0, 332459761d0,   !   2214  1.00000
     . 84022384d0,  47548244d0,-375973610d0, -66344979d0, 332088421d0,   !   2215  1.00000
     . 83980651d0,  47774315d0,-375624516d0, -66658525d0, 331715349d0,   !   2216  1.00000
     . 83938721d0,  48000284d0,-375273795d0, -66971899d0, 331340546d0,   !   2217  1.00000
     . 83896595d0,  48226149d0,-374921447d0, -67285101d0, 330964012d0,   !   2218  1.00000
     . 83854274d0,  48451910d0,-374567475d0, -67598131d0, 330585749d0,   !   2219  1.00000
     . 83811756d0,  48677567d0,-374211878d0, -67910986d0, 330205759d0,   !   2220  1.00000
     . 83769043d0,  48903119d0,-373854659d0, -68223667d0, 329824043d0,   !   2221  1.00000
     . 83726135d0,  49128566d0,-373495818d0, -68536173d0, 329440602d0,   !   2222  1.00000
     . 83683031d0,  49353907d0,-373135357d0, -68848502d0, 329055438d0,   !   2223  1.00000
     . 83639732d0,  49579142d0,-372773276d0, -69160654d0, 328668553d0,   !   2224  1.00000
     . 83596237d0,  49804269d0,-372409578d0, -69472628d0, 328279946d0,   !   2225  1.00000
     . 83552548d0,  50029290d0,-372044263d0, -69784424d0, 327889621d0,   !   2226  1.00000
     . 83508663d0,  50254203d0,-371677332d0, -70096040d0, 327497578d0,   !   2227  1.00000
     . 83464584d0,  50479007d0,-371308787d0, -70407475d0, 327103819d0,   !   2228  1.00000
     . 83420311d0,  50703703d0,-370938628d0, -70718729d0, 326708345d0,   !   2229  1.00000
     . 83375842d0,  50928289d0,-370566858d0, -71029801d0, 326311158d0,   !   2230  1.00000
     . 83331179d0,  51152765d0,-370193477d0, -71340690d0, 325912258d0,   !   2231  1.00000
     . 83286322d0,  51377131d0,-369818487d0, -71651395d0, 325511649d0,   !   2232  1.00000
     . 83241271d0,  51601386d0,-369441889d0, -71961915d0, 325109330d0,   !   2233  1.00000
     . 83196026d0,  51825530d0,-369063683d0, -72272250d0, 324705304d0,   !   2234  1.00000
     . 83150586d0,  52049562d0,-368683872d0, -72582399d0, 324299572d0,   !   2235  1.00000
     . 83104953d0,  52273482d0,-368302457d0, -72892361d0, 323892135d0,   !   2236  1.00000
     . 83059127d0,  52497289d0,-367919438d0, -73202134d0, 323482995d0,   !   2237  1.00000
     . 83013106d0,  52720982d0,-367534818d0, -73511719d0, 323072154d0,   !   2238  1.00000
     . 82966893d0,  52944562d0,-367148597d0, -73821114d0, 322659612d0,   !   2239  1.00000
     . 82920486d0,  53168027d0,-366760777d0, -74130319d0, 322245372d0,   !   2240  1.00000
     . 82873885d0,  53391378d0,-366371359d0, -74439332d0, 321829435d0,   !   2241  1.00000
     . 82827092d0,  53614613d0,-365980345d0, -74748153d0, 321411802d0,   !   2242  1.00000
     . 82780106d0,  53837732d0,-365587735d0, -75056782d0, 320992476d0,   !   2243  1.00000
     . 82732927d0,  54060735d0,-365193531d0, -75365216d0, 320571456d0,   !   2244  1.00000
     . 82685556d0,  54283621d0,-364797734d0, -75673456d0, 320148746d0,   !   2245  1.00000
     . 82637992d0,  54506390d0,-364400346d0, -75981500d0, 319724346d0,   !   2246  1.00000
     . 82590235d0,  54729041d0,-364001368d0, -76289348d0, 319298258d0,   !   2247  1.00000
     . 82542287d0,  54951573d0,-363600801d0, -76596998d0, 318870484d0,   !   2248  1.00000
     . 82494146d0,  55173987d0,-363198646d0, -76904451d0, 318441026d0,   !   2249  1.00000
     . 82445813d0,  55396281d0,-362794906d0, -77211705d0, 318009884d0,   !   2250  1.00000
     . 82397289d0,  55618455d0,-362389581d0, -77518759d0, 317577060d0,   !   2251  1.00000
     . 82348573d0,  55840509d0,-361982673d0, -77825613d0, 317142557d0,   !   2252  1.00000
     . 82299665d0,  56062442d0,-361574182d0, -78132265d0, 316706374d0,   !   2253  1.00000
     . 82250566d0,  56284254d0,-361164111d0, -78438715d0, 316268516d0,   !   2254  1.00000
     . 82201275d0,  56505944d0,-360752461d0, -78744963d0, 315828981d0,   !   2255  1.00000
     . 82151793d0,  56727512d0,-360339233d0, -79051006d0, 315387774d0,   !   2256  1.00000
     . 82102121d0,  56948956d0,-359924429d0, -79356844d0, 314944894d0,   !   2257  1.00000
     . 82052257d0,  57170277d0,-359508049d0, -79662478d0, 314500344d0,   !   2258  1.00000
     . 82002203d0,  57391475d0,-359090096d0, -79967904d0, 314054124d0,   !   2259  1.00000
     . 81951958d0,  57612548d0,-358670570d0, -80273124d0, 313606238d0,   !   2260  1.00000
     . 81901523d0,  57833496d0,-358249474d0, -80578135d0, 313156687d0,   !   2261  1.00000
     . 81850897d0,  58054319d0,-357826808d0, -80882938d0, 312705471d0,   !   2262  1.00000
     . 81800082d0,  58275016d0,-357402574d0, -81187531d0, 312252593d0,   !   2263  1.00000
     . 81749076d0,  58495586d0,-356976774d0, -81491913d0, 311798055d0,   !   2264  1.00000
     . 81697880d0,  58716030d0,-356549408d0, -81796084d0, 311341858d0,   !   2265  1.00000
     . 81646495d0,  58936346d0,-356120479d0, -82100043d0, 310884004d0,   !   2266  1.00000
     . 81594920d0,  59156534d0,-355689987d0, -82403789d0, 310424494d0,   !   2267  1.00000
     . 81543155d0,  59376594d0,-355257935d0, -82707321d0, 309963331d0,   !   2268  1.00000
     . 81491202d0,  59596525d0,-354824323d0, -83010638d0, 309500515d0,   !   2269  1.00000
     . 81439059d0,  59816327d0,-354389153d0, -83313740d0, 309036049d0,   !   2270  1.00000
     . 81386727d0,  60035999d0,-353952427d0, -83616625d0, 308569934d0,   !   2271  1.00000
     . 81334206d0,  60255540d0,-353514145d0, -83919293d0, 308102172d0,   !   2272  1.00000
     . 81281497d0,  60474951d0,-353074311d0, -84221743d0, 307632765d0,   !   2273  1.00000
     . 81228598d0,  60694230d0,-352632924d0, -84523974d0, 307161714d0,   !   2274  1.00000
     . 81175512d0,  60913377d0,-352189986d0, -84825985d0, 306689022d0,   !   2275  1.00000
     . 81122237d0,  61132392d0,-351745500d0, -85127776d0, 306214689d0,   !   2276  1.00000
     . 81068774d0,  61351273d0,-351299466d0, -85429346d0, 305738718d0,   !   2277  1.00000
     . 81015123d0,  61570022d0,-350851886d0, -85730693d0, 305261111d0,   !   2278  1.00000
     . 80961284d0,  61788636d0,-350402761d0, -86031817d0, 304781869d0,   !   2279  1.00000
     . 80907258d0,  62007117d0,-349952094d0, -86332717d0, 304300993d0,   !   2280  1.00000
     . 80853044d0,  62225462d0,-349499884d0, -86633393d0, 303818487d0,   !   2281  1.00000
     . 80798643d0,  62443672d0,-349046135d0, -86933843d0, 303334351d0,   !   2282  1.00000
     . 80744054d0,  62661746d0,-348590848d0, -87234066d0, 302848587d0,   !   2283  1.00000
     . 80689279d0,  62879683d0,-348134023d0, -87534062d0, 302361198d0,   !   2284  1.00000
     . 80634317d0,  63097484d0,-347675663d0, -87833831d0, 301872184d0,   !   2285  1.00000
     . 80579167d0,  63315147d0,-347215770d0, -88133370d0, 301381548d0,   !   2286  1.00000
     . 80523832d0,  63532672d0,-346754344d0, -88432679d0, 300889291d0,   !   2287  1.00000
     . 80468310d0,  63750059d0,-346291387d0, -88731758d0, 300395416d0,   !   2288  1.00000
     . 80412601d0,  63967307d0,-345826902d0, -89030606d0, 299899924d0,   !   2289  1.00000
     . 80356707d0,  64184415d0,-345360888d0, -89329221d0, 299402817d0,   !   2290  1.00000
     . 80300626d0,  64401384d0,-344893349d0, -89627603d0, 298904097d0,   !   2291  1.00000
     . 80244360d0,  64618212d0,-344424286d0, -89925751d0, 298403766d0,   !   2292  1.00000
     . 80187908d0,  64834899d0,-343953699d0, -90223665d0, 297901825d0,   !   2293  1.00000
     . 80131271d0,  65051445d0,-343481592d0, -90521343d0, 297398276d0,   !   2294  1.00000
     . 80074448d0,  65267849d0,-343007964d0, -90818784d0, 296893122d0,   !   2295  1.00000
     . 80017440d0,  65484110d0,-342532819d0, -91115988d0, 296386364d0,   !   2296  1.00000
     . 79960248d0,  65700228d0,-342056158d0, -91412954d0, 295878004d0,   !   2297  1.00000
     . 79902870d0,  65916203d0,-341577981d0, -91709681d0, 295368044d0,   !   2298  1.00000
     . 79845308d0,  66132035d0,-341098292d0, -92006169d0, 294856486d0,   !   2299  1.00000
     . 79787561d0,  66347721d0,-340617091d0, -92302416d0, 294343331d0,   !   2300  1.00000
     . 79729629d0,  66563263d0,-340134380d0, -92598421d0, 293828582d0,   !   2301  1.00000
     . 79671514d0,  66778659d0,-339650161d0, -92894184d0, 293312241d0,   !   2302  1.00000
     . 79613215d0,  66993909d0,-339164435d0, -93189704d0, 292794309d0,   !   2303  1.00000
     . 79554731d0,  67209013d0,-338677204d0, -93484980d0, 292274788d0,   !   2304  1.00000
     . 79496064d0,  67423970d0,-338188470d0, -93780011d0, 291753681d0,   !   2305  1.00000
     . 79437214d0,  67638780d0,-337698235d0, -94074796d0, 291230989d0,   !   2306  1.00000
     . 79378180d0,  67853441d0,-337206499d0, -94369335d0, 290706714d0,   !   2307  1.00000
     . 79318963d0,  68067954d0,-336713265d0, -94663627d0, 290180858d0,   !   2308  1.00000
     . 79259563d0,  68282318d0,-336218535d0, -94957671d0, 289653424d0,   !   2309  1.00000
     . 79199980d0,  68496533d0,-335722310d0, -95251465d0, 289124412d0,   !   2310  1.00000
     . 79140215d0,  68710598d0,-335224591d0, -95545010d0, 288593826d0,   !   2311  1.00000
     . 79080266d0,  68924512d0,-334725381d0, -95838304d0, 288061666d0,   !   2312  1.00000
     . 79020136d0,  69138275d0,-334224681d0, -96131347d0, 287527936d0,   !   2313  1.00000
     . 78959823d0,  69351887d0,-333722493d0, -96424138d0, 286992636d0,   !   2314  1.00000
     . 78899329d0,  69565347d0,-333218819d0, -96716675d0, 286455770d0,   !   2315  1.00000
     . 78838652d0,  69778654d0,-332713660d0, -97008958d0, 285917339d0,   !   2316  1.00000
     . 78777794d0,  69991809d0,-332207018d0, -97300987d0, 285377344d0,   !   2317  1.00000
     . 78716755d0,  70204810d0,-331698895d0, -97592760d0, 284835789d0,   !   2318  1.00000
     . 78655534d0,  70417657d0,-331189292d0, -97884276d0, 284292674d0,   !   2319  1.00000
     . 78594132d0,  70630349d0,-330678212d0, -98175535d0, 283748003d0,   !   2320  1.00000
     . 78532549d0,  70842887d0,-330165655d0, -98466536d0, 283201777d0,   !   2321  1.00000
     . 78470785d0,  71055269d0,-329651625d0, -98757278d0, 282653998d0,   !   2322  1.00000
     . 78408841d0,  71267495d0,-329136122d0, -99047760d0, 282104669d0,   !   2323  1.00000
     . 78346716d0,  71479565d0,-328619148d0, -99337981d0, 281553790d0,   !   2324  1.00000
     . 78284411d0,  71691477d0,-328100705d0, -99627941d0, 281001365d0,   !   2325  1.00000
     . 78221926d0,  71903233d0,-327580795d0, -99917639d0, 280447396d0,   !   2326  1.00000
     . 78159261d0,  72114830d0,-327059420d0,-100207073d0, 279891884d0,   !   2327  1.00000
     . 78096416d0,  72326269d0,-326536581d0,-100496244d0, 279334831d0,   !   2328  1.00000
     . 78033392d0,  72537548d0,-326012281d0,-100785149d0, 278776240d0,   !   2329  1.00000
     . 77970188d0,  72748668d0,-325486520d0,-101073789d0, 278216113d0,   !   2330  1.00000
     . 77906806d0,  72959629d0,-324959301d0,-101362163d0, 277654452d0,   !   2331  1.00000
     . 77843244d0,  73170429d0,-324430626d0,-101650269d0, 277091259d0,   !   2332  1.00000
     . 77779503d0,  73381067d0,-323900496d0,-101938107d0, 276526536d0,   !   2333  1.00000
     . 77715584d0,  73591545d0,-323368914d0,-102225675d0, 275960285d0,   !   2334  1.00000
     . 77651486d0,  73801860d0,-322835880d0,-102512974d0, 275392509d0,   !   2335  1.00000
     . 77587210d0,  74012013d0,-322301398d0,-102800003d0, 274823209d0,   !   2336  1.00000
     . 77522756d0,  74222003d0,-321765468d0,-103086759d0, 274252388d0,   !   2337  1.00000
     . 77458125d0,  74431829d0,-321228093d0,-103373244d0, 273680048d0,   !   2338  1.00000
     . 77393315d0,  74641491d0,-320689274d0,-103659455d0, 273106190d0,   !   2339  1.00000
     . 77328328d0,  74850989d0,-320149014d0,-103945392d0, 272530818d0,   !   2340  1.00000
     . 77263163d0,  75060322d0,-319607314d0,-104231055d0, 271953933d0,   !   2341  1.00000
     . 77197822d0,  75269489d0,-319064176d0,-104516441d0, 271375537d0,   !   2342  1.00000
     . 77132303d0,  75478490d0,-318519601d0,-104801551d0, 270795634d0,   !   2343  1.00000
     . 77066608d0,  75687325d0,-317973593d0,-105086384d0, 270214224d0,   !   2344  1.00000
     . 77000736d0,  75895993d0,-317426152d0,-105370939d0, 269631310d0,   !   2345  1.00000
     . 76934687d0,  76104493d0,-316877281d0,-105655214d0, 269046894d0,   !   2346  1.00000
     . 76868463d0,  76312826d0,-316326981d0,-105939210d0, 268460979d0,   !   2347  1.00000
     . 76802062d0,  76520989d0,-315775255d0,-106222925d0, 267873566d0,   !   2348  1.00000
     . 76735486d0,  76728984d0,-315222104d0,-106506358d0, 267284659d0,   !   2349  1.00000
     . 76668734d0,  76936809d0,-314667530d0,-106789509d0, 266694258d0,   !   2350  1.00000
     . 76601806d0,  77144464d0,-314111535d0,-107072377d0, 266102367d0,   !   2351  1.00000
     . 76534703d0,  77351949d0,-313554121d0,-107354960d0, 265508988d0,   !   2352  1.00000
     . 76467426d0,  77559263d0,-312995291d0,-107637259d0, 264914122d0,   !   2353  1.00000
     . 76399973d0,  77766405d0,-312435045d0,-107919272d0, 264317772d0,   !   2354  1.00000
     . 76332345d0,  77973375d0,-311873386d0,-108200999d0, 263719941d0,   !   2355  1.00000
     . 76264544d0,  78180173d0,-311310315d0,-108482438d0, 263120631d0,   !   2356  1.00000
     . 76196568d0,  78386797d0,-310745836d0,-108763589d0, 262519843d0,   !   2357  1.00000
     . 76128417d0,  78593248d0,-310179949d0,-109044450d0, 261917581d0,   !   2358  1.00000
     . 76060093d0,  78799525d0,-309612657d0,-109325022d0, 261313846d0,   !   2359  1.00000
     . 75991595d0,  79005628d0,-309043961d0,-109605304d0, 260708641d0,   !   2360  1.00000
     . 75922924d0,  79211555d0,-308473864d0,-109885293d0, 260101969d0,   !   2361  1.00000
     . 75854080d0,  79417307d0,-307902368d0,-110164990d0, 259493830d0,   !   2362  1.00000
     . 75785062d0,  79622883d0,-307329474d0,-110444394d0, 258884228d0,   !   2363  1.00000
     . 75715871d0,  79828282d0,-306755185d0,-110723504d0, 258273166d0,   !   2364  1.00000
     . 75646508d0,  80033505d0,-306179502d0,-111002319d0, 257660645d0,   !   2365  1.00000
     . 75576972d0,  80238550d0,-305602428d0,-111280838d0, 257046667d0,   !   2366  1.00000
     . 75507264d0,  80443416d0,-305023964d0,-111559061d0, 256431236d0,   !   2367  1.00000
     . 75437384d0,  80648104d0,-304444113d0,-111836986d0, 255814353d0,   !   2368  1.00000
     . 75367332d0,  80852614d0,-303862877d0,-112114613d0, 255196020d0,   !   2369  1.00000
     . 75297108d0,  81056943d0,-303280257d0,-112391941d0, 254576241d0,   !   2370  1.00000
     . 75226713d0,  81261093d0,-302696256d0,-112668969d0, 253955018d0,   !   2371  1.00000
     . 75156146d0,  81465062d0,-302110876d0,-112945696d0, 253332352d0,   !   2372  1.00000
     . 75085409d0,  81668850d0,-301524118d0,-113222121d0, 252708246d0,   !   2373  1.00000
     . 75014500d0,  81872456d0,-300935985d0,-113498244d0, 252082704d0,   !   2374  1.00000
     . 74943421d0,  82075880d0,-300346479d0,-113774064d0, 251455726d0,   !   2375  1.00000
     . 74872172d0,  82279122d0,-299755602d0,-114049580d0, 250827315d0,   !   2376  1.00000
     . 74800752d0,  82482181d0,-299163356d0,-114324791d0, 250197475d0,   !   2377  1.00000
     . 74729162d0,  82685056d0,-298569743d0,-114599696d0, 249566207d0,   !   2378  1.00000
     . 74657402d0,  82887747d0,-297974765d0,-114874294d0, 248933513d0,   !   2379  1.00000
     . 74585473d0,  83090254d0,-297378424d0,-115148585d0, 248299397d0,   !   2380  1.00000
     . 74513374d0,  83292575d0,-296780723d0,-115422568d0, 247663860d0,   !   2381  1.00000
     . 74441106d0,  83494711d0,-296181663d0,-115696241d0, 247026906d0,   !   2382  1.00000
     . 74368669d0,  83696661d0,-295581246d0,-115969605d0, 246388535d0,   !   2383  1.00000
     . 74296063d0,  83898424d0,-294979475d0,-116242658d0, 245748752d0,   !   2384  1.00000
     . 74223289d0,  84100001d0,-294376352d0,-116515399d0, 245107558d0,   !   2385  1.00000
     . 74150346d0,  84301389d0,-293771878d0,-116787828d0, 244464956d0,   !   2386  1.00000
     . 74077235d0,  84502590d0,-293166057d0,-117059943d0, 243820948d0,   !   2387  1.00000
     . 74003956d0,  84703602d0,-292558890d0,-117331744d0, 243175537d0,   !   2388  1.00000
     . 73930509d0,  84904425d0,-291950378d0,-117603231d0, 242528725d0,   !   2389  1.00000
     . 73856895d0,  85105058d0,-291340525d0,-117874401d0, 241880515d0,   !   2390  1.00000
     . 73783114d0,  85305501d0,-290729333d0,-118145255d0, 241230909d0,   !   2391  1.00000
     . 73709165d0,  85505754d0,-290116803d0,-118415791d0, 240579910d0,   !   2392  1.00000
     . 73635050d0,  85705816d0,-289502938d0,-118686009d0, 239927520d0,   !   2393  1.00000
     . 73560768d0,  85905686d0,-288887740d0,-118955908d0, 239273743d0,   !   2394  1.00000
     . 73486319d0,  86105364d0,-288271212d0,-119225487d0, 238618579d0,   !   2395  1.00000
     . 73411705d0,  86304849d0,-287653354d0,-119494745d0, 237962033d0,   !   2396  1.00000
     . 73336924d0,  86504141d0,-287034170d0,-119763682d0, 237304105d0,   !   2397  1.00000
     . 73261978d0,  86703240d0,-286413662d0,-120032296d0, 236644800d0,   !   2398  1.00000
     . 73186865d0,  86902145d0,-285791831d0,-120300587d0, 235984120d0,   !   2399  1.00000
     . 73111588d0,  87100855d0,-285168680d0,-120568553d0, 235322066d0,   !   2400  1.00000
     . 73036146d0,  87299369d0,-284544212d0,-120836195d0, 234658642d0,   !   2401  1.00000
     . 72960538d0,  87497689d0,-283918428d0,-121103511d0, 233993850d0,   !   2402  1.00000
     . 72884766d0,  87695812d0,-283291331d0,-121370500d0, 233327694d0,   !   2403  1.00000
     . 72808829d0,  87893739d0,-282662923d0,-121637161d0, 232660174d0,   !   2404  1.00000
     . 72732729d0,  88091468d0,-282033206d0,-121903494d0, 231991295d0,   !   2405  1.00000
     . 72656464d0,  88289000d0,-281402182d0,-122169498d0, 231321058d0,   !   2406  1.00000
     . 72580035d0,  88486334d0,-280769854d0,-122435172d0, 230649466d0,   !   2407  1.00000
     . 72503443d0,  88683469d0,-280136224d0,-122700516d0, 229976522d0,   !   2408  1.00000
     . 72426687d0,  88880405d0,-279501294d0,-122965527d0, 229302229d0,   !   2409  1.00000
     . 72349769d0,  89077142d0,-278865066d0,-123230206d0, 228626589d0,   !   2410  1.00000
     . 72272687d0,  89273678d0,-278227542d0,-123494552d0, 227949604d0,   !   2411  1.00000
     . 72195443d0,  89470014d0,-277588726d0,-123758563d0, 227271278d0,   !   2412  1.00000
     . 72118036d0,  89666148d0,-276948619d0,-124022240d0, 226591612d0,   !   2413  1.00000
     . 72040468d0,  89862082d0,-276307223d0,-124285580d0, 225910611d0,   !   2414  1.00000
     . 71962737d0,  90057812d0,-275664540d0,-124548584d0, 225228275d0,   !   2415  1.00000
     . 71884844d0,  90253341d0,-275020574d0,-124811251d0, 224544609d0,   !   2416  1.00000
     . 71806790d0,  90448666d0,-274375326d0,-125073579d0, 223859614d0,   !   2417  1.00000
     . 71728575d0,  90643788d0,-273728799d0,-125335568d0, 223173294d0,   !   2418  1.00000
     . 71650198d0,  90838705d0,-273080994d0,-125597217d0, 222485650d0,   !   2419  1.00000
     . 71571661d0,  91033418d0,-272431914d0,-125858524d0, 221796686d0,   !   2420  1.00000
     . 71492963d0,  91227926d0,-271781562d0,-126119491d0, 221106405d0,   !   2421  1.00000
     . 71414105d0,  91422228d0,-271129940d0,-126380114d0, 220414809d0,   !   2422  1.00000
     . 71335086d0,  91616324d0,-270477050d0,-126640395d0, 219721900d0,   !   2423  1.00000
     . 71255908d0,  91810213d0,-269822895d0,-126900331d0, 219027683d0,   !   2424  1.00000
     . 71176570d0,  92003895d0,-269167476d0,-127159922d0, 218332158d0,   !   2425  1.00000
     . 71097073d0,  92197370d0,-268510796d0,-127419167d0, 217635329d0,   !   2426  1.00000
     . 71017416d0,  92390636d0,-267852858d0,-127678066d0, 216937199d0,   !   2427  1.00000
     . 70937601d0,  92583694d0,-267193664d0,-127936617d0, 216237771d0,   !   2428  1.00000
     . 70857626d0,  92776543d0,-266533216d0,-128194819d0, 215537047d0,   !   2429  1.00000
     . 70777494d0,  92969182d0,-265871516d0,-128452673d0, 214835030d0,   !   2430  1.00000
     . 70697202d0,  93161611d0,-265208568d0,-128710176d0, 214131722d0,   !   2431  1.00000
     . 70616753d0,  93353829d0,-264544372d0,-128967328d0, 213427127d0,   !   2432  1.00000
     . 70536147d0,  93545836d0,-263878933d0,-129224129d0, 212721248d0,   !   2433  1.00000
     . 70455382d0,  93737631d0,-263212251d0,-129480577d0, 212014086d0,   !   2434  1.00000
     . 70374461d0,  93929215d0,-262544330d0,-129736672d0, 211305646d0,   !   2435  1.00000
     . 70293382d0,  94120585d0,-261875172d0,-129992412d0, 210595929d0,   !   2436  1.00000
     . 70212146d0,  94311743d0,-261204780d0,-130247798d0, 209884939d0,   !   2437  1.00000
     . 70130754d0,  94502687d0,-260533155d0,-130502827d0, 209172677d0,   !   2438  1.00000
     . 70049206d0,  94693417d0,-259860300d0,-130757500d0, 208459148d0,   !   2439  1.00000
     . 69967501d0,  94883932d0,-259186217d0,-131011815d0, 207744354d0,   !   2440  1.00000
     . 69885641d0,  95074232d0,-258510910d0,-131265772d0, 207028298d0,   !   2441  1.00000
     . 69803625d0,  95264316d0,-257834380d0,-131519369d0, 206310983d0,   !   2442  1.00000
     . 69721454d0,  95454184d0,-257156630d0,-131772607d0, 205592410d0,   !   2443  1.00000
     . 69639127d0,  95643835d0,-256477662d0,-132025483d0, 204872585d0,   !   2444  1.00000
     . 69556646d0,  95833270d0,-255797479d0,-132277998d0, 204151508d0,   !   2445  1.00000
     . 69474010d0,  96022486d0,-255116083d0,-132530150d0, 203429183d0,   !   2446  1.00000
     . 69391220d0,  96211485d0,-254433477d0,-132781939d0, 202705614d0,   !   2447  1.00000
     . 69308275d0,  96400265d0,-253749663d0,-133033364d0, 201980802d0,   !   2448  1.00000
     . 69225177d0,  96588825d0,-253064643d0,-133284423d0, 201254750d0,   !   2449  1.00000
     . 69141925d0,  96777166d0,-252378421d0,-133535116d0, 200527463d0,   !   2450  1.00000
     . 69058519d0,  96965287d0,-251690998d0,-133785443d0, 199798942d0,   !   2451  1.00000
     . 68974961d0,  97153187d0,-251002378d0,-134035402d0, 199069190d0,   !   2452  1.00000
     . 68891249d0,  97340866d0,-250312561d0,-134284993d0, 198338210d0,   !   2453  1.00000
     . 68807385d0,  97528323d0,-249621552d0,-134534214d0, 197606006d0,   !   2454  1.00000
     . 68723369d0,  97715558d0,-248929353d0,-134783066d0, 196872580d0,   !   2455  1.00000
     . 68639200d0,  97902571d0,-248235965d0,-135031546d0, 196137934d0,   !   2456  1.00000
     . 68554880d0,  98089360d0,-247541393d0,-135279655d0, 195402073d0,   !   2457  1.00000
     . 68470408d0,  98275926d0,-246845637d0,-135527391d0, 194664999d0,   !   2458  1.00000
     . 68385785d0,  98462267d0,-246148701d0,-135774753d0, 193926715d0,   !   2459  1.00000
     . 68301010d0,  98648384d0,-245450587d0,-136021741d0, 193187223d0,   !   2460  1.00000
     . 68216085d0,  98834275d0,-244751298d0,-136268355d0, 192446527d0,   !   2461  1.00000
     . 68131008d0,  99019941d0,-244050836d0,-136514592d0, 191704630d0,   !   2462  1.00000
     . 68045782d0,  99205381d0,-243349204d0,-136760452d0, 190961535d0,   !   2463  1.00000
     . 67960406d0,  99390594d0,-242646405d0,-137005935d0, 190217245d0,   !   2464  1.00000
     . 67874879d0,  99575579d0,-241942440d0,-137251039d0, 189471762d0,   !   2465  1.00000
     . 67789203d0,  99760337d0,-241237313d0,-137495764d0, 188725090d0,   !   2466  1.00000
     . 67703378d0,  99944867d0,-240531026d0,-137740109d0, 187977232d0,   !   2467  1.00000
     . 67617404d0, 100129168d0,-239823582d0,-137984073d0, 187228190d0,   !   2468  1.00000
     . 67531280d0, 100313240d0,-239114983d0,-138227655d0, 186477969d0,   !   2469  1.00000
     . 67445009d0, 100497083d0,-238405232d0,-138470854d0, 185726570d0,   !   2470  1.00000
     . 67358589d0, 100680695d0,-237694331d0,-138713670d0, 184973996d0,   !   2471  1.00000
     . 67272021d0, 100864076d0,-236982284d0,-138956102d0, 184220252d0,   !   2472  1.00000
     . 67185305d0, 101047226d0,-236269092d0,-139198148d0, 183465340d0,   !   2473  1.00000
     . 67098441d0, 101230145d0,-235554758d0,-139439809d0, 182709262d0,   !   2474  1.00000
     . 67011431d0, 101412831d0,-234839286d0,-139681083d0, 181952023d0,   !   2475  1.00000
     . 66924273d0, 101595285d0,-234122677d0,-139921969d0, 181193624d0,   !   2476  1.00000
     . 66836969d0, 101777505d0,-233404934d0,-140162466d0, 180434070d0,   !   2477  1.00000
     . 66749518d0, 101959492d0,-232686060d0,-140402575d0, 179673363d0,   !   2478  1.00000
     . 66661921d0, 102141244d0,-231966057d0,-140642293d0, 178911506d0,   !   2479  1.00000
     . 66574178d0, 102322762d0,-231244929d0,-140881620d0, 178148503d0,   !   2480  1.00000
     . 66486290d0, 102504044d0,-230522677d0,-141120555d0, 177384356d0,   !   2481  1.00000
     . 66398256d0, 102685091d0,-229799305d0,-141359098d0, 176619069d0,   !   2482  1.00000
     . 66310077d0, 102865902d0,-229074814d0,-141597247d0, 175852644d0,   !   2483  1.00000
     . 66221753d0, 103046475d0,-228349209d0,-141835002d0, 175085085d0,   !   2484  1.00000
     . 66133284d0, 103226812d0,-227622491d0,-142072361d0, 174316395d0,   !   2485  1.00000
     . 66044671d0, 103406911d0,-226894663d0,-142309325d0, 173546578d0,   !   2486  1.00000
     . 65955915d0, 103586772d0,-226165728d0,-142545892d0, 172775635d0,   !   2487  1.00000
     . 65867014d0, 103766393d0,-225435689d0,-142782061d0, 172003571d0,   !   2488  1.00000
     . 65777970d0, 103945776d0,-224704547d0,-143017831d0, 171230388d0,   !   2489  1.00000
     . 65688782d0, 104124919d0,-223972307d0,-143253202d0, 170456090d0,   !   2490  1.00000
     . 65599452d0, 104303822d0,-223238970d0,-143488173d0, 169680679d0,   !   2491  1.00000
     . 65509978d0, 104482484d0,-222504540d0,-143722743d0, 168904160d0,   !   2492  1.00000
     . 65420363d0, 104660905d0,-221769019d0,-143956912d0, 168126535d0,   !   2493  1.00000
     . 65330605d0, 104839084d0,-221032409d0,-144190677d0, 167347807d0,   !   2494  1.00000
     . 65240706d0, 105017021d0,-220294714d0,-144424039d0, 166567979d0,   !   2495  1.00000
     . 65150664d0, 105194715d0,-219555936d0,-144656997d0, 165787056d0,   !   2496  1.00000
     . 65060482d0, 105372165d0,-218816079d0,-144889549d0, 165005039d0,   !   2497  1.00000
     . 64970158d0, 105549373d0,-218075144d0,-145121695d0, 164221932d0,   !   2498  1.00000
     . 64879694d0, 105726335d0,-217333135d0,-145353435d0, 163437739d0,   !   2499  1.00000
     . 64789089d0, 105903054d0,-216590054d0,-145584767d0, 162652462d0,   !   2500  1.00000
     . 64698344d0, 106079526d0,-215845904d0,-145815690d0, 161866105d0,   !   2501  1.00000
     . 64607458d0, 106255754d0,-215100687d0,-146046204d0, 161078671d0,   !   2502  1.00000
     . 64516434d0, 106431734d0,-214354408d0,-146276307d0, 160290164d0,   !   2503  1.00000
     . 64425269d0, 106607468d0,-213607068d0,-146506000d0, 159500586d0,   !   2504  1.00000
     . 64333966d0, 106782955d0,-212858670d0,-146735281d0, 158709940d0,   !   2505  1.00000
     . 64242524d0, 106958194d0,-212109217d0,-146964149d0, 157918231d0,   !   2506  1.00000
     . 64150943d0, 107133185d0,-211358712d0,-147192603d0, 157125461d0,   !   2507  1.00000
     . 64059224d0, 107307927d0,-210607157d0,-147420643d0, 156331634d0,   !   2508  1.00000
     . 63967367d0, 107482420d0,-209854556d0,-147648268d0, 155536752d0,   !   2509  1.00000
     . 63875373d0, 107656663d0,-209100911d0,-147875477d0, 154740820d0,   !   2510  1.00000
     . 63783240d0, 107830656d0,-208346226d0,-148102269d0, 153943840d0,   !   2511  1.00000
     . 63690971d0, 108004398d0,-207590502d0,-148328643d0, 153145815d0,   !   2512  1.00000
     . 63598565d0, 108177889d0,-206833743d0,-148554599d0, 152346750d0,   !   2513  1.00000
     . 63506022d0, 108351128d0,-206075951d0,-148780135d0, 151546647d0,   !   2514  1.00000
     . 63413343d0, 108524115d0,-205317130d0,-149005251d0, 150745510d0,   !   2515  1.00000
     . 63320528d0, 108696849d0,-204557282d0,-149229946d0, 149943342d0,   !   2516  1.00000
     . 63227578d0, 108869329d0,-203796410d0,-149454219d0, 149140146d0,   !   2517  1.00000
     . 63134491d0, 109041556d0,-203034518d0,-149678069d0, 148335926d0,   !   2518  1.00000
     . 63041270d0, 109213529d0,-202271607d0,-149901496d0, 147530685d0,   !   2519  1.00000
     . 62947914d0, 109385247d0,-201507681d0,-150124498d0, 146724427d0,   !   2520  1.00000
     . 62854423d0, 109556710d0,-200742742d0,-150347075d0, 145917154d0,   !   2521  1.00000
     . 62760798d0, 109727916d0,-199976794d0,-150569226d0, 145108870d0,   !   2522  1.00000
     . 62667039d0, 109898867d0,-199209840d0,-150790950d0, 144299579d0,   !   2523  1.00000
     . 62573146d0, 110069561d0,-198441881d0,-151012247d0, 143489283d0,   !   2524  1.00000
     . 62479120d0, 110239997d0,-197672922d0,-151233114d0, 142677987d0,   !   2525  1.00000
     . 62384961d0, 110410176d0,-196902965d0,-151453553d0, 141865694d0,   !   2526  1.00000
     . 62290669d0, 110580096d0,-196132013d0,-151673561d0, 141052406d0,   !   2527  1.00000
     . 62196244d0, 110749758d0,-195360070d0,-151893137d0, 140238128d0,   !   2528  1.00000
     . 62101687d0, 110919160d0,-194587136d0,-152112282d0, 139422863d0,   !   2529  1.00000
     . 62006998d0, 111088303d0,-193813217d0,-152330995d0, 138606614d0,   !   2530  1.00000
     . 61912178d0, 111257185d0,-193038315d0,-152549273d0, 137789385d0,   !   2531  1.00000
     . 61817226d0, 111425807d0,-192262432d0,-152767118d0, 136971179d0,   !   2532  1.00000
     . 61722143d0, 111594167d0,-191485572d0,-152984526d0, 136152000d0,   !   2533  1.00000
     . 61626930d0, 111762265d0,-190707738d0,-153201499d0, 135331850d0,   !   2534  1.00000
     . 61531585d0, 111930101d0,-189928932d0,-153418035d0, 134510734d0,   !   2535  1.00000
     . 61436111d0, 112097675d0,-189149157d0,-153634133d0, 133688656d0,   !   2536  1.00000
     . 61340507d0, 112264985d0,-188368417d0,-153849793d0, 132865617d0,   !   2537  1.00000
     . 61244773d0, 112432031d0,-187586715d0,-154065013d0, 132041623d0,   !   2538  1.00000
     . 61148910d0, 112598813d0,-186804053d0,-154279792d0, 131216675d0,   !   2539  1.00000
     . 61052918d0, 112765330d0,-186020435d0,-154494131d0, 130390779d0,   !   2540  1.00000
     . 60956797d0, 112931582d0,-185235863d0,-154708028d0, 129563937d0,   !   2541  1.00000
     . 60860548d0, 113097568d0,-184450340d0,-154921481d0, 128736153d0,   !   2542  1.00000
     . 60764171d0, 113263287d0,-183663870d0,-155134492d0, 127907430d0,   !   2543  1.00000
     . 60667666d0, 113428740d0,-182876456d0,-155347058d0, 127077772d0,   !   2544  1.00000
     . 60571033d0, 113593926d0,-182088100d0,-155559178d0, 126247183d0,   !   2545  1.00000
     . 60474273d0, 113758843d0,-181298805d0,-155770852d0, 125415665d0,   !   2546  1.00000
     . 60377387d0, 113923493d0,-180508575d0,-155982080d0, 124583223d0,   !   2547  1.00000
     . 60280374d0, 114087873d0,-179717413d0,-156192859d0, 123749860d0,   !   2548  1.00000
     . 60183234d0, 114251985d0,-178925321d0,-156403190d0, 122915579d0,   !   2549  1.00000
     . 60085969d0, 114415826d0,-178132303d0,-156613072d0, 122080385d0,   !   2550  1.00000
     . 59988578d0, 114579397d0,-177338362d0,-156822503d0, 121244280d0,   !   2551  1.00000
     . 59891061d0, 114742697d0,-176543500d0,-157031483d0, 120407268d0,   !   2552  1.00000
     . 59793420d0, 114905726d0,-175747722d0,-157240011d0, 119569353d0,   !   2553  1.00000
     . 59695654d0, 115068483d0,-174951029d0,-157448086d0, 118730538d0,   !   2554  1.00000
     . 59597763d0, 115230968d0,-174153425d0,-157655708d0, 117890828d0,   !   2555  1.00000
     . 59499748d0, 115393180d0,-173354913d0,-157862875d0, 117050225d0,   !   2556  1.00000
     . 59401610d0, 115555119d0,-172555497d0,-158069586d0, 116208733d0,   !   2557  1.00000
     . 59303348d0, 115716784d0,-171755178d0,-158275842d0, 115366356d0,   !   2558  1.00000
     . 59204963d0, 115878174d0,-170953961d0,-158481640d0, 114523097d0,   !   2559  1.00000
     . 59106455d0, 116039290d0,-170151849d0,-158686981d0, 113678960d0,   !   2560  1.00000
     . 59007824d0, 116200130d0,-169348844d0,-158891863d0, 112833948d0,   !   2561  1.00000
     . 58909071d0, 116360695d0,-168544949d0,-159096286d0, 111988066d0,   !   2562  1.00000
     . 58810196d0, 116520983d0,-167740169d0,-159300248d0, 111141317d0,   !   2563  1.00000
     . 58711199d0, 116680994d0,-166934505d0,-159503749d0, 110293705d0,   !   2564  1.00000
     . 58612082d0, 116840728d0,-166127962d0,-159706788d0, 109445232d0,   !   2565  1.00000
     . 58512843d0, 117000184d0,-165320541d0,-159909364d0, 108595903d0,   !   2566  1.00000
     . 58413483d0, 117159362d0,-164512247d0,-160111477d0, 107745722d0,   !   2567  1.00000
     . 58314003d0, 117318261d0,-163703083d0,-160313125d0, 106894692d0,   !   2568  1.00000
     . 58214403d0, 117476880d0,-162893051d0,-160514308d0, 106042817d0,   !   2569  1.00000
     . 58114683d0, 117635220d0,-162082155d0,-160715024d0, 105190100d0,   !   2570  1.00000
     . 58014844d0, 117793279d0,-161270397d0,-160915274d0, 104336546d0,   !   2571  1.00000
     . 57914886d0, 117951058d0,-160457782d0,-161115055d0, 103482157d0,   !   2572  1.00000
     . 57814809d0, 118108555d0,-159644313d0,-161314368d0, 102626938d0,   !   2573  1.00000
     . 57714613d0, 118265770d0,-158829991d0,-161513212d0, 101770892d0,   !   2574  1.00000
     . 57614299d0, 118422703d0,-158014822d0,-161711585d0, 100914024d0,   !   2575  1.00000
     . 57513868d0, 118579353d0,-157198807d0,-161909487d0, 100056336d0,   !   2576  1.00000
     . 57413319d0, 118735720d0,-156381950d0,-162106916d0,  99197832d0,   !   2577  1.00000
     . 57312652d0, 118891803d0,-155564255d0,-162303873d0,  98338517d0,   !   2578  1.00000
     . 57211869d0, 119047601d0,-154745724d0,-162500357d0,  97478394d0,   !   2579  1.00000
     . 57110969d0, 119203115d0,-153926361d0,-162696365d0,  96617467d0,   !   2580  1.00000
     . 57009953d0, 119358343d0,-153106169d0,-162891899d0,  95755739d0,   !   2581  1.00000
     . 56908822d0, 119513286d0,-152285151d0,-163086956d0,  94893214d0,   !   2582  1.00000
     . 56807574d0, 119667942d0,-151463310d0,-163281536d0,  94029896d0,   !   2583  1.00000
     . 56706211d0, 119822311d0,-150640650d0,-163475638d0,  93165790d0,   !   2584  1.00000
     . 56604734d0, 119976392d0,-149817173d0,-163669261d0,  92300897d0,   !   2585  1.00000
     . 56503141d0, 120130186d0,-148992884d0,-163862405d0,  91435223d0,   !   2586  1.00000
     . 56401435d0, 120283692d0,-148167786d0,-164055069d0,  90568772d0,   !   2587  1.00000
     . 56299614d0, 120436908d0,-147341880d0,-164247251d0,  89701546d0,   !   2588  1.00000
     . 56197680d0, 120589836d0,-146515172d0,-164438951d0,  88833550d0,   !   2589  1.00000
     . 56095633d0, 120742473d0,-145687664d0,-164630169d0,  87964787d0,   !   2590  1.00000
     . 55993472d0, 120894820d0,-144859360d0,-164820902d0,  87095262d0,   !   2591  1.00000
     . 55891199d0, 121046876d0,-144030262d0,-165011152d0,  86224979d0,   !   2592  1.00000
     . 55788814d0, 121198641d0,-143200374d0,-165200915d0,  85353940d0,   !   2593  1.00000
     . 55686316d0, 121350114d0,-142369699d0,-165390193d0,  84482150d0,   !   2594  1.00000
     . 55583707d0, 121501294d0,-141538242d0,-165578983d0,  83609613d0,   !   2595  1.00000
     . 55480987d0, 121652182d0,-140706004d0,-165767285d0,  82736333d0,   !   2596  1.00000
     . 55378155d0, 121802776d0,-139872989d0,-165955099d0,  81862313d0,   !   2597  1.00000
     . 55275213d0, 121953076d0,-139039201d0,-166142423d0,  80987558d0,   !   2598  1.00000
     . 55172161d0, 122103082d0,-138204643d0,-166329257d0,  80112070d0,   !   2599  1.00000
     . 55068998d0, 122252793d0,-137369318d0,-166515600d0,  79235855d0,   !   2600  1.00000
     . 54965726d0, 122402208d0,-136533229d0,-166701450d0,  78358916d0,   !   2601  1.00000
     . 54862345d0, 122551328d0,-135696381d0,-166886807d0,  77481257d0,   !   2602  1.00000
     . 54758854d0, 122700151d0,-134858775d0,-167071671d0,  76602882d0,   !   2603  1.00000
     . 54655255d0, 122848678d0,-134020417d0,-167256040d0,  75723794d0,   !   2604  1.00000
     . 54551548d0, 122996906d0,-133181308d0,-167439914d0,  74843998d0,   !   2605  1.00000
     . 54447732d0, 123144838d0,-132341452d0,-167623291d0,  73963498d0,   !   2606  1.00000
     . 54343809d0, 123292470d0,-131500854d0,-167806172d0,  73082297d0,   !   2607  1.00000
     . 54239779d0, 123439804d0,-130659515d0,-167988554d0,  72200399d0,   !   2608  1.00000
     . 54135642d0, 123586838d0,-129817440d0,-168170438d0,  71317809d0,   !   2609  1.00000
     . 54031398d0, 123733573d0,-128974631d0,-168351822d0,  70434530d0,   !   2610  1.00000
     . 53927047d0, 123880007d0,-128131093d0,-168532705d0,  69550567d0,   !   2611  1.00000
     . 53822591d0, 124026140d0,-127286828d0,-168713088d0,  68665922d0,   !   2612  1.00000
     . 53718029d0, 124171971d0,-126441841d0,-168892968d0,  67780601d0,   !   2613  1.00000
     . 53613362d0, 124317501d0,-125596134d0,-169072345d0,  66894607d0,   !   2614  1.00000
     . 53508590d0, 124462728d0,-124749710d0,-169251219d0,  66007945d0,   !   2615  1.00000
     . 53403713d0, 124607653d0,-123902574d0,-169429588d0,  65120617d0,   !   2616  1.00000
     . 53298732d0, 124752274d0,-123054729d0,-169607451d0,  64232629d0,   !   2617  1.00000
     . 53193648d0, 124896591d0,-122206178d0,-169784808d0,  63343983d0,   !   2618  1.00000
     . 53088459d0, 125040603d0,-121356924d0,-169961658d0,  62454685d0,   !   2619  1.00000
     . 52983168d0, 125184311d0,-120506971d0,-170138000d0,  61564738d0,   !   2620  1.00000
     . 52877773d0, 125327713d0,-119656323d0,-170313834d0,  60674146d0,   !   2621  1.00000
     . 52772276d0, 125470810d0,-118804983d0,-170489157d0,  59782914d0,   !   2622  1.00000
     . 52666677d0, 125613600d0,-117952954d0,-170663970d0,  58891044d0,   !   2623  1.00000
     . 52560976d0, 125756083d0,-117100240d0,-170838272d0,  57998542d0,   !   2624  1.00000
     . 52455173d0, 125898258d0,-116246844d0,-171012062d0,  57105412d0,   !   2625  1.00000
     . 52349270d0, 126040126d0,-115392770d0,-171185339d0,  56211656d0,   !   2626  1.00000
     . 52243265d0, 126181685d0,-114538021d0,-171358102d0,  55317281d0,   !   2627  1.00000
     . 52137160d0, 126322935d0,-113682601d0,-171530350d0,  54422288d0,   !   2628  1.00000
     . 52030955d0, 126463876d0,-112826513d0,-171702083d0,  53526683d0,   !   2629  1.00000
     . 51924651d0, 126604507d0,-111969761d0,-171873299d0,  52630470d0,   !   2630  1.00000
     . 51818247d0, 126744827d0,-111112349d0,-172043999d0,  51733653d0,   !   2631  1.00000
     . 51711743d0, 126884837d0,-110254279d0,-172214180d0,  50836236d0,   !   2632  1.00000
     . 51605141d0, 127024535d0,-109395555d0,-172383843d0,  49938222d0,   !   2633  1.00000
     . 51498441d0, 127163921d0,-108536181d0,-172552986d0,  49039617d0,   !   2634  1.00000
     . 51391643d0, 127302994d0,-107676160d0,-172721608d0,  48140424d0,   !   2635  1.00000
     . 51284747d0, 127441755d0,-106815497d0,-172889709d0,  47240647d0,   !   2636  1.00000
     . 51177754d0, 127580202d0,-105954193d0,-173057288d0,  46340291d0,   !   2637  1.00000
     . 51070664d0, 127718335d0,-105092254d0,-173224344d0,  45439359d0,   !   2638  1.00000
     . 50963478d0, 127856154d0,-104229682d0,-173390876d0,  44537856d0,   !   2639  1.00000
     . 50856195d0, 127993658d0,-103366481d0,-173556884d0,  43635785d0,   !   2640  1.00000
     . 50748816d0, 128130846d0,-102502655d0,-173722366d0,  42733152d0,   !   2641  1.00000
     . 50641342d0, 128267719d0,-101638207d0,-173887322d0,  41829960d0,   !   2642  1.00000
     . 50533773d0, 128404275d0,-100773141d0,-174051750d0,  40926213d0,   !   2643  1.00000
     . 50426109d0, 128540514d0, -99907460d0,-174215650d0,  40021916d0,   !   2644  1.00000
     . 50318350d0, 128676435d0, -99041168d0,-174379022d0,  39117072d0,   !   2645  1.00000
     . 50210498d0, 128812039d0, -98174269d0,-174541864d0,  38211687d0,   !   2646  1.00000
     . 50102552d0, 128947324d0, -97306765d0,-174704175d0,  37305763d0,   !   2647  1.00000
     . 49994512d0, 129082290d0, -96438662d0,-174865955d0,  36399306d0,   !   2648  1.00000
     . 49886379d0, 129216936d0, -95569961d0,-175027203d0,  35492319d0,   !   2649  1.00000
     . 49778154d0, 129351263d0, -94700668d0,-175187918d0,  34584807d0,   !   2650  1.00000
     . 49669837d0, 129485269d0, -93830785d0,-175348099d0,  33676773d0,   !   2651  1.00000
     . 49561427d0, 129618955d0, -92960316d0,-175507746d0,  32768223d0,   !   2652  1.00000
     . 49452926d0, 129752318d0, -92089265d0,-175666856d0,  31859160d0,   !   2653  1.00000
     . 49344334d0, 129885360d0, -91217636d0,-175825431d0,  30949589d0,   !   2654  1.00000
     . 49235652d0, 130018079d0, -90345432d0,-175983468d0,  30039514d0,   !   2655  1.00000
     . 49126878d0, 130150476d0, -89472656d0,-176140967d0,  29128939d0,   !   2656  1.00000
     . 49018015d0, 130282548d0, -88599313d0,-176297927d0,  28217868d0,   !   2657  1.00000
     . 48909062d0, 130414297d0, -87725406d0,-176454348d0,  27306306d0,   !   2658  1.00000
     . 48800019d0, 130545722d0, -86850938d0,-176610228d0,  26394257d0,   !   2659  1.00000
     . 48690888d0, 130676821d0, -85975914d0,-176765567d0,  25481725d0,   !   2660  1.00000
     . 48581668d0, 130807595d0, -85100337d0,-176920363d0,  24568714d0,   !   2661  1.00000
     . 48472360d0, 130938043d0, -84224211d0,-177074617d0,  23655229d0,   !   2662  1.00000
     . 48362964d0, 131068164d0, -83347539d0,-177228326d0,  22741274d0,   !   2663  1.00000
     . 48253480d0, 131197958d0, -82470325d0,-177381491d0,  21826854d0,   !   2664  1.00000
     . 48143909d0, 131327425d0, -81592573d0,-177534110d0,  20911972d0,   !   2665  1.00000
     . 48034252d0, 131456564d0, -80714287d0,-177686183d0,  19996633d0,   !   2666  1.00000
     . 47924508d0, 131585374d0, -79835469d0,-177837709d0,  19080842d0,   !   2667  1.00000
     . 47814678d0, 131713856d0, -78956125d0,-177988687d0,  18164602d0,   !   2668  1.00000
     . 47704762d0, 131842008d0, -78076258d0,-178139115d0,  17247919d0,   !   2669  1.00000
     . 47594762d0, 131969830d0, -77195870d0,-178288994d0,  16330795d0,   !   2670  1.00000
     . 47484676d0, 132097321d0, -76314967d0,-178438323d0,  15413236d0,   !   2671  1.00000
     . 47374506d0, 132224481d0, -75433552d0,-178587100d0,  14495247d0,   !   2672  1.00000
     . 47264251d0, 132351310d0, -74551628d0,-178735325d0,  13576831d0,   !   2673  1.00000
     . 47153913d0, 132477807d0, -73669200d0,-178882996d0,  12657992d0,   !   2674  1.00000
     . 47043491d0, 132603971d0, -72786271d0,-179030114d0,  11738736d0,   !   2675  1.00000
     . 46932987d0, 132729803d0, -71902844d0,-179176677d0,  10819066d0,   !   2676  1.00000
     . 46822400d0, 132855300d0, -71018925d0,-179322685d0,   9898987d0,   !   2677  1.00000
     . 46711730d0, 132980464d0, -70134515d0,-179468136d0,   8978504d0,   !   2678  1.00000
     . 46600979d0, 133105294d0, -69249620d0,-179613030d0,   8057620d0,   !   2679  1.00000
     . 46490146d0, 133229788d0, -68364243d0,-179757365d0,   7136341d0,   !   2680  1.00000
     . 46379232d0, 133353947d0, -67478387d0,-179901142d0,   6214670d0,   !   2681  1.00000
     . 46268237d0, 133477769d0, -66592058d0,-180044359d0,   5292612d0,   !   2682  1.00000
     . 46157162d0, 133601255d0, -65705257d0,-180187016d0,   4370171d0,   !   2683  1.00000
     . 46046007d0, 133724405d0, -64817990d0,-180329111d0,   3447353d0,   !   2684  1.00000
     . 45934772d0, 133847216d0, -63930259d0,-180470644d0,   2524160d0,   !   2685  1.00000
     . 45823458d0, 133969690d0, -63042070d0,-180611613d0,   1600599d0,   !   2686  1.00000
     . 45712065d0, 134091825d0, -62153425d0,-180752019d0,    676673d0,   !   2687  1.00000
     . 45600594d0, 134213621d0, -61264328d0,-180891860d0,   -247613d0,   !   2688  1.00000
     . 45489044d0, 134335077d0, -60374784d0,-181031135d0,  -1172256d0,   !   2689  1.00000
     . 45377417d0, 134456194d0, -59484795d0,-181169844d0,  -2097249d0,   !   2690  1.00000
     . 45265713d0, 134576970d0, -58594367d0,-181307986d0,  -3022590d0,   !   2691  1.00000
     . 45153931d0, 134697405d0, -57703503d0,-181445560d0,  -3948273d0,   !   2692  1.00000
     . 45042073d0, 134817498d0, -56812206d0,-181582564d0,  -4874294d0,   !   2693  1.00000
     . 44930139d0, 134937249d0, -55920480d0,-181718999d0,  -5800648d0,   !   2694  1.00000
     . 44818130d0, 135056658d0, -55028330d0,-181854863d0,  -6727332d0,   !   2695  1.00000
     . 44706044d0, 135175724d0, -54135760d0,-181990156d0,  -7654339d0,   !   2696  1.00000
     . 44593884d0, 135294446d0, -53242772d0,-182124876d0,  -8581667d0,   !   2697  1.00000
     . 44481649d0, 135412824d0, -52349372d0,-182259024d0,  -9509309d0,   !   2698  1.00000
     . 44369341d0, 135530858d0, -51455562d0,-182392597d0, -10437263d0,   !   2699  1.00000
     . 44256958d0, 135648546d0, -50561348d0,-182525596d0, -11365523d0,   !   2700  1.00000
     . 44144501d0, 135765889d0, -49666732d0,-182658019d0, -12294084d0,   !   2701  1.00000
     . 44031972d0, 135882886d0, -48771718d0,-182789865d0, -13222942d0,   !   2702  1.00000
     . 43919370d0, 135999537d0, -47876311d0,-182921134d0, -14152094d0,   !   2703  1.00000
     . 43806696d0, 136115840d0, -46980515d0,-183051826d0, -15081533d0,   !   2704  1.00000
     . 43693950d0, 136231796d0, -46084333d0,-183181938d0, -16011255d0,   !   2705  1.00000
     . 43581132d0, 136347404d0, -45187769d0,-183311470d0, -16941256d0,   !   2706  1.00000
     . 43468244d0, 136462663d0, -44290828d0,-183440422d0, -17871532d0,   !   2707  1.00000
     . 43355285d0, 136577573d0, -43393512d0,-183568792d0, -18802077d0,   !   2708  1.00000
     . 43242255d0, 136692133d0, -42495827d0,-183696580d0, -19732887d0,   !   2709  1.00000
     . 43129156d0, 136806344d0, -41597776d0,-183823784d0, -20663958d0,   !   2710  1.00000
     . 43015987d0, 136920204d0, -40699363d0,-183950405d0, -21595284d0,   !   2711  1.00000
     . 42902748d0, 137033713d0, -39800592d0,-184076441d0, -22526861d0,   !   2712  1.00000
     . 42789442d0, 137146870d0, -38901467d0,-184201891d0, -23458686d0,   !   2713  1.00000
     . 42676067d0, 137259675d0, -38001992d0,-184326755d0, -24390752d0,   !   2714  1.00000
     . 42562624d0, 137372128d0, -37102171d0,-184451031d0, -25323055d0,   !   2715  1.00000
     . 42449113d0, 137484227d0, -36202008d0,-184574720d0, -26255591d0,   !   2716  1.00000
     . 42335535d0, 137595973d0, -35301507d0,-184697819d0, -27188355d0,   !   2717  1.00000
     . 42221891d0, 137707365d0, -34400672d0,-184820328d0, -28121342d0,   !   2718  1.00000
     . 42108181d0, 137818403d0, -33499506d0,-184942247d0, -29054548d0,   !   2719  1.00000
     . 41994404d0, 137929085d0, -32598015d0,-185063574d0, -29987968d0,   !   2720  1.00000
     . 41880562d0, 138039411d0, -31696201d0,-185184309d0, -30921597d0,   !   2721  1.00000
     . 41766655d0, 138149382d0, -30794070d0,-185304451d0, -31855430d0,   !   2722  1.00000
     . 41652683d0, 138258996d0, -29891624d0,-185423999d0, -32789464d0,   !   2723  1.00000
     . 41538647d0, 138368253d0, -28988869d0,-185542952d0, -33723693d0,   !   2724  1.00000
     . 41424547d0, 138477152d0, -28085807d0,-185661309d0, -34658112d0,   !   2725  1.00000
     . 41310383d0, 138585693d0, -27182444d0,-185779070d0, -35592717d0,   !   2726  1.00000
     . 41196157d0, 138693876d0, -26278783d0,-185896233d0, -36527503d0,   !   2727  1.00000
     . 41081867d0, 138801700d0, -25374828d0,-186012798d0, -37462465d0,   !   2728  1.00000
     . 40967516d0, 138909164d0, -24470584d0,-186128765d0, -38397599d0,   !   2729  1.00000
     . 40853102d0, 139016267d0, -23566053d0,-186244131d0, -39332900d0,   !   2730  1.00000
     . 40738628d0, 139123011d0, -22661241d0,-186358896d0, -40268363d0,   !   2731  1.00000
     . 40624092d0, 139229393d0, -21756152d0,-186473060d0, -41203983d0,   !   2732  1.00000
     . 40509495d0, 139335413d0, -20850789d0,-186586622d0, -42139756d0,   !   2733  1.00000
     . 40394838d0, 139441072d0, -19945157d0,-186699581d0, -43075676d0,   !   2734  1.00000
     . 40280122d0, 139546368d0, -19039260d0,-186811935d0, -44011740d0,   !   2735  1.00000
     . 40165346d0, 139651300d0, -18133101d0,-186923685d0, -44947942d0,   !   2736  1.00000
     . 40050511d0, 139755870d0, -17226686d0,-187034829d0, -45884277d0,   !   2737  1.00000
     . 39935617d0, 139860075d0, -16320017d0,-187145366d0, -46820741d0,   !   2738  1.00000
     . 39820665d0, 139963915d0, -15413100d0,-187255296d0, -47757329d0,   !   2739  1.00000
     . 39705656d0, 140067390d0, -14505938d0,-187364618d0, -48694036d0,   !   2740  1.00000
     . 39590589d0, 140170500d0, -13598536d0,-187473331d0, -49630857d0,   !   2741  1.00000
     . 39475465d0, 140273244d0, -12690897d0,-187581433d0, -50567788d0,   !   2742  1.00000
     . 39360285d0, 140375621d0, -11783026d0,-187688926d0, -51504823d0,   !   2743  1.00000
     . 39245048d0, 140477631d0, -10874927d0,-187795806d0, -52441958d0,   !   2744  1.00000
     . 39129756d0, 140579273d0,  -9966603d0,-187902074d0, -53379187d0,   !   2745  1.00000
     . 39014408d0, 140680547d0,  -9058060d0,-188007729d0, -54316507d0,   !   2746  1.00000
     . 38899006d0, 140781452d0,  -8149302d0,-188112770d0, -55253912d0,   !   2747  1.00000
     . 38783549d0, 140881989d0,  -7240332d0,-188217196d0, -56191398d0,   !   2748  1.00000
     . 38668039d0, 140982156d0,  -6331154d0,-188321006d0, -57128959d0,   !   2749  1.00000
     . 38552474d0, 141081952d0,  -5421773d0,-188424200d0, -58066590d0,   !   2750  1.00000
     . 38436857d0, 141181378d0,  -4512194d0,-188526776d0, -59004288d0,   !   2751  1.00000
     . 38321186d0, 141280433d0,  -3602420d0,-188628734d0, -59942046d0,   !   2752  1.00000
     . 38205464d0, 141379116d0,  -2692455d0,-188730073d0, -60879860d0,   !   2753  1.00000
     . 38089689d0, 141477427d0,  -1782303d0,-188830793d0, -61817725d0,   !   2754  1.00000
     . 37973863d0, 141575366d0,   -871970d0,-188930891d0, -62755636d0,   !   2755  1.00000
     . 37857986d0, 141672931d0,     38540d0,-189030368d0, -63693589d0,   !   2756  1.00000
     . 37742058d0, 141770123d0,    949225d0,-189129222d0, -64631578d0,   !   2757  1.00000
     . 37626080d0, 141866940d0,   1860080d0,-189227453d0, -65569598d0,   !   2758  1.00000
     . 37510052d0, 141963383d0,   2771100d0,-189325060d0, -66507644d0,   !   2759  1.00000
     . 37393974d0, 142059451d0,   3682280d0,-189422042d0, -67445713d0,   !   2760  1.00000
     . 37277848d0, 142155143d0,   4593618d0,-189518399d0, -68383797d0,   !   2761  1.00000
     . 37161673d0, 142250459d0,   5505108d0,-189614128d0, -69321894d0,   !   2762  1.00000
     . 37045450d0, 142345399d0,   6416746d0,-189709231d0, -70259997d0,   !   2763  1.00000
     . 36929180d0, 142439961d0,   7328529d0,-189803705d0, -71198101d0,   !   2764  1.00000
     . 36812862d0, 142534145d0,   8240451d0,-189897550d0, -72136203d0,   !   2765  1.00000
     . 36696497d0, 142627952d0,   9152509d0,-189990765d0, -73074296d0,   !   2766  1.00000
     . 36580086d0, 142721380d0,  10064697d0,-190083349d0, -74012376d0,   !   2767  1.00000
     . 36463628d0, 142814428d0,  10977013d0,-190175301d0, -74950437d0,   !   2768  1.00000
     . 36347126d0, 142907097d0,  11889452d0,-190266622d0, -75888475d0,   !   2769  1.00000
     . 36230578d0, 142999386d0,  12802008d0,-190357309d0, -76826485d0,   !   2770  1.00000
     . 36113985d0, 143091294d0,  13714679d0,-190447361d0, -77764462d0,   !   2771  1.00000
     . 35997349d0, 143182821d0,  14627460d0,-190536779d0, -78702401d0,   !   2772  1.00000
     . 35880668d0, 143273966d0,  15540346d0,-190625561d0, -79640296d0,   !   2773  1.00000
     . 35763944d0, 143364729d0,  16453333d0,-190713707d0, -80578143d0,   !   2774  1.00000
     . 35647177d0, 143455110d0,  17366417d0,-190801215d0, -81515936d0,   !   2775  1.00000
     . 35530367d0, 143545107d0,  18279593d0,-190888085d0, -82453671d0,   !   2776  1.00000
     . 35413516d0, 143634721d0,  19192857d0,-190974315d0, -83391343d0,   !   2777  1.00000
     . 35296622d0, 143723950d0,  20106206d0,-191059906d0, -84328946d0,   !   2778  1.00000
     . 35179688d0, 143812795d0,  21019633d0,-191144856d0, -85266475d0,   !   2779  1.00000
     . 35062712d0, 143901255d0,  21933136d0,-191229164d0, -86203926d0,   !   2780  1.00000
     . 34945697d0, 143989329d0,  22846710d0,-191312830d0, -87141293d0,   !   2781  1.00000
     . 34828641d0, 144077017d0,  23760350d0,-191395852d0, -88078571d0,   !   2782  1.00000
     . 34711546d0, 144164318d0,  24674051d0,-191478230d0, -89015755d0,   !   2783  1.00000
     . 34594412d0, 144251232d0,  25587811d0,-191559964d0, -89952840d0,   !   2784  1.00000
     . 34477239d0, 144337758d0,  26501624d0,-191641051d0, -90889821d0,   !   2785  1.00000
     . 34360028d0, 144423896d0,  27415486d0,-191721492d0, -91826693d0,   !   2786  1.00000
     . 34242779d0, 144509645d0,  28329392d0,-191801286d0, -92763450d0,   !   2787  1.00000
     . 34125493d0, 144595006d0,  29243338d0,-191880431d0, -93700088d0,   !   2788  1.00000
     . 34008170d0, 144679976d0,  30157320d0,-191958927d0, -94636601d0,   !   2789  1.00000
     . 33890811d0, 144764557d0,  31071333d0,-192036773d0, -95572984d0,   !   2790  1.00000
     . 33773416d0, 144848746d0,  31985374d0,-192113968d0, -96509233d0,   !   2791  1.00000
     . 33655985d0, 144932545d0,  32899436d0,-192190511d0, -97445341d0,   !   2792  1.00000
     . 33538519d0, 145015951d0,  33813517d0,-192266402d0, -98381304d0,   !   2793  1.00000
     . 33421019d0, 145098966d0,  34727611d0,-192341640d0, -99317117d0,   !   2794  1.00000
     . 33303484d0, 145181588d0,  35641715d0,-192416223d0,-100252774d0,   !   2795  1.00000
     . 33185916d0, 145263817d0,  36555823d0,-192490152d0,-101188271d0,   !   2796  1.00000
     . 33068314d0, 145345652d0,  37469931d0,-192563425d0,-102123601d0,   !   2797  1.00000
     . 32950679d0, 145427092d0,  38384035d0,-192636041d0,-103058761d0,   !   2798  1.00000
     . 32833012d0, 145508138d0,  39298131d0,-192707999d0,-103993744d0,   !   2799  1.00000
     . 32715313d0, 145588789d0,  40212214d0,-192779299d0,-104928546d0,   !   2800  1.00000
     . 32597583d0, 145669045d0,  41126279d0,-192849940d0,-105863161d0,   !   2801  1.00000
     . 32479821d0, 145748904d0,  42040322d0,-192919921d0,-106797584d0,   !   2802  1.00000
     . 32362029d0, 145828366d0,  42954338d0,-192989241d0,-107731809d0,   !   2803  1.00000
     . 32244206d0, 145907431d0,  43868324d0,-193057900d0,-108665833d0,   !   2804  1.00000
     . 32126354d0, 145986098d0,  44782274d0,-193125896d0,-109599648d0,   !   2805  1.00000
     . 32008472d0, 146064367d0,  45696183d0,-193193228d0,-110533251d0,   !   2806  1.00000
     . 31890562d0, 146142238d0,  46610049d0,-193259897d0,-111466636d0,   !   2807  1.00000
     . 31772623d0, 146219709d0,  47523865d0,-193325900d0,-112399797d0,   !   2808  1.00000
     . 31654657d0, 146296781d0,  48437628d0,-193391238d0,-113332730d0,   !   2809  1.00000
     . 31536663d0, 146373452d0,  49351333d0,-193455908d0,-114265429d0,   !   2810  1.00000
     . 31418642d0, 146449722d0,  50264976d0,-193519912d0,-115197888d0,   !   2811  1.00000
     . 31300594d0, 146525592d0,  51178551d0,-193583247d0,-116130103d0,   !   2812  1.00000
     . 31182520d0, 146601059d0,  52092055d0,-193645912d0,-117062068d0,   !   2813  1.00000
     . 31064421d0, 146676125d0,  53005482d0,-193707908d0,-117993778d0,   !   2814  1.00000
     . 30946296d0, 146750787d0,  53918829d0,-193769233d0,-118925228d0,   !   2815  1.00000
     . 30828146d0, 146825047d0,  54832090d0,-193829886d0,-119856412d0,   !   2816  1.00000
     . 30709973d0, 146898903d0,  55745262d0,-193889867d0,-120787325d0,   !   2817  1.00000
     . 30591775d0, 146972354d0,  56658339d0,-193949174d0,-121717962d0,   !   2818  1.00000
     . 30473554d0, 147045401d0,  57571318d0,-194007807d0,-122648317d0,   !   2819  1.00000
     . 30355310d0, 147118042d0,  58484192d0,-194065765d0,-123578385d0,   !   2820  1.00000
     . 30237044d0, 147190278d0,  59396959d0,-194123047d0,-124508161d0,   !   2821  1.00000
     . 30118755d0, 147262108d0,  60309613d0,-194179652d0,-125437639d0,   !   2822  1.00000
     . 30000446d0, 147333531d0,  61222149d0,-194235580d0,-126366814d0,   !   2823  1.00000
     . 29882115d0, 147404546d0,  62134564d0,-194290829d0,-127295681d0,   !   2824  1.00000
     . 29763763d0, 147475154d0,  63046852d0,-194345399d0,-128224233d0,   !   2825  1.00000
     . 29645391d0, 147545354d0,  63959009d0,-194399289d0,-129152467d0,   !   2826  1.00000
     . 29527000d0, 147615144d0,  64871031d0,-194452499d0,-130080376d0,   !   2827  1.00000
     . 29408589d0, 147684526d0,  65782911d0,-194505026d0,-131007955d0,   !   2828  1.00000
     . 29290159d0, 147753498d0,  66694647d0,-194556871d0,-131935199d0,   !   2829  1.00000
     . 29171712d0, 147822060d0,  67606234d0,-194608033d0,-132862102d0,   !   2830  1.00000
     . 29053246d0, 147890210d0,  68517666d0,-194658510d0,-133788659d0,   !   2831  1.00000
     . 28934763d0, 147957950d0,  69428939d0,-194708302d0,-134714865d0,   !   2832  1.00000
     . 28816263d0, 148025278d0,  70340048d0,-194757408d0,-135640713d0,   !   2833  1.00000
     . 28697746d0, 148092193d0,  71250990d0,-194805828d0,-136566199d0,   !   2834  1.00000
     . 28579214d0, 148158696d0,  72161758d0,-194853560d0,-137491317d0,   !   2835  1.00000
     . 28460666d0, 148224786d0,  73072349d0,-194900603d0,-138416062d0,   !   2836  1.00000
     . 28342103d0, 148290462d0,  73982758d0,-194946958d0,-139340428d0,   !   2837  1.00000
     . 28223525d0, 148355723d0,  74892980d0,-194992622d0,-140264410d0,   !   2838  1.00000
     . 28104934d0, 148420570d0,  75803010d0,-195037595d0,-141188003d0,   !   2839  1.00000
     . 27986328d0, 148485002d0,  76712844d0,-195081877d0,-142111200d0,   !   2840  1.00000
     . 27867710d0, 148549017d0,  77622477d0,-195125465d0,-143033997d0,   !   2841  1.00000
     . 27749079d0, 148612617d0,  78531904d0,-195168361d0,-143956387d0,   !   2842  1.00000
     . 27630435d0, 148675800d0,  79441121d0,-195210562d0,-144878366d0,   !   2843  1.00000
     . 27511780d0, 148738565d0,  80350123d0,-195252068d0,-145799928d0,   !   2844  1.00000
     . 27393113d0, 148800913d0,  81258905d0,-195292878d0,-146721067d0,   !   2845  1.00000
     . 27274436d0, 148862842d0,  82167463d0,-195332992d0,-147641779d0,   !   2846  1.00000
     . 27155748d0, 148924353d0,  83075791d0,-195372408d0,-148562056d0,   !   2847  1.00000
     . 27037051d0, 148985444d0,  83983886d0,-195411125d0,-149481895d0,   !   2848  1.00000
     . 26918344d0, 149046116d0,  84891741d0,-195449143d0,-150401289d0,   !   2849  1.00000
     . 26799628d0, 149106367d0,  85799353d0,-195486461d0,-151320233d0,   !   2850  1.00000
     . 26680903d0, 149166198d0,  86706717d0,-195523078d0,-152238721d0,   !   2851  1.00000
     . 26562171d0, 149225607d0,  87613828d0,-195558993d0,-153156747d0,   !   2852  1.00000
     . 26443431d0, 149284595d0,  88520681d0,-195594206d0,-154074307d0,   !   2853  1.00000
     . 26324684d0, 149343160d0,  89427272d0,-195628715d0,-154991395d0,   !   2854  1.00000
     . 26205931d0, 149401302d0,  90333595d0,-195662520d0,-155908004d0,   !   2855  1.00000
     . 26087172d0, 149459022d0,  91239646d0,-195695620d0,-156824131d0,   !   2856  1.00000
     . 25968407d0, 149516317d0,  92145421d0,-195728014d0,-157739768d0,   !   2857  1.00000
     . 25849636d0, 149573188d0,  93050913d0,-195759701d0,-158654910d0,   !   2858  1.00000
     . 25730862d0, 149629634d0,  93956120d0,-195790680d0,-159569552d0,   !   2859  1.00000
     . 25612083d0, 149685655d0,  94861035d0,-195820951d0,-160483688d0,   !   2860  1.00000
     . 25493300d0, 149741250d0,  95765653d0,-195850513d0,-161397313d0,   !   2861  1.00000
     . 25374514d0, 149796419d0,  96669971d0,-195879364d0,-162310420d0,   !   2862  1.00000
     . 25255726d0, 149851161d0,  97573984d0,-195907505d0,-163223005d0,   !   2863  1.00000
     . 25136935d0, 149905476d0,  98477686d0,-195934933d0,-164135061d0,   !   2864  1.00000
     . 25018143d0, 149959363d0,  99381072d0,-195961649d0,-165046584d0,   !   2865  1.00000
     . 24899349d0, 150012822d0, 100284138d0,-195987652d0,-165957567d0,   !   2866  1.00000
     . 24780555d0, 150065852d0, 101186879d0,-196012940d0,-166868004d0,   !   2867  1.00000
     . 24661760d0, 150118452d0, 102089291d0,-196037513d0,-167777891d0,   !   2868  1.00000
     . 24542966d0, 150170623d0, 102991367d0,-196061370d0,-168687221d0,   !   2869  1.00000
     . 24424172d0, 150222363d0, 103893104d0,-196084511d0,-169595988d0,   !   2870  1.00000
     . 24305379d0, 150273673d0, 104794497d0,-196106933d0,-170504188d0,   !   2871  1.00000
     . 24186588d0, 150324551d0, 105695540d0,-196128638d0,-171411814d0,   !   2872  1.00000
     . 24067800d0, 150374998d0, 106596229d0,-196149622d0,-172318861d0,   !   2873  1.00000
     . 23949014d0, 150425012d0, 107496559d0,-196169887d0,-173225323d0,   !   2874  1.00000
     . 23830231d0, 150474593d0, 108396524d0,-196189431d0,-174131194d0,   !   2875  1.00000
     . 23711451d0, 150523742d0, 109296121d0,-196208253d0,-175036469d0,   !   2876  1.00000
     . 23592676d0, 150572456d0, 110195344d0,-196226352d0,-175941142d0,   !   2877  1.00000
     . 23473905d0, 150620736d0, 111094189d0,-196243727d0,-176845207d0,   !   2878  1.00000
     . 23355140d0, 150668581d0, 111992650d0,-196260379d0,-177748658d0,   !   2879  1.00000
     . 23236380d0, 150715991d0, 112890722d0,-196276305d0,-178651490d0,   !   2880  1.00000
     . 23117626d0, 150762965d0, 113788401d0,-196291505d0,-179553697d0,   !   2881  1.00000
     . 22998878d0, 150809502d0, 114685682d0,-196305978d0,-180455273d0,   !   2882  1.00000
     . 22880138d0, 150855603d0, 115582559d0,-196319724d0,-181356213d0,   !   2883  1.00000
     . 22761405d0, 150901266d0, 116479028d0,-196332741d0,-182256510d0,   !   2884  1.00000
     . 22642681d0, 150946492d0, 117375084d0,-196345029d0,-183156159d0,   !   2885  1.00000
     . 22523965d0, 150991279d0, 118270721d0,-196356586d0,-184055155d0,   !   2886  1.00000
     . 22405258d0, 151035628d0, 119165936d0,-196367413d0,-184953490d0,   !   2887  1.00000
     . 22286560d0, 151079537d0, 120060723d0,-196377507d0,-185851161d0,   !   2888  1.00000
     . 22167873d0, 151123006d0, 120955076d0,-196386869d0,-186748160d0,   !   2889  1.00000
     . 22049196d0, 151166035d0, 121848992d0,-196395498d0,-187644482d0,   !   2890  1.00000
     . 21930530d0, 151208623d0, 122742464d0,-196403392d0,-188540121d0,   !   2891  1.00000
     . 21811875d0, 151250770d0, 123635489d0,-196410551d0,-189435072d0,   !   2892  1.00000
     . 21693233d0, 151292475d0, 124528060d0,-196416973d0,-190329329d0,   !   2893  1.00000
     . 21574603d0, 151333738d0, 125420174d0,-196422659d0,-191222885d0,   !   2894  1.00000
     . 21455986d0, 151374557d0, 126311825d0,-196427608d0,-192115735d0,   !   2895  1.00000
     . 21337383d0, 151414934d0, 127203008d0,-196431817d0,-193007873d0,   !   2896  1.00000
     . 21218794d0, 151454866d0, 128093717d0,-196435288d0,-193899293d0,   !   2897  1.00000
     . 21100219d0, 151494354d0, 128983949d0,-196438018d0,-194789990d0,   !   2898  1.00000
     . 20981660d0, 151533398d0, 129873698d0,-196440007d0,-195679957d0,   !   2899  1.00000
     . 20863116d0, 151571996d0, 130762959d0,-196441254d0,-196569189d0,   !   2900  1.00000
     . 20744588d0, 151610148d0, 131651727d0,-196441758d0,-197457679d0,   !   2901  1.00000
     . 20626076d0, 151647854d0, 132539996d0,-196441519d0,-198345423d0,   !   2902  1.00000
     . 20507582d0, 151685113d0, 133427763d0,-196440536d0,-199232414d0,   !   2903  1.00000
     . 20389105d0, 151721924d0, 134315021d0,-196438807d0,-200118645d0,   !   2904  1.00000
     . 20270646d0, 151758288d0, 135201766d0,-196436332d0,-201004112d0,   !   2905  1.00000
     . 20152206d0, 151794203d0, 136087992d0,-196433111d0,-201888808d0,   !   2906  1.00000
     . 20033785d0, 151829669d0, 136973695d0,-196429141d0,-202772728d0,   !   2907  1.00000
     . 19915383d0, 151864687d0, 137858870d0,-196424423d0,-203655865d0,   !   2908  1.00000
     . 19797002d0, 151899254d0, 138743511d0,-196418956d0,-204538214d0,   !   2909  1.00000
     . 19678641d0, 151933371d0, 139627614d0,-196412738d0,-205419768d0,   !   2910  1.00000
     . 19560301d0, 151967037d0, 140511172d0,-196405769d0,-206300522d0,   !   2911  1.00000
     . 19441983d0, 152000251d0, 141394182d0,-196398049d0,-207180470d0,   !   2912  1.00000
     . 19323687d0, 152033014d0, 142276638d0,-196389575d0,-208059605d0,   !   2913  1.00000
     . 19205414d0, 152065324d0, 143158535d0,-196380348d0,-208937922d0,   !   2914  1.00000
     . 19087164d0, 152097182d0, 144039868d0,-196370367d0,-209815416d0,   !   2915  1.00000
     . 18968937d0, 152128586d0, 144920631d0,-196359630d0,-210692079d0,   !   2916  1.00000
     . 18850735d0, 152159536d0, 145800820d0,-196348138d0,-211567906d0,   !   2917  1.00000
     . 18732557d0, 152190032d0, 146680430d0,-196335888d0,-212442891d0,   !   2918  1.00000
     . 18614405d0, 152220073d0, 147559454d0,-196322881d0,-213317028d0,   !   2919  1.00000
     . 18496278d0, 152249658d0, 148437889d0,-196309115d0,-214190311d0,   !   2920  1.00000
     . 18378178d0, 152278788d0, 149315729d0,-196294589d0,-215062734d0,   !   2921  1.00000
     . 18260104d0, 152307461d0, 150192969d0,-196279304d0,-215934290d0,   !   2922  1.00000
     . 18142058d0, 152335678d0, 151069604d0,-196263257d0,-216804975d0,   !   2923  1.00000
     . 18024039d0, 152363437d0, 151945628d0,-196246448d0,-217674782d0,   !   2924  1.00000
     . 17906049d0, 152390738d0, 152821036d0,-196228877d0,-218543704d0,   !   2925  1.00000
     . 17788088d0, 152417581d0, 153695824d0,-196210542d0,-219411736d0,   !   2926  1.00000
     . 17670156d0, 152443964d0, 154569985d0,-196191442d0,-220278872d0,   !   2927  1.00000
     . 17552253d0, 152469889d0, 155443516d0,-196171577d0,-221145106d0,   !   2928  1.00000
     . 17434382d0, 152495353d0, 156316410d0,-196150946d0,-222010432d0,   !   2929  1.00000
     . 17316541d0, 152520357d0, 157188662d0,-196129549d0,-222874843d0,   !   2930  1.00000
     . 17198731d0, 152544901d0, 158060268d0,-196107383d0,-223738333d0,   !   2931  1.00000
     . 17080954d0, 152568983d0, 158931222d0,-196084449d0,-224600898d0,   !   2932  1.00000
     . 16963209d0, 152592602d0, 159801518d0,-196060745d0,-225462529d0,   !   2933  1.00000
     . 16845496d0, 152615760d0, 160671152d0,-196036271d0,-226323222d0,   !   2934  1.00000
     . 16727818d0, 152638455d0, 161540119d0,-196011026d0,-227182970d0,   !   2935  1.00000
     . 16610173d0, 152660686d0, 162408412d0,-195985009d0,-228041767d0,   !   2936  1.00000
     . 16492563d0, 152682453d0, 163276027d0,-195958219d0,-228899607d0,   !   2937  1.00000
     . 16374988d0, 152703756d0, 164142959d0,-195930655d0,-229756484d0,   !   2938  1.00000
     . 16257449d0, 152724594d0, 165009203d0,-195902317d0,-230612391d0,   !   2939  1.00000
     . 16139946d0, 152744966d0, 165874752d0,-195873204d0,-231467324d0,   !   2940  1.00000
     . 16022480d0, 152764873d0, 166739602d0,-195843314d0,-232321274d0,   !   2941  1.00000
     . 15905050d0, 152784313d0, 167603748d0,-195812648d0,-233174237d0,   !   2942  1.00000
     . 15787659d0, 152803286d0, 168467184d0,-195781204d0,-234026207d0,   !   2943  1.00000
     . 15670305d0, 152821792d0, 169329905d0,-195748981d0,-234877176d0,   !   2944  1.00000
     . 15552991d0, 152839830d0, 170191906d0,-195715979d0,-235727139d0,   !   2945  1.00000
     . 15435716d0, 152857400d0, 171053182d0,-195682196d0,-236576090d0,   !   2946  1.00000
     . 15318480d0, 152874500d0, 171913726d0,-195647632d0,-237424022d0,   !   2947  1.00000
     . 15201285d0, 152891132d0, 172773535d0,-195612286d0,-238270930d0,   !   2948  1.00000
     . 15084131d0, 152907293d0, 173632602d0,-195576158d0,-239116807d0,   !   2949  1.00000
     . 14967019d0, 152922984d0, 174490922d0,-195539245d0,-239961647d0,   !   2950  1.00000
     . 14849948d0, 152938203d0, 175348490d0,-195501549d0,-240805443d0,   !   2951  1.00000
     . 14732920d0, 152952952d0, 176205300d0,-195463066d0,-241648190d0,   !   2952  1.00000
     . 14615934d0, 152967228d0, 177061348d0,-195423798d0,-242489882d0,   !   2953  1.00000
     . 14498993d0, 152981032d0, 177916628d0,-195383743d0,-243330511d0,   !   2954  1.00000
     . 14382095d0, 152994363d0, 178771135d0,-195342900d0,-244170072d0,   !   2955  1.00000
     . 14265242d0, 153007221d0, 179624862d0,-195301268d0,-245008559d0,   !   2956  1.00000
     . 14148435d0, 153019605d0, 180477806d0,-195258847d0,-245845966d0,   !   2957  1.00000
     . 14031672d0, 153031514d0, 181329960d0,-195215635d0,-246682285d0,   !   2958  1.00000
     . 13914957d0, 153042949d0, 182181319d0,-195171632d0,-247517512d0,   !   2959  1.00000
     . 13798288d0, 153053907d0, 183031878d0,-195126837d0,-248351639d0,   !   2960  1.00000
     . 13681666d0, 153064390d0, 183881631d0,-195081250d0,-249184660d0,   !   2961  1.00000
     . 13565092d0, 153074397d0, 184730574d0,-195034868d0,-250016570d0,   !   2962  1.00000
     . 13448567d0, 153083926d0, 185578700d0,-194987693d0,-250847361d0,   !   2963  1.00000
     . 13332090d0, 153092979d0, 186426004d0,-194939721d0,-251677028d0,   !   2964  1.00000
     . 13215663d0, 153101553d0, 187272481d0,-194890954d0,-252505564d0,   !   2965  1.00000
     . 13099286d0, 153109648d0, 188118126d0,-194841389d0,-253332963d0,   !   2966  1.00000
     . 12982959d0, 153117265d0, 188962932d0,-194791027d0,-254159218d0,   !   2967  1.00000
     . 12866684d0, 153124402d0, 189806896d0,-194739866d0,-254984324d0,   !   2968  1.00000
     . 12750460d0, 153131060d0, 190650010d0,-194687906d0,-255808274d0,   !   2969  1.00000
     . 12634289d0, 153137237d0, 191492270d0,-194635145d0,-256631061d0,   !   2970  1.00000
     . 12518170d0, 153142932d0, 192333671d0,-194581583d0,-257452680d0,   !   2971  1.00000
     . 12402104d0, 153148147d0, 193174206d0,-194527219d0,-258273124d0,   !   2972  1.00000
     . 12286093d0, 153152879d0, 194013871d0,-194472052d0,-259092386d0,   !   2973  1.00000
     . 12170135d0, 153157129d0, 194852660d0,-194416081d0,-259910461d0,   !   2974  1.00000
     . 12054233d0, 153160896d0, 195690567d0,-194359306d0,-260727341d0,   !   2975  1.00000
     . 11938386d0, 153164180d0, 196527588d0,-194301725d0,-261543021d0,   !   2976  1.00000
     . 11822594d0, 153166980d0, 197363716d0,-194243338d0,-262357495d0,   !   2977  1.00000
     . 11706860d0, 153169295d0, 198198946d0,-194184144d0,-263170755d0,   !   2978  1.00000
     . 11591183d0, 153171125d0, 199033273d0,-194124142d0,-263982795d0,   !   2979  1.00000
     . 11475563d0, 153172469d0, 199866692d0,-194063332d0,-264793609d0,   !   2980  1.00000
     . 11360001d0, 153173328d0, 200699195d0,-194001712d0,-265603191d0,   !   2981  1.00000
     . 11244498d0, 153173700d0, 201530780d0,-193939281d0,-266411535d0,   !   2982  1.00000
     . 11129055d0, 153173586d0, 202361438d0,-193876040d0,-267218633d0,   !   2983  1.00000
     . 11013671d0, 153172983d0, 203191166d0,-193811986d0,-268024479d0,   !   2984  1.00000
     . 10898348d0, 153171893d0, 204019958d0,-193747119d0,-268829067d0,   !   2985  1.00000
     . 10783085d0, 153170314d0, 204847808d0,-193681439d0,-269632391d0,   !   2986  1.00000
     . 10667885d0, 153168247d0, 205674710d0,-193614944d0,-270434444d0,   !   2987  1.00000
     . 10552746d0, 153165689d0, 206500660d0,-193547634d0,-271235219d0,   !   2988  1.00000
     . 10437669d0, 153162642d0, 207325651d0,-193479507d0,-272034711d0,   !   2989  1.00000
     . 10322656d0, 153159104d0, 208149679d0,-193410563d0,-272832912d0,   !   2990  1.00000
     . 10207707d0, 153155076d0, 208972736d0,-193340802d0,-273629816d0,   !   2991  1.00000
     . 10092822d0, 153150556d0, 209794819d0,-193270221d0,-274425418d0,   !   2992  1.00000
     .  9978002d0, 153145543d0, 210615921d0,-193198821d0,-275219709d0,   !   2993  1.00000
     .  9863247d0, 153140039d0, 211436037d0,-193126601d0,-276012684d0,   !   2994  1.00000
     .  9748558d0, 153134041d0, 212255161d0,-193053559d0,-276804337d0,   !   2995  1.00000
     .  9633936d0, 153127550d0, 213073288d0,-192979696d0,-277594660d0,   !   2996  1.00000
     .  9519380d0, 153120565d0, 213890412d0,-192905009d0,-278383648d0,   !   2997  1.00000
     .  9404892d0, 153113085d0, 214706527d0,-192829499d0,-279171293d0,   !   2998  1.00000
     .  9290473d0, 153105110d0, 215521629d0,-192753164d0,-279957590d0,   !   2999  1.00000
     .  9176122d0, 153096640d0, 216335711d0,-192676003d0,-280742532d0,   !   3000  1.00000
     .  9061840d0, 153087674d0, 217148767d0,-192598017d0,-281526111d0,   !   3001  1.00000
     .  8947629d0, 153078211d0, 217960793d0,-192519203d0,-282308323d0,   !   3002  1.00000
     .  8833487d0, 153068251d0, 218771782d0,-192439561d0,-283089159d0,   !   3003  1.00000
     .  8719417d0, 153057794d0, 219581729d0,-192359091d0,-283868614d0,   !   3004  1.00000
     .  8605418d0, 153046839d0, 220390628d0,-192277791d0,-284646681d0,   !   3005  1.00000
     .  8491492d0, 153035385d0, 221198474d0,-192195660d0,-285423353d0,   !   3006  1.00000
     .  8377638d0, 153023432d0, 222005261d0,-192112699d0,-286198624d0,   !   3007  1.00000
     .  8263857d0, 153010980d0, 222810983d0,-192028905d0,-286972488d0,   !   3008  1.00000
     .  8150150d0, 152998028d0, 223615635d0,-191944278d0,-287744937d0,   !   3009  1.00000
     .  8036517d0, 152984575d0, 224419211d0,-191858817d0,-288515965d0,   !   3010  1.00000
     .  7922960d0, 152970621d0, 225221705d0,-191772522d0,-289285566d0,   !   3011  1.00000
     .  7809477d0, 152956166d0, 226023112d0,-191685391d0,-290053733d0,   !   3012  1.00000
     .  7696071d0, 152941208d0, 226823427d0,-191597424d0,-290820459d0,   !   3013  1.00000
     .  7582742d0, 152925748d0, 227622642d0,-191508620d0,-291585738d0,   !   3014  1.00000
     .  7469490d0, 152909785d0, 228420753d0,-191418978d0,-292349562d0,   !   3015  1.00000
     .  7356315d0, 152893319d0, 229217755d0,-191328497d0,-293111926d0,   !   3016  1.00000
     .  7243219d0, 152876348d0, 230013640d0,-191237176d0,-293872823d0,   !   3017  1.00000
     .  7130202d0, 152858873d0, 230808404d0,-191145015d0,-294632246d0,   !   3018  1.00000
     .  7017264d0, 152840893d0, 231602042d0,-191052013d0,-295390189d0,   !   3019  1.00000
     .  6904406d0, 152822407d0, 232394546d0,-190958169d0,-296146644d0,   !   3020  1.00000
     .  6791629d0, 152803416d0, 233185912d0,-190863481d0,-296901605d0,   !   3021  1.00000
     .  6678933d0, 152783917d0, 233976133d0,-190767950d0,-297655066d0,   !   3022  1.00000
     .  6566319d0, 152763912d0, 234765205d0,-190671574d0,-298407020d0,   !   3023  1.00000
     .  6453787d0, 152743399d0, 235553121d0,-190574353d0,-299157460d0,   !   3024  1.00000
     .  6341338d0, 152722378d0, 236339875d0,-190476285d0,-299906379d0,   !   3025  1.00000
     .  6228973d0, 152700849d0, 237125463d0,-190377370d0,-300653771d0,   !   3026  1.00000
     .  6116692d0, 152678810d0, 237909877d0,-190277607d0,-301399629d0,   !   3027  1.00000
     .  6004495d0, 152656262d0, 238693113d0,-190176996d0,-302143946d0,   !   3028  1.00000
     .  5892384d0, 152633204d0, 239475164d0,-190075534d0,-302886716d0,   !   3029  1.00000
     .  5780358d0, 152609636d0, 240256025d0,-189973223d0,-303627932d0,   !   3030  1.00000
     .  5668419d0, 152585556d0, 241035690d0,-189870060d0,-304367587d0,   !   3031  1.00000
     .  5556566d0, 152560965d0, 241814153d0,-189766044d0,-305105674d0,   !   3032  1.00000
     .  5444802d0, 152535862d0, 242591409d0,-189661176d0,-305842187d0,   !   3033  1.00000
     .  5333125d0, 152510246d0, 243367452d0,-189555454d0,-306577119d0,   !   3034  1.00000
     .  5221537d0, 152484117d0, 244142275d0,-189448877d0,-307310463d0,   !   3035  1.00000
     .  5110038d0, 152457474d0, 244915873d0,-189341445d0,-308042212d0,   !   3036  1.00000
     .  4998629d0, 152430318d0, 245688241d0,-189233156d0,-308772360d0,   !   3037  1.00000
     .  4887311d0, 152402647d0, 246459372d0,-189124010d0,-309500900d0,   !   3038  1.00000
     .  4776083d0, 152374461d0, 247229261d0,-189014006d0,-310227825d0,   !   3039  1.00000
     .  4664948d0, 152345759d0, 247997901d0,-188903144d0,-310953128d0,   !   3040  1.00000
     .  4553904d0, 152316541d0, 248765288d0,-188791421d0,-311676802d0,   !   3041  1.00000
     .  4442953d0, 152286807d0, 249531414d0,-188678838d0,-312398841d0,   !   3042  1.00000
     .  4332095d0, 152256555d0, 250296275d0,-188565394d0,-313119238d0,   !   3043  1.00000
     .  4221332d0, 152225786d0, 251059865d0,-188451088d0,-313837987d0,   !   3044  1.00000
     .  4110663d0, 152194499d0, 251822177d0,-188335918d0,-314555079d0,   !   3045  1.00000
     .  4000089d0, 152162693d0, 252583206d0,-188219885d0,-315270508d0,   !   3046  1.00000
     .  3889610d0, 152130369d0, 253342945d0,-188102987d0,-315984269d0,   !   3047  1.00000
     .  3779228d0, 152097524d0, 254101390d0,-187985223d0,-316696352d0,   !   3048  1.00000
     .  3668943d0, 152064160d0, 254858534d0,-187866593d0,-317406753d0,   !   3049  1.00000
     .  3558756d0, 152030275d0, 255614371d0,-187747096d0,-318115464d0,   !   3050  1.00000
     .  3448666d0, 151995869d0, 256368896d0,-187626731d0,-318822478d0,   !   3051  1.00000
     .  3338675d0, 151960941d0, 257122102d0,-187505497d0,-319527788d0,   !   3052  1.00000
     .  3228783d0, 151925491d0, 257873983d0,-187383393d0,-320231387d0,   !   3053  1.00000
     .  3118991d0, 151889519d0, 258624534d0,-187260419d0,-320933269d0,   !   3054  1.00000
     .  3009300d0, 151853024d0, 259373749d0,-187136573d0,-321633427d0,   !   3055  1.00000
     .  2899709d0, 151816004d0, 260121622d0,-187011856d0,-322331853d0,   !   3056  1.00000
     .  2790220d0, 151778461d0, 260868147d0,-186886265d0,-323028541d0,   !   3057  1.00000
     .  2680833d0, 151740394d0, 261613318d0,-186759800d0,-323723484d0,   !   3058  1.00000
     .  2571549d0, 151701801d0, 262357129d0,-186632460d0,-324416675d0,   !   3059  1.00000
     .  2462369d0, 151662682d0, 263099574d0,-186504245d0,-325108107d0,   !   3060  1.00000
     .  2353292d0, 151623038d0, 263840647d0,-186375154d0,-325797773d0,   !   3061  1.00000
     .  2244320d0, 151582867d0, 264580343d0,-186245185d0,-326485666d0,   !   3062  1.00000
     .  2135453d0, 151542169d0, 265318655d0,-186114338d0,-327171779d0,   !   3063  1.00000
     .  2026691d0, 151500943d0, 266055577d0,-185982612d0,-327856106d0,   !   3064  1.00000
     .  1918037d0, 151459189d0, 266791103d0,-185850007d0,-328538639d0,   !   3065  1.00000
     .  1809489d0, 151416906d0, 267525228d0,-185716521d0,-329219372d0,   !   3066  1.00000
     .  1701048d0, 151374095d0, 268257946d0,-185582153d0,-329898297d0,   !   3067  1.00000
     .  1592716d0, 151330753d0, 268989250d0,-185446904d0,-330575407d0,   !   3068  1.00000
     .  1484492d0, 151286882d0, 269719134d0,-185310771d0,-331250696d0,   !   3069  1.00000
     .  1376378d0, 151242480d0, 270447593d0,-185173754d0,-331924156d0,   !   3070  1.00000
     .  1268374d0, 151197547d0, 271174620d0,-185035853d0,-332595781d0,   !   3071  1.00000
     .  1160480d0, 151152082d0, 271900210d0,-184897066d0,-333265563d0,   !   3072  1.00000
     .  1052698d0, 151106085d0, 272624356d0,-184757393d0,-333933496d0,   !   3073  1.00000
     .   945027d0, 151059556d0, 273347053d0,-184616832d0,-334599572d0,   !   3074  1.00000
     .   837468d0, 151012493d0, 274068294d0,-184475383d0,-335263784d0,   !   3075  1.00000
     .   730023d0, 150964897d0, 274788073d0,-184333046d0,-335926126d0,   !   3076  1.00000
     .   622691d0, 150916767d0, 275506385d0,-184189818d0,-336586590d0,   !   3077  1.00000
     .   515473d0, 150868102d0, 276223223d0,-184045700d0,-337245170d0,   !   3078  1.00000
     .   408370d0, 150818902d0, 276938582d0,-183900691d0,-337901858d0,   !   3079  1.00000
     .   301382d0, 150769166d0, 277652454d0,-183754789d0,-338556647d0,   !   3080  1.00000
     .   194511d0, 150718894d0, 278364835d0,-183607994d0,-339209530d0,   !   3081  1.00000
     .    87756d0, 150668086d0, 279075718d0,-183460305d0,-339860500d0,   !   3082  1.00000
     .   -18881d0, 150616740d0, 279785097d0,-183311722d0,-340509550d0,   !   3083  1.00000
     .  -125401d0, 150564857d0, 280492966d0,-183162243d0,-341156674d0,   !   3084  1.00000
     .  -231803d0, 150512435d0, 281199319d0,-183011867d0,-341801863d0,   !   3085  1.00000
     .  -338086d0, 150459475d0, 281904149d0,-182860594d0,-342445110d0,   !   3086  1.00000
     .  -444249d0, 150405976d0, 282607451d0,-182708423d0,-343086410d0,   !   3087  1.00000
     .  -550292d0, 150351937d0, 283309219d0,-182555353d0,-343725754d0,   !   3088  1.00000
     .  -656214d0, 150297358d0, 284009447d0,-182401384d0,-344363135d0,   !   3089  1.00000
     .  -762014d0, 150242239d0, 284708128d0,-182246513d0,-344998546d0,   !   3090  1.00000
     .  -867693d0, 150186578d0, 285405256d0,-182090742d0,-345631981d0,   !   3091  1.00000
     .  -973249d0, 150130376d0, 286100825d0,-181934067d0,-346263432d0,   !   3092  1.00000
     . -1078681d0, 150073631d0, 286794829d0,-181776490d0,-346892892d0,   !   3093  1.00000
     . -1183990d0, 150016344d0, 287487263d0,-181618009d0,-347520353d0,   !   3094  1.00000
     . -1289174d0, 149958514d0, 288178119d0,-181458623d0,-348145809d0,   !   3095  1.00000
     . -1394233d0, 149900140d0, 288867391d0,-181298331d0,-348769253d0,   !   3096  1.00000
     . -1499166d0, 149841221d0, 289555075d0,-181137133d0,-349390676d0,   !   3097  1.00000
     . -1603973d0, 149781759d0, 290241162d0,-180975028d0,-350010073d0,   !   3098  1.00000
     . -1708653d0, 149721751d0, 290925648d0,-180812014d0,-350627436d0,   !   3099  1.00000
     . -1813206d0, 149661197d0, 291608526d0,-180648092d0,-351242757d0,   !   3100  1.00000
     . -1917630d0, 149600097d0, 292289789d0,-180483260d0,-351856030d0,   !   3101  1.00000
     . -2021925d0, 149538451d0, 292969433d0,-180317517d0,-352467247d0,   !   3102  1.00000
     . -2126091d0, 149476257d0, 293647449d0,-180150862d0,-353076402d0,   !   3103  1.00000
     . -2230127d0, 149413516d0, 294323834d0,-179983295d0,-353683486d0,   !   3104  1.00000
     . -2334032d0, 149350226d0, 294998579d0,-179814815d0,-354288493d0,   !   3105  1.00000
     . -2437805d0, 149286388d0, 295671679d0,-179645422d0,-354891415d0,   !   3106  1.00000
     . -2541447d0, 149222001d0, 296343128d0,-179475113d0,-355492246d0,   !   3107  1.00000
     . -2644957d0, 149157064d0, 297012919d0,-179303888d0,-356090978d0,   !   3108  1.00000
     . -2748333d0, 149091577d0, 297681047d0,-179131748d0,-356687603d0,   !   3109  1.00000
     . -2851575d0, 149025539d0, 298347505d0,-178958689d0,-357282115d0,   !   3110  1.00000
     . -2954683d0, 148958950d0, 299012287d0,-178784713d0,-357874506d0,   !   3111  1.00000
     . -3057656d0, 148891809d0, 299675386d0,-178609818d0,-358464769d0,   !   3112  1.00000
     . -3160493d0, 148824117d0, 300336796d0,-178434002d0,-359052897d0,   !   3113  1.00000
     . -3263194d0, 148755871d0, 300996512d0,-178257266d0,-359638882d0,   !   3114  1.00000
     . -3365758d0, 148687072d0, 301654527d0,-178079609d0,-360222717d0,   !   3115  1.00000
     . -3468185d0, 148617720d0, 302310834d0,-177901029d0,-360804395d0,   !   3116  1.00000
     . -3570473d0, 148547813d0, 302965428d0,-177721526d0,-361383908d0,   !   3117  1.00000
     . -3672623d0, 148477352d0, 303618302d0,-177541099d0,-361961250d0,   !   3118  1.00000
     . -3774633d0, 148406336d0, 304269450d0,-177359748d0,-362536413d0,   !   3119  1.00000
     . -3876503d0, 148334764d0, 304918865d0,-177177470d0,-363109389d0,   !   3120  1.00000
     . -3978232d0, 148262635d0, 305566541d0,-176994266d0,-363680172d0,   !   3121  1.00000
     . -4079820d0, 148189950d0, 306212473d0,-176810135d0,-364248753d0,   !   3122  1.00000
     . -4181266d0, 148116708d0, 306856653d0,-176625075d0,-364815126d0,   !   3123  1.00000
     . -4282570d0, 148042908d0, 307499076d0,-176439087d0,-365379283d0,   !   3124  1.00000
     . -4383730d0, 147968550d0, 308139734d0,-176252168d0,-365941218d0,   !   3125  1.00000
     . -4484747d0, 147893633d0, 308778623d0,-176064319d0,-366500921d0,   !   3126  1.00000
     . -4585619d0, 147818157d0, 309415735d0,-175875538d0,-367058387d0,   !   3127  1.00000
     . -4686346d0, 147742121d0, 310051064d0,-175685825d0,-367613608d0,   !   3128  1.00000
     . -4786927d0, 147665526d0, 310684604d0,-175495179d0,-368166576d0,   !   3129  1.00000
     . -4887362d0, 147588369d0, 311316348d0,-175303599d0,-368717284d0,   !   3130  1.00000
     . -4987650d0, 147510651d0, 311946291d0,-175111084d0,-369265725d0,   !   3131  1.00000
     . -5087791d0, 147432372d0, 312574425d0,-174917633d0,-369811891d0,   !   3132  1.00000
     . -5187783d0, 147353530d0, 313200745d0,-174723245d0,-370355775d0,   !   3133  1.00000
     . -5287626d0, 147274125d0, 313825244d0,-174527921d0,-370897369d0,   !   3134  1.00000
     . -5387320d0, 147194158d0, 314447916d0,-174331658d0,-371436666d0,   !   3135  1.00000
     . -5486863d0, 147113626d0, 315068754d0,-174134455d0,-371973659d0,   !   3136  1.00000
     . -5586256d0, 147032531d0, 315687751d0,-173936314d0,-372508340d0,   !   3137  1.00000
     . -5685497d0, 146950870d0, 316304903d0,-173737231d0,-373040701d0,   !   3138  1.00000
     . -5784586d0, 146868645d0, 316920202d0,-173537207d0,-373570736d0,   !   3139  1.00000
     . -5883523d0, 146785853d0, 317533642d0,-173336240d0,-374098437d0,   !   3140  1.00000
     . -5982306d0, 146702496d0, 318145216d0,-173134330d0,-374623796d0,   !   3141  1.00000
     . -6080935d0, 146618572d0, 318754918d0,-172931476d0,-375146806d0,   !   3142  1.00000
     . -6179409d0, 146534080d0, 319362742d0,-172727677d0,-375667459d0,   !   3143  1.00000
     . -6277728d0, 146449021d0, 319968681d0,-172522933d0,-376185749d0,   !   3144  1.00000
     . -6375891d0, 146363394d0, 320572728d0,-172317242d0,-376701666d0,   !   3145  1.00000
     . -6473898d0, 146277198d0, 321174878d0,-172110603d0,-377215205d0,   !   3146  1.00000
     . -6571747d0, 146190432d0, 321775125d0,-171903016d0,-377726357d0,   !   3147  1.00000
     . -6669439d0, 146103097d0, 322373460d0,-171694480d0,-378235115d0,   !   3148  1.00000
     . -6766972d0, 146015192d0, 322969879d0,-171484995d0,-378741472d0,   !   3149  1.00000
     . -6864346d0, 145926716d0, 323564374d0,-171274558d0,-379245419d0,   !   3150  1.00000
     . -6961560d0, 145837669d0, 324156940d0,-171063170d0,-379746950d0,   !   3151  1.00000
     . -7058613d0, 145748050d0, 324747569d0,-170850830d0,-380246057d0,   !   3152  1.00000
     . -7155506d0, 145657858d0, 325336255d0,-170637536d0,-380742732d0,   !   3153  1.00000
     . -7252237d0, 145567094d0, 325922992d0,-170423288d0,-381236968d0,   !   3154  1.00000
     . -7348805d0, 145475757d0, 326507774d0,-170208085d0,-381728757d0,   !   3155  1.00000
     . -7445211d0, 145383846d0, 327090593d0,-169991926d0,-382218092d0,   !   3156  1.00000
     . -7541452d0, 145291361d0, 327671444d0,-169774811d0,-382704965d0,   !   3157  1.00000
     . -7637530d0, 145198301d0, 328250319d0,-169556739d0,-383189368d0,   !   3158  1.00000
     . -7733442d0, 145104665d0, 328827213d0,-169337708d0,-383671295d0,   !   3159  1.00000
     . -7829189d0, 145010454d0, 329402119d0,-169117718d0,-384150737d0,   !   3160  1.00000
     . -7924769d0, 144915667d0, 329975030d0,-168896768d0,-384627686d0,   !   3161  1.00000
     . -8020183d0, 144820303d0, 330545940d0,-168674857d0,-385102136d0,   !   3162  1.00000
     . -8115429d0, 144724362d0, 331114842d0,-168451985d0,-385574079d0,   !   3163  1.00000
     . -8210507d0, 144627843d0, 331681730d0,-168228150d0,-386043506d0,   !   3164  1.00000
     . -8305415d0, 144530746d0, 332246597d0,-168003352d0,-386510411d0,   !   3165  1.00000
     . -8400155d0, 144433070d0, 332809437d0,-167777590d0,-386974786d0,   !   3166  1.00000
     . -8494724d0, 144334814d0, 333370243d0,-167550863d0,-387436622d0,   !   3167  1.00000
     . -8589122d0, 144235979d0, 333929008d0,-167323170d0,-387895913d0,   !   3168  1.00000
     . -8683349d0, 144136564d0, 334485727d0,-167094511d0,-388352651d0,   !   3169  1.00000
     . -8777403d0, 144036568d0, 335040392d0,-166864884d0,-388806829d0,   !   3170  1.00000
     . -8871285d0, 143935990d0, 335592997d0,-166634288d0,-389258437d0,   !   3171  1.00000
     . -8964993d0, 143834831d0, 336143535d0,-166402724d0,-389707470d0,   !   3172  1.00000
     . -9058527d0, 143733090d0, 336692001d0,-166170190d0,-390153919d0,   !   3173  1.00000
     . -9151887d0, 143630766d0, 337238386d0,-165936685d0,-390597776d0,   !   3174  1.00000
     . -9245070d0, 143527858d0, 337782685d0,-165702208d0,-391039034d0,   !   3175  1.00000
     . -9338078d0, 143424367d0, 338324891d0,-165466759d0,-391477685d0,   !   3176  1.00000
     . -9430909d0, 143320291d0, 338864997d0,-165230336d0,-391913722d0,   !   3177  1.00000
     . -9523562d0, 143215631d0, 339402998d0,-164992939d0,-392347137d0,   !   3178  1.00000
     . -9616038d0, 143110386d0, 339938885d0,-164754568d0,-392777921d0,   !   3179  1.00000
     . -9708335d0, 143004554d0, 340472653d0,-164515220d0,-393206068d0,   !   3180  1.00000
     . -9800452d0, 142898137d0, 341004295d0,-164274896d0,-393631569d0,   !   3181  1.00000
     . -9892389d0, 142791132d0, 341533804d0,-164033594d0,-394054418d0,   !   3182  1.00000
     . -9984145d0, 142683541d0, 342061174d0,-163791314d0,-394474605d0,   !   3183  1.00000
     .-10075720d0, 142575361d0, 342586398d0,-163548055d0,-394892124d0,   !   3184  1.00000
     .-10167113d0, 142466593d0, 343109470d0,-163303816d0,-395306967d0,   !   3185  1.00000
     .-10258323d0, 142357237d0, 343630382d0,-163058596d0,-395719125d0,   !   3186  1.00000
     .-10349350d0, 142247291d0, 344149128d0,-162812395d0,-396128592d0,   !   3187  1.00000
     .-10440192d0, 142136755d0, 344665702d0,-162565211d0,-396535358d0,   !   3188  1.00000
     .-10530850d0, 142025629d0, 345180096d0,-162317043d0,-396939418d0,   !   3189  1.00000
     .-10621323d0, 141913912d0, 345692305d0,-162067892d0,-397340762d0,   !   3190  1.00000
     .-10711609d0, 141801604d0, 346202321d0,-161817755d0,-397739383d0,   !   3191  1.00000
     .-10801709d0, 141688705d0, 346710137d0,-161566633d0,-398135274d0,   !   3192  1.00000
     .-10891621d0, 141575212d0, 347215748d0,-161314524d0,-398528426d0,   !   3193  1.00000
     .-10981345d0, 141461127d0, 347719146d0,-161061427d0,-398918831d0,   !   3194  1.00000
     .-11070880d0, 141346449d0, 348220325d0,-160807342d0,-399306483d0,   !   3195  1.00000
     .-11160226d0, 141231177d0, 348719278d0,-160552268d0,-399691372d0,   !   3196  1.00000
     .-11249382d0, 141115311d0, 349215997d0,-160296204d0,-400073491d0,   !   3197  1.00000
     .-11338347d0, 140998850d0, 349710478d0,-160039149d0,-400452833d0,   !   3198  1.00000
     .-11427120d0, 140881793d0, 350202712d0,-159781102d0,-400829389d0,   !   3199  1.00000
     .-11515702d0, 140764141d0, 350692693d0,-159522063d0,-401203152d0,   !   3200  1.00000
     .-11604091d0, 140645892d0, 351180415d0,-159262031d0,-401574114d0,   !   3201  1.00000
     .-11692286d0, 140527047d0, 351665870d0,-159001004d0,-401942266d0,   !   3202  1.00000
     .-11780287d0, 140407604d0, 352149052d0,-158738982d0,-402307602d0,   !   3203  1.00000
     .-11868093d0, 140287563d0, 352629954d0,-158475964d0,-402670113d0,   !   3204  1.00000
     .-11955704d0, 140166925d0, 353108569d0,-158211950d0,-403029790d0,   !   3205  1.00000
     .-12043118d0, 140045687d0, 353584891d0,-157946938d0,-403386628d0,   !   3206  1.00000
     .-12130336d0, 139923850d0, 354058912d0,-157680928d0,-403740617d0,   !   3207  1.00000
     .-12217356d0, 139801413d0, 354530627d0,-157413918d0,-404091749d0,   !   3208  1.00000
     .-12304178d0, 139678376d0, 355000028d0,-157145909d0,-404440017d0,   !   3209  1.00000
     .-12390801d0, 139554738d0, 355467108d0,-156876899d0,-404785413d0,   !   3210  1.00000
     .-12477224d0, 139430498d0, 355931861d0,-156606887d0,-405127929d0,   !   3211  1.00000
     .-12563448d0, 139305657d0, 356394280d0,-156335872d0,-405467557d0,   !   3212  1.00000
     .-12649470d0, 139180214d0, 356854358d0,-156063854d0,-405804288d0,   !   3213  1.00000
     .-12735291d0, 139054167d0, 357312087d0,-155790832d0,-406138116d0,   !   3214  1.00000
     .-12820909d0, 138927518d0, 357767463d0,-155516804d0,-406469032d0,   !   3215  1.00000
     .-12906324d0, 138800264d0, 358220477d0,-155241771d0,-406797028d0,   !   3216  1.00000
     .-12991536d0, 138672406d0, 358671123d0,-154965731d0,-407122096d0,   !   3217  1.00000
     .-13076544d0, 138543944d0, 359119393d0,-154688684d0,-407444228d0,   !   3218  1.00000
     .-13161346d0, 138414876d0, 359565282d0,-154410628d0,-407763417d0,   !   3219  1.00000
     .-13245943d0, 138285202d0, 360008782d0,-154131563d0,-408079654d0,   !   3220  1.00000
     .-13330333d0, 138154922d0, 360449887d0,-153851488d0,-408392931d0,   !   3221  1.00000
     .-13414516d0, 138024034d0, 360888589d0,-153570402d0,-408703240d0,   !   3222  1.00000
     .-13498492d0, 137892540d0, 361324881d0,-153288304d0,-409010574d0,   !   3223  1.00000
     .-13582259d0, 137760438d0, 361758758d0,-153005194d0,-409314924d0,   !   3224  1.00000
     .-13665816d0, 137627727d0, 362190211d0,-152721070d0,-409616282d0,   !   3225  1.00000
     .-13749164d0, 137494408d0, 362619235d0,-152435932d0,-409914641d0,   !   3226  1.00000
     .-13832302d0, 137360479d0, 363045821d0,-152149779d0,-410209992d0,   !   3227  1.00000
     .-13915228d0, 137225940d0, 363469964d0,-151862610d0,-410502327d0,   !   3228  1.00000
     .-13997942d0, 137090791d0, 363891657d0,-151574424d0,-410791638d0,   !   3229  1.00000
     .-14080444d0, 136955031d0, 364310891d0,-151285221d0,-411077917d0,   !   3230  1.00000
     .-14162733d0, 136818660d0, 364727662d0,-150994999d0,-411361156d0,   !   3231  1.00000
     .-14244807d0, 136681677d0, 365141961d0,-150703758d0,-411641347d0,   !   3232  1.00000
     .-14326667d0, 136544081d0, 365553782d0,-150411497d0,-411918483d0,   !   3233  1.00000
     .-14408312d0, 136405873d0, 365963118d0,-150118215d0,-412192554d0,   !   3234  1.00000
     .-14489740d0, 136267051d0, 366369961d0,-149823911d0,-412463552d0,   !   3235  1.00000
     .-14570952d0, 136127616d0, 366774306d0,-149528585d0,-412731471d0,   !   3236  1.00000
     .-14651946d0, 135987566d0, 367176145d0,-149232235d0,-412996301d0,   !   3237  1.00000
     .-14732722d0, 135846901d0, 367575470d0,-148934861d0,-413258034d0,   !   3238  1.00000
     .-14813280d0, 135705621d0, 367972276d0,-148636462d0,-413516663d0,   !   3239  1.00000
     .-14893617d0, 135563725d0, 368366555d0,-148337036d0,-413772179d0,   !   3240  1.00000
     .-14973735d0, 135421212d0, 368758301d0,-148036585d0,-414024575d0,   !   3241  1.00000
     .-15053632d0, 135278083d0, 369147505d0,-147735105d0,-414273841d0,   !   3242  1.00000
     .-15133307d0, 135134337d0, 369534162d0,-147432597d0,-414519970d0,   !   3243  1.00000
     .-15212759d0, 134989972d0, 369918264d0,-147129060d0,-414762954d0,   !   3244  1.00000
     .-15291989d0, 134844989d0, 370299805d0,-146824492d0,-415002785d0,   !   3245  1.00000
     .-15370995d0, 134699387d0, 370678776d0,-146518894d0,-415239454d0,   !   3246  1.00000
     .-15449777d0, 134553166d0, 371055173d0,-146212264d0,-415472953d0,   !   3247  1.00000
     .-15528333d0, 134406325d0, 371428986d0,-145904601d0,-415703275d0,   !   3248  1.00000
     .-15606664d0, 134258863d0, 371800210d0,-145595905d0,-415930410d0,   !   3249  1.00000
     .-15684768d0, 134110781d0, 372168837d0,-145286175d0,-416154352d0,   !   3250  1.00000
     .-15762645d0, 133962077d0, 372534860d0,-144975409d0,-416375091d0,   !   3251  1.00000
     .-15840294d0, 133812752d0, 372898273d0,-144663607d0,-416592619d0,   !   3252  1.00000
     .-15917714d0, 133662804d0, 373259068d0,-144350769d0,-416806928d0,   !   3253  1.00000
     .-15994906d0, 133512233d0, 373617238d0,-144036893d0,-417018011d0,   !   3254  1.00000
     .-16071867d0, 133361038d0, 373972776d0,-143721978d0,-417225858d0,   !   3255  1.00000
     .-16148597d0, 133209220d0, 374325675d0,-143406025d0,-417430462d0,   !   3256  1.00000
     .-16225096d0, 133056777d0, 374675928d0,-143089031d0,-417631814d0,   !   3257  1.00000
     .-16301362d0, 132903710d0, 375023529d0,-142770996d0,-417829906d0,   !   3258  1.00000
     .-16377396d0, 132750017d0, 375368469d0,-142451919d0,-418024731d0,   !   3259  1.00000
     .-16453197d0, 132595698d0, 375710742d0,-142131800d0,-418216278d0,   !   3260  1.00000
     .-16528763d0, 132440753d0, 376050341d0,-141810637d0,-418404542d0,   !   3261  1.00000
     .-16604094d0, 132285181d0, 376387258d0,-141488430d0,-418589512d0,   !   3262  1.00000
     .-16679189d0, 132128981d0, 376721487d0,-141165177d0,-418771181d0,   !   3263  1.00000
     .-16754048d0, 131972154d0, 377053021d0,-140840879d0,-418949541d0,   !   3264  1.00000
     .-16828670d0, 131814698d0, 377381852d0,-140515534d0,-419124583d0,   !   3265  1.00000
     .-16903054d0, 131656613d0, 377707973d0,-140189141d0,-419296299d0,   !   3266  1.00000
     .-16977199d0, 131497899d0, 378031378d0,-139861700d0,-419464681d0,   !   3267  1.00000
     .-17051105d0, 131338555d0, 378352059d0,-139533209d0,-419629720d0,   !   3268  1.00000
     .-17124772d0, 131178580d0, 378670008d0,-139203668d0,-419791409d0,   !   3269  1.00000
     .-17198197d0, 131017975d0, 378985220d0,-138873076d0,-419949738d0,   !   3270  1.00000
     .-17271381d0, 130856738d0, 379297686d0,-138541433d0,-420104700d0,   !   3271  1.00000
     .-17344324d0, 130694869d0, 379607400d0,-138208736d0,-420256286d0,   !   3272  1.00000
     .-17417023d0, 130532368d0, 379914354d0,-137874986d0,-420404488d0,   !   3273  1.00000
     .-17489479d0, 130369234d0, 380218541d0,-137540182d0,-420549297d0,   !   3274  1.00000
     .-17561690d0, 130205466d0, 380519955d0,-137204323d0,-420690705d0,   !   3275  1.00000
     .-17633657d0, 130041065d0, 380818587d0,-136867407d0,-420828705d0,   !   3276  1.00000
     .-17705378d0, 129876029d0, 381114432d0,-136529435d0,-420963286d0,   !   3277  1.00000
     .-17776852d0, 129710358d0, 381407481d0,-136190405d0,-421094442d0,   !   3278  1.00000
     .-17848080d0, 129544051d0, 381697727d0,-135850317d0,-421222164d0,   !   3279  1.00000
     .-17919059d0, 129377109d0, 381985164d0,-135509169d0,-421346443d0,   !   3280  1.00000
     .-17989790d0, 129209531d0, 382269783d0,-135166961d0,-421467271d0,   !   3281  1.00000
     .-18060272d0, 129041315d0, 382551579d0,-134823692d0,-421584639d0,   !   3282  1.00000
     .-18130504d0, 128872462d0, 382830543d0,-134479361d0,-421698540d0,   !   3283  1.00000
     .-18200485d0, 128702971d0, 383106669d0,-134133967d0,-421808965d0,   !   3284  1.00000
     .-18270214d0, 128532842d0, 383379949d0,-133787510d0,-421915904d0,   !   3285  1.00000
     .-18339692d0, 128362073d0, 383650375d0,-133439989d0,-422019351d0,   !   3286  1.00000
     .-18408916d0, 128190665d0, 383917942d0,-133091402d0,-422119296d0,   !   3287  1.00000
     .-18477887d0, 128018618d0, 384182642d0,-132741749d0,-422215732d0,   !   3288  1.00000
     .-18546604d0, 127845929d0, 384444466d0,-132391030d0,-422308649d0,   !   3289  1.00000
     .-18615066d0, 127672600d0, 384703409d0,-132039242d0,-422398039d0,   !   3290  1.00000
     .-18683272d0, 127498630d0, 384959463d0,-131686386d0,-422483894d0,   !   3291  1.00000
     .-18751221d0, 127324017d0, 385212620d0,-131332461d0,-422566206d0,   !   3292  1.00000
     .-18818913d0, 127148762d0, 385462874d0,-130977466d0,-422644965d0,   !   3293  1.00000
     .-18886347d0, 126972864d0, 385710216d0,-130621399d0,-422720164d0,   !   3294  1.00000
     .-18953523d0, 126796323d0, 385954641d0,-130264261d0,-422791793d0,   !   3295  1.00000
     .-19020439d0, 126619137d0, 386196140d0,-129906049d0,-422859845d0,   !   3296  1.00000
     .-19087095d0, 126441307d0, 386434706d0,-129546765d0,-422924311d0,   !   3297  1.00000
     .-19153490d0, 126262832d0, 386670332d0,-129186406d0,-422985182d0,   !   3298  1.00000
     .-19219623d0, 126083712d0, 386903011d0,-128824971d0,-423042451d0,   !   3299  1.00000
     .-19285494d0, 125903946d0, 387132736d0,-128462461d0,-423096107d0,   !   3300  1.00000
     .-19351102d0, 125723533d0, 387359498d0,-128098873d0,-423146144d0,   !   3301  1.00000
     .-19416447d0, 125542474d0, 387583291d0,-127734208d0,-423192552d0,   !   3302  1.00000
     .-19481526d0, 125360766d0, 387804108d0,-127368465d0,-423235323d0,   !   3303  1.00000
     .-19546341d0, 125178411d0, 388021941d0,-127001642d0,-423274448d0,   !   3304  1.00000
     .-19610889d0, 124995408d0, 388236783d0,-126633738d0,-423309919d0,   !   3305  1.00000
     .-19675171d0, 124811755d0, 388448626d0,-126264754d0,-423341727d0,   !   3306  1.00000
     .-19739186d0, 124627453d0, 388657463d0,-125894688d0,-423369864d0,   !   3307  1.00000
     .-19802932d0, 124442501d0, 388863287d0,-125523539d0,-423394321d0,   !   3308  1.00000
     .-19866409d0, 124256899d0, 389066091d0,-125151306d0,-423415090d0,   !   3309  1.00000
     .-19929617d0, 124070645d0, 389265866d0,-124777989d0,-423432161d0,   !   3310  1.00000
     .-19992554d0, 123883740d0, 389462607d0,-124403587d0,-423445527d0,   !   3311  1.00000
     .-20055220d0, 123696183d0, 389656304d0,-124028098d0,-423455179d0,   !   3312  1.00000
     .-20117614d0, 123507974d0, 389846952d0,-123651523d0,-423461108d0,   !   3313  1.00000
     .-20179736d0, 123319111d0, 390034543d0,-123273860d0,-423463305d0,   !   3314  1.00000
     .-20241584d0, 123129596d0, 390219068d0,-122895108d0,-423461763d0,   !   3315  1.00000
     .-20303158d0, 122939426d0, 390400522d0,-122515267d0,-423456472d0,   !   3316  1.00000
     .-20364458d0, 122748601d0, 390578896d0,-122134336d0,-423447424d0,   !   3317  1.00000
     .-20425481d0, 122557122d0, 390754183d0,-121752313d0,-423434609d0,   !   3318  1.00000
     .-20486229d0, 122364987d0, 390926376d0,-121369199d0,-423418021d0,   !   3319  1.00000
     .-20546699d0, 122172196d0, 391095467d0,-120984991d0,-423397649d0,   !   3320  1.00000
     .-20606892d0, 121978749d0, 391261449d0,-120599691d0,-423373485d0,   !   3321  1.00000
     .-20666806d0, 121784644d0, 391424314d0,-120213295d0,-423345521d0,   !   3322  1.00000
     .-20726441d0, 121589883d0, 391584055d0,-119825805d0,-423313748d0,   !   3323  1.00000
     .-20785796d0, 121394463d0, 391740665d0,-119437218d0,-423278157d0,   !   3324  1.00000
     .-20844870d0, 121198384d0, 391894136d0,-119047534d0,-423238740d0,   !   3325  1.00000
     .-20903663d0, 121001647d0, 392044460d0,-118656753d0,-423195487d0,   !   3326  1.00000
     .-20962173d0, 120804250d0, 392191631d0,-118264873d0,-423148391d0,   !   3327  1.00000
     .-21020400d0, 120606193d0, 392335641d0,-117871893d0,-423097442d0,   !   3328  1.00000
     .-21078344d0, 120407476d0, 392476482d0,-117477814d0,-423042632d0,   !   3329  1.00000
     .-21136003d0, 120208097d0, 392614146d0,-117082633d0,-422983952d0,   !   3330  1.00000
     .-21193377d0, 120008058d0, 392748628d0,-116686350d0,-422921394d0,   !   3331  1.00000
     .-21250465d0, 119807356d0, 392879918d0,-116288964d0,-422854948d0,   !   3332  1.00000
     .-21307266d0, 119605991d0, 393008009d0,-115890475d0,-422784606d0,   !   3333  1.00000
     .-21363779d0, 119403964d0, 393132894d0,-115490881d0,-422710359d0,   !   3334  1.00000
     .-21420005d0, 119201273d0, 393254566d0,-115090182d0,-422632199d0,   !   3335  1.00000
     .-21475941d0, 118997918d0, 393373017d0,-114688376d0,-422550116d0,   !   3336  1.00000
     .-21531588d0, 118793898d0, 393488239d0,-114285464d0,-422464102d0,   !   3337  1.00000
     .-21586944d0, 118589214d0, 393600225d0,-113881444d0,-422374149d0,   !   3338  1.00000
     .-21642009d0, 118383864d0, 393708967d0,-113476315d0,-422280246d0,   !   3339  1.00000
     .-21696782d0, 118177848d0, 393814458d0,-113070077d0,-422182387d0,   !   3340  1.00000
     .-21751262d0, 117971165d0, 393916691d0,-112662729d0,-422080561d0,   !   3341  1.00000
     .-21805449d0, 117763815d0, 394015657d0,-112254269d0,-421974761d0,   !   3342  1.00000
     .-21859341d0, 117555798d0, 394111350d0,-111844697d0,-421864977d0,   !   3343  1.00000
     .-21912939d0, 117347113d0, 394203761d0,-111434013d0,-421751200d0,   !   3344  1.00000
     .-21966241d0, 117137759d0, 394292883d0,-111022215d0,-421633422d0,   !   3345  1.00000
     .-22019246d0, 116927736d0, 394378709d0,-110609302d0,-421511634d0,   !   3346  1.00000
     .-22071954d0, 116717044d0, 394461231d0,-110195274d0,-421385827d0,   !   3347  1.00000
     .-22124363d0, 116505682d0, 394540442d0,-109780130d0,-421255992d0,   !   3348  1.00000
     .-22176475d0, 116293649d0, 394616333d0,-109363869d0,-421122120d0,   !   3349  1.00000
     .-22228286d0, 116080945d0, 394688898d0,-108946490d0,-420984203d0,   !   3350  1.00000
     .-22279797d0, 115867569d0, 394758128d0,-108527992d0,-420842232d0,   !   3351  1.00000
     .-22331007d0, 115653522d0, 394824017d0,-108108375d0,-420696197d0,   !   3352  1.00000
     .-22381916d0, 115438802d0, 394886556d0,-107687637d0,-420546091d0,   !   3353  1.00000
     .-22432521d0, 115223409d0, 394945738d0,-107265778d0,-420391904d0,   !   3354  1.00000
     .-22482824d0, 115007342d0, 395001555d0,-106842798d0,-420233627d0,   !   3355  1.00000
     .-22532822d0, 114790601d0, 395054000d0,-106418694d0,-420071251d0,   !   3356  1.00000
     .-22582515d0, 114573186d0, 395103065d0,-105993467d0,-419904768d0,   !   3357  1.00000
     .-22631903d0, 114355096d0, 395148743d0,-105567115d0,-419734168d0,   !   3358  1.00000
     .-22680984d0, 114136330d0, 395191025d0,-105139638d0,-419559443d0,   !   3359  1.00000
     .-22729759d0, 113916888d0, 395229905d0,-104711034d0,-419380584d0,   !   3360  1.00000
     .-22778225d0, 113696769d0, 395265374d0,-104281304d0,-419197582d0,   !   3361  1.00000
     .-22826382d0, 113475974d0, 395297425d0,-103850446d0,-419010428d0,   !   3362  1.00000
     .-22874230d0, 113254501d0, 395326050d0,-103418459d0,-418819112d0,   !   3363  1.00000
     .-22921768d0, 113032350d0, 395351242d0,-102985342d0,-418623627d0,   !   3364  1.00000
     .-22968995d0, 112809520d0, 395372993d0,-102551095d0,-418423963d0,   !   3365  1.00000
     .-23015910d0, 112586011d0, 395391295d0,-102115717d0,-418220111d0,   !   3366  1.00000
     .-23062513d0, 112361823d0, 395406140d0,-101679207d0,-418012063d0,   !   3367  1.00000
     .-23108802d0, 112136954d0, 395417522d0,-101241564d0,-417799809d0,   !   3368  1.00000
     .-23154777d0, 111911405d0, 395425431d0,-100802788d0,-417583340d0,   !   3369  1.00000
     .-23200437d0, 111685175d0, 395429862d0,-100362876d0,-417362648d0,   !   3370  1.00000
     .-23245782d0, 111458264d0, 395430805d0, -99921830d0,-417137722d0,   !   3371  1.00000
     .-23290810d0, 111230670d0, 395428253d0, -99479647d0,-416908556d0,   !   3372  1.00000
     .-23335522d0, 111002394d0, 395422198d0, -99036327d0,-416675139d0,   !   3373  1.00000
     .-23379915d0, 110773434d0, 395412633d0, -98591869d0,-416437462d0,   !   3374  1.00000
     .-23423989d0, 110543791d0, 395399550d0, -98146272d0,-416195516d0,   !   3375  1.00000
     .-23467744d0, 110313464d0, 395382942d0, -97699536d0,-415949293d0,   !   3376  1.00000
     .-23511179d0, 110082453d0, 395362799d0, -97251660d0,-415698784d0,   !   3377  1.00000
     .-23554293d0, 109850756d0, 395339116d0, -96802642d0,-415443978d0,   !   3378  1.00000
     .-23597085d0, 109618374d0, 395311884d0, -96352482d0,-415184868d0,   !   3379  1.00000
     .-23639555d0, 109385305d0, 395281094d0, -95901180d0,-414921445d0,   !   3380  1.00000
     .-23681701d0, 109151550d0, 395246741d0, -95448733d0,-414653698d0,   !   3381  1.00000
     .-23723523d0, 108917108d0, 395208815d0, -94995142d0,-414381620d0,   !   3382  1.00000
     .-23765020d0, 108681979d0, 395167309d0, -94540406d0,-414105201d0,   !   3383  1.00000
     .-23806191d0, 108446161d0, 395122215d0, -94084523d0,-413824432d0,   !   3384  1.00000
     .-23847036d0, 108209654d0, 395073526d0, -93627493d0,-413539304d0,   !   3385  1.00000
     .-23887554d0, 107972459d0, 395021233d0, -93169316d0,-413249808d0,   !   3386  1.00000
     .-23927744d0, 107734574d0, 394965329d0, -92709989d0,-412955935d0,   !   3387  1.00000
     .-23967605d0, 107495999d0, 394905806d0, -92249513d0,-412657676d0,   !   3388  1.00000
     .-24007136d0, 107256733d0, 394842657d0, -91787887d0,-412355021d0,   !   3389  1.00000
     .-24046337d0, 107016776d0, 394775872d0, -91325109d0,-412047961d0,   !   3390  1.00000
     .-24085207d0, 106776128d0, 394705446d0, -90861179d0,-411736489d0,   !   3391  1.00000
     .-24123745d0, 106534787d0, 394631369d0, -90396096d0,-411420593d0,   !   3392  1.00000
     .-24161950d0, 106292754d0, 394553634d0, -89929860d0,-411100266d0,   !   3393  1.00000
     .-24199822d0, 106050028d0, 394472234d0, -89462469d0,-410775498d0,   !   3394  1.00000
     .-24237359d0, 105806608d0, 394387159d0, -88993922d0,-410446279d0,   !   3395  1.00000
     .-24274562d0, 105562494d0, 394298404d0, -88524220d0,-410112602d0,   !   3396  1.00000
     .-24311428d0, 105317685d0, 394205958d0, -88053360d0,-409774456d0,   !   3397  1.00000
     .-24347958d0, 105072182d0, 394109816d0, -87581342d0,-409431833d0,   !   3398  1.00000
     .-24384151d0, 104825982d0, 394009969d0, -87108165d0,-409084723d0,   !   3399  1.00000
     .-24420005d0, 104579087d0, 393906408d0, -86633829d0,-408733117d0,   !   3400  1.00000
     .-24455521d0, 104331495d0, 393799127d0, -86158332d0,-408377006d0,   !   3401  1.00000
     .-24490696d0, 104083206d0, 393688117d0, -85681674d0,-408016381d0,   !   3402  1.00000
     .-24525531d0, 103834219d0, 393573370d0, -85203854d0,-407651232d0,   !   3403  1.00000
     .-24560025d0, 103584534d0, 393454879d0, -84724871d0,-407281551d0,   !   3404  1.00000
     .-24594177d0, 103334151d0, 393332636d0, -84244724d0,-406907328d0,   !   3405  1.00000
     .-24627985d0, 103083068d0, 393206633d0, -83763413d0,-406528554d0,   !   3406  1.00000
     .-24661450d0, 102831286d0, 393076861d0, -83280936d0,-406145220d0,   !   3407  1.00000
     .-24694571d0, 102578803d0, 392943313d0, -82797293d0,-405757316d0,   !   3408  1.00000
     .-24727346d0, 102325620d0, 392805981d0, -82312482d0,-405364834d0,   !   3409  1.00000
     .-24759775d0, 102071736d0, 392664858d0, -81826504d0,-404967764d0,   !   3410  1.00000
     .-24791857d0, 101817150d0, 392519934d0, -81339357d0,-404566096d0,   !   3411  1.00000
     .-24823591d0, 101561862d0, 392371203d0, -80851040d0,-404159822d0,   !   3412  1.00000
     .-24854977d0, 101305871d0, 392218656d0, -80361553d0,-403748933d0,   !   3413  1.00000
     .-24886014d0, 101049178d0, 392062286d0, -79870894d0,-403333419d0,   !   3414  1.00000
     .-24916701d0, 100791780d0, 391902083d0, -79379064d0,-402913270d0,   !   3415  1.00000
     .-24947036d0, 100533679d0, 391738042d0, -78886060d0,-402488478d0,   !   3416  1.00000
     .-24977021d0, 100274873d0, 391570153d0, -78391883d0,-402059033d0,   !   3417  1.00000
     .-25006652d0, 100015361d0, 391398408d0, -77896530d0,-401624926d0,   !   3418  1.00000
     .-25035931d0,  99755144d0, 391222800d0, -77400003d0,-401186148d0,   !   3419  1.00000
     .-25064855d0,  99494221d0, 391043320d0, -76902299d0,-400742689d0,   !   3420  1.00000
     .-25093424d0,  99232592d0, 390859961d0, -76403418d0,-400294541d0,   !   3421  1.00000
     .-25121638d0,  98970255d0, 390672715d0, -75903359d0,-399841693d0,   !   3422  1.00000
     .-25149496d0,  98707210d0, 390481573d0, -75402121d0,-399384136d0,   !   3423  1.00000
     .-25176996d0,  98443457d0, 390286528d0, -74899703d0,-398921861d0,   !   3424  1.00000
     .-25204138d0,  98178996d0, 390087572d0, -74396105d0,-398454860d0,   !   3425  1.00000
     .-25230921d0,  97913825d0, 389884696d0, -73891326d0,-397983121d0,   !   3426  1.00000
     .-25257345d0,  97647945d0, 389677893d0, -73385364d0,-397506637d0,   !   3427  1.00000
     .-25283408d0,  97381355d0, 389467154d0, -72878219d0,-397025397d0,   !   3428  1.00000
     .-25309110d0,  97114054d0, 389252472d0, -72369891d0,-396539393d0,   !   3429  1.00000
     .-25334450d0,  96846042d0, 389033838d0, -71860378d0,-396048615d0,   !   3430  1.00000
     .-25359427d0,  96577318d0, 388811245d0, -71349679d0,-395553053d0,   !   3431  1.00000
     .-25384040d0,  96307882d0, 388584685d0, -70837794d0,-395052699d0,   !   3432  1.00000
     .-25408289d0,  96037733d0, 388354148d0, -70324722d0,-394547542d0,   !   3433  1.00000
     .-25432172d0,  95766871d0, 388119628d0, -69810462d0,-394037574d0,   !   3434  1.00000
     .-25455690d0,  95495296d0, 387881116d0, -69295013d0,-393522785d0,   !   3435  1.00000
     .-25478840d0,  95223006d0, 387638605d0, -68778374d0,-393003166d0,   !   3436  1.00000
     .-25501623d0,  94950002d0, 387392085d0, -68260545d0,-392478706d0,   !   3437  1.00000
     .-25524037d0,  94676282d0, 387141550d0, -67741524d0,-391949398d0,   !   3438  1.00000
     .-25546082d0,  94401847d0, 386886990d0, -67221312d0,-391415231d0,   !   3439  1.00000
     .-25567757d0,  94126696d0, 386628398d0, -66699906d0,-390876195d0,   !   3440  1.00000
     .-25589061d0,  93850828d0, 386365766d0, -66177306d0,-390332283d0,   !   3441  1.00000
     .-25609993d0,  93574242d0, 386099086d0, -65653512d0,-389783483d0,   !   3442  1.00000
     .-25630553d0,  93296940d0, 385828349d0, -65128522d0,-389229787d0,   !   3443  1.00000
     .-25650739d0,  93018918d0, 385553548d0, -64602335d0,-388671185d0,   !   3444  1.00000
     .-25670551d0,  92740179d0, 385274674d0, -64074952d0,-388107668d0,   !   3445  1.00000
     .-25689988d0,  92460720d0, 384991719d0, -63546370d0,-387539226d0,   !   3446  1.00000
     .-25709049d0,  92180541d0, 384704675d0, -63016590d0,-386965849d0,   !   3447  1.00000
     .-25727734d0,  91899642d0, 384413534d0, -62485610d0,-386387529d0,   !   3448  1.00000
     .-25746041d0,  91618023d0, 384118288d0, -61953429d0,-385804256d0,   !   3449  1.00000
     .-25763970d0,  91335682d0, 383818929d0, -61420047d0,-385216020d0,   !   3450  1.00000
     .-25781520d0,  91052619d0, 383515448d0, -60885463d0,-384622811d0,   !   3451  1.00000
     .-25798689d0,  90768835d0, 383207837d0, -60349675d0,-384024621d0,   !   3452  1.00000
     .-25815479d0,  90484327d0, 382896089d0, -59812684d0,-383421439d0,   !   3453  1.00000
     .-25831886d0,  90199097d0, 382580194d0, -59274488d0,-382813257d0,   !   3454  1.00000
     .-25847912d0,  89913143d0, 382260145d0, -58735086d0,-382200064d0,   !   3455  1.00000
     .-25863554d0,  89626464d0, 381935934d0, -58194478d0,-381581851d0,   !   3456  1.00000
     .-25878812d0,  89339061d0, 381607553d0, -57652663d0,-380958609d0,   !   3457  1.00000
     .-25893685d0,  89050933d0, 381274992d0, -57109640d0,-380330328d0,   !   3458  1.00000
     .-25908173d0,  88762079d0, 380938245d0, -56565408d0,-379696998d0,   !   3459  1.00000
     .-25922275d0,  88472498d0, 380597303d0, -56019966d0,-379058610d0,   !   3460  1.00000
     .-25935989d0,  88182191d0, 380252157d0, -55473313d0,-378415154d0,   !   3461  1.00000
     .-25949315d0,  87891157d0, 379902799d0, -54925450d0,-377766621d0,   !   3462  1.00000
     .-25962252d0,  87599395d0, 379549222d0, -54376373d0,-377113002d0,   !   3463  1.00000
     .-25974799d0,  87306905d0, 379191417d0, -53826084d0,-376454285d0,   !   3464  1.00000
     .-25986956d0,  87013686d0, 378829375d0, -53274581d0,-375790463d0,   !   3465  1.00000
     .-25998721d0,  86719738d0, 378463089d0, -52721863d0,-375121525d0,   !   3466  1.00000
     .-26010094d0,  86425060d0, 378092551d0, -52167930d0,-374447462d0,   !   3467  1.00000
     .-26021075d0,  86129652d0, 377717751d0, -51612780d0,-373768263d0,   !   3468  1.00000
     .-26031661d0,  85833514d0, 377338682d0, -51056413d0,-373083921d0,   !   3469  1.00000
     .-26041853d0,  85536644d0, 376955336d0, -50498827d0,-372394424d0,   !   3470  1.00000
     .-26051649d0,  85239042d0, 376567703d0, -49940023d0,-371699763d0,   !   3471  1.00000
     .-26061048d0,  84940708d0, 376175777d0, -49379999d0,-370999929d0,   !   3472  1.00000
     .-26070051d0,  84641641d0, 375779549d0, -48818755d0,-370294911d0,   !   3473  1.00000
     .-26078656d0,  84341842d0, 375379009d0, -48256289d0,-369584701d0,   !   3474  1.00000
     .-26086862d0,  84041308d0, 374974151d0, -47692600d0,-368869289d0,   !   3475  1.00000
     .-26094668d0,  83740040d0, 374564966d0, -47127689d0,-368148664d0,   !   3476  1.00000
     .-26102073d0,  83438038d0, 374151445d0, -46561554d0,-367422818d0,   !   3477  1.00000
     .-26109078d0,  83135300d0, 373733581d0, -45994194d0,-366691740d0,   !   3478  1.00000
     .-26115680d0,  82831826d0, 373311364d0, -45425608d0,-365955421d0,   !   3479  1.00000
     .-26121879d0,  82527617d0, 372884787d0, -44855796d0,-365213851d0,   !   3480  1.00000
     .-26127675d0,  82222670d0, 372453841d0, -44284756d0,-364467020d0,   !   3481  1.00000
     .-26133066d0,  81916986d0, 372018518d0, -43712489d0,-363714920d0,   !   3482  1.00000
     .-26138051d0,  81610565d0, 371578810d0, -43138992d0,-362957539d0,   !   3483  1.00000
     .-26142630d0,  81303405d0, 371134708d0, -42564266d0,-362194868d0,   !   3484  1.00000
     .-26146802d0,  80995507d0, 370686204d0, -41988309d0,-361426898d0,   !   3485  1.00000
     .-26150566d0,  80686869d0, 370233289d0, -41411120d0,-360653619d0,   !   3486  1.00000
     .-26153921d0,  80377491d0, 369775956d0, -40832700d0,-359875021d0,   !   3487  1.00000
     .-26156866d0,  80067374d0, 369314195d0, -40253046d0,-359091095d0,   !   3488  1.00000
     .-26159401d0,  79756515d0, 368847999d0, -39672158d0,-358301829d0,   !   3489  1.00000
     .-26161525d0,  79444915d0, 368377359d0, -39090036d0,-357507216d0,   !   3490  1.00000
     .-26163236d0,  79132574d0, 367902267d0, -38506677d0,-356707245d0,   !   3491  1.00000
     .-26164535d0,  78819490d0, 367422714d0, -37922083d0,-355901906d0,   !   3492  1.00000
     .-26165419d0,  78505664d0, 366938692d0, -37336251d0,-355091189d0,   !   3493  1.00000
     .-26165889d0,  78191094d0, 366450193d0, -36749180d0,-354275085d0,   !   3494  1.00000
     .-26165943d0,  77875780d0, 365957207d0, -36160871d0,-353453584d0,   !   3495  1.00000
     .-26165581d0,  77559722d0, 365459728d0, -35571322d0,-352626676d0,   !   3496  1.00000
     .-26164802d0,  77242920d0, 364957746d0, -34980533d0,-351794351d0,   !   3497  1.00000
     .-26163605d0,  76925372d0, 364451252d0, -34388501d0,-350956599d0,   !   3498  1.00000
     .-26161989d0,  76607078d0, 363940239d0, -33795228d0,-350113412d0,   !   3499  1.00000
     .-26159953d0,  76288038d0, 363424699d0, -33200711d0,-349264778d0,   !   3500  1.00000
     .-26157496d0,  75968252d0, 362904621d0, -32604951d0,-348410688d0,   !   3501  1.00000
     .-26154619d0,  75647718d0, 362379999d0, -32007945d0,-347551132d0,   !   3502  1.00000
     .-26151318d0,  75326436d0, 361850824d0, -31409694d0,-346686100d0,   !   3503  1.00000
     .-26147595d0,  75004406d0, 361317087d0, -30810197d0,-345815583d0,   !   3504  1.00000
     .-26143448d0,  74681627d0, 360778780d0, -30209452d0,-344939570d0,   !   3505  1.00000
     .-26138876d0,  74358099d0, 360235894d0, -29607459d0,-344058052d0,   !   3506  1.00000
     .-26133878d0,  74033822d0, 359688421d0, -29004217d0,-343171018d0,   !   3507  1.00000
     .-26128454d0,  73708793d0, 359136352d0, -28399725d0,-342278460d0,   !   3508  1.00000
     .-26122602d0,  73383015d0, 358579679d0, -27793982d0,-341380367d0,   !   3509  1.00000
     .-26116322d0,  73056484d0, 358018394d0, -27186988d0,-340476728d0,   !   3510  1.00000
     .-26109613d0,  72729202d0, 357452487d0, -26578741d0,-339567535d0,   !   3511  1.00000
     .-26102475d0,  72401168d0, 356881951d0, -25969242d0,-338652777d0,   !   3512  1.00000
     .-26094905d0,  72072381d0, 356306777d0, -25358488d0,-337732445d0,   !   3513  1.00000
     .-26086904d0,  71742841d0, 355726957d0, -24746480d0,-336806528d0,   !   3514  1.00000
     .-26078471d0,  71412547d0, 355142481d0, -24133215d0,-335875016d0,   !   3515  1.00000
     .-26069604d0,  71081498d0, 354553342d0, -23518695d0,-334937900d0,   !   3516  1.00000
     .-26060303d0,  70749695d0, 353959530d0, -22902916d0,-333995170d0,   !   3517  1.00000
     .-26050567d0,  70417137d0, 353361038d0, -22285880d0,-333046816d0,   !   3518  1.00000
     .-26040395d0,  70083822d0, 352757857d0, -21667585d0,-332092827d0,   !   3519  1.00000
     .-26029786d0,  69749752d0, 352149977d0, -21048030d0,-331133194d0,   !   3520  1.00000
     .-26018740d0,  69414924d0, 351537392d0, -20427214d0,-330167906d0,   !   3521  1.00000
     .-26007255d0,  69079340d0, 350920092d0, -19805136d0,-329196955d0,   !   3522  1.00000
     .-25995331d0,  68742997d0, 350298068d0, -19181796d0,-328220329d0,   !   3523  1.00000
     .-25982967d0,  68405896d0, 349671312d0, -18557193d0,-327238020d0,   !   3524  1.00000
     .-25970162d0,  68068037d0, 349039816d0, -17931326d0,-326250016d0,   !   3525  1.00000
     .-25956915d0,  67729418d0, 348403570d0, -17304194d0,-325256308d0,   !   3526  1.00000
     .-25943225d0,  67390039d0, 347762567d0, -16675796d0,-324256886d0,   !   3527  1.00000
     .-25929092d0,  67049900d0, 347116798d0, -16046132d0,-323251739d0,   !   3528  1.00000
     .-25914514d0,  66709000d0, 346466253d0, -15415201d0,-322240859d0,   !   3529  1.00000
     .-25899491d0,  66367339d0, 345810925d0, -14783001d0,-321224234d0,   !   3530  1.00000
     .-25884021d0,  66024916d0, 345150805d0, -14149532d0,-320201855d0,   !   3531  1.00000
     .-25868105d0,  65681730d0, 344485884d0, -13514793d0,-319173712d0,   !   3532  1.00000
     .-25851741d0,  65337782d0, 343816154d0, -12878784d0,-318139794d0,   !   3533  1.00000
     .-25834928d0,  64993071d0, 343141606d0, -12241502d0,-317100092d0,   !   3534  1.00000
     .-25817665d0,  64647596d0, 342462232d0, -11602949d0,-316054596d0,   !   3535  1.00000
     .-25799952d0,  64301356d0, 341778022d0, -10963122d0,-315003295d0,   !   3536  1.00000
     .-25781788d0,  63954351d0, 341088968d0, -10322021d0,-313946179d0,   !   3537  1.00000
     .-25763171d0,  63606582d0, 340395061d0,  -9679646d0,-312883239d0,   !   3538  1.00000
     .-25744102d0,  63258046d0, 339696294d0,  -9035994d0,-311814464d0,   !   3539  1.00000
     .-25724578d0,  62908744d0, 338992657d0,  -8391066d0,-310739844d0,   !   3540  1.00000
     .-25704600d0,  62558675d0, 338284141d0,  -7744861d0,-309659369d0,   !   3541  1.00000
     .-25684166d0,  62207839d0, 337570738d0,  -7097377d0,-308573029d0,   !   3542  1.00000
     .-25663275d0,  61856235d0, 336852439d0,  -6448614d0,-307480814d0,   !   3543  1.00000
     .-25641927d0,  61503863d0, 336129236d0,  -5798572d0,-306382713d0,   !   3544  1.00000
     .-25620121d0,  61150722d0, 335401119d0,  -5147248d0,-305278717d0,   !   3545  1.00000
     .-25597855d0,  60796812d0, 334668081d0,  -4494643d0,-304168816d0,   !   3546  1.00000
     .-25575130d0,  60442131d0, 333930112d0,  -3840756d0,-303052998d0,   !   3547  1.00000
     .-25551943d0,  60086681d0, 333187203d0,  -3185585d0,-301931255d0,   !   3548  1.00000
     .-25528295d0,  59730459d0, 332439347d0,  -2529130d0,-300803576d0,   !   3549  1.00000
     .-25504184d0,  59373466d0, 331686534d0,  -1871390d0,-299669951d0,   !   3550  1.00000
     .-25479610d0,  59015702d0, 330928755d0,  -1212365d0,-298530369d0,   !   3551  1.00000
     .-25454571d0,  58657165d0, 330166003d0,   -552053d0,-297384821d0,   !   3552  1.00000
     .-25429067d0,  58297855d0, 329398267d0,    109546d0,-296233296d0,   !   3553  1.00000
     .-25403097d0,  57937771d0, 328625540d0,    772434d0,-295075784d0,   !   3554  1.00000
     .-25376660d0,  57576914d0, 327847813d0,   1436610d0,-293912275d0,   !   3555  1.00000
     .-25349755d0,  57215283d0, 327065076d0,   2102077d0,-292742759d0,   !   3556  1.00000
     .-25322381d0,  56852877d0, 326277321d0,   2768835d0,-291567225d0,   !   3557  1.00000
     .-25294538d0,  56489695d0, 325484540d0,   3436884d0,-290385663d0,   !   3558  1.00000
     .-25266224d0,  56125737d0, 324686724d0,   4106226d0,-289198064d0,   !   3559  1.00000
     .-25237439d0,  55761004d0, 323883863d0,   4776861d0,-288004416d0,   !   3560  1.00000
     .-25208182d0,  55395493d0, 323075949d0,   5448790d0,-286804710d0,   !   3561  1.00000
     .-25178451d0,  55029205d0, 322262974d0,   6122014d0,-285598936d0,   !   3562  1.00000
     .-25148246d0,  54662139d0, 321444928d0,   6796534d0,-284387082d0,   !   3563  1.00000
     .-25117567d0,  54294294d0, 320621803d0,   7472351d0,-283169140d0,   !   3564  1.00000
     .-25086411d0,  53925671d0, 319793589d0,   8149465d0,-281945097d0,   !   3565  1.00000
     .-25054779d0,  53556269d0, 318960279d0,   8827877d0,-280714946d0,   !   3566  1.00000
     .-25022670d0,  53186086d0, 318121863d0,   9507588d0,-279478674d0,   !   3567  1.00000
     .-24990082d0,  52815123d0, 317278333d0,  10188600d0,-278236272d0,   !   3568  1.00000
     .-24957014d0,  52443380d0, 316429679d0,  10870912d0,-276987730d0,   !   3569  1.00000
     .-24923467d0,  52070855d0, 315575893d0,  11554526d0,-275733037d0,   !   3570  1.00000
     .-24889438d0,  51697548d0, 314716966d0,  12239442d0,-274472182d0,   !   3571  1.00000
     .-24854928d0,  51323458d0, 313852889d0,  12925661d0,-273205157d0,   !   3572  1.00000
     .-24819935d0,  50948586d0, 312983653d0,  13613185d0,-271931949d0,   !   3573  1.00000
     .-24784457d0,  50572931d0, 312109250d0,  14302014d0,-270652550d0,   !   3574  1.00000
     .-24748496d0,  50196491d0, 311229670d0,  14992149d0,-269366948d0,   !   3575  1.00000
     .-24712049d0,  49819267d0, 310344905d0,  15683590d0,-268075133d0,   !   3576  1.00000
     .-24675115d0,  49441258d0, 309454946d0,  16376338d0,-266777096d0,   !   3577  1.00000
     .-24637694d0,  49062464d0, 308559784d0,  17070395d0,-265472825d0,   !   3578  1.00000
     .-24599785d0,  48682884d0, 307659410d0,  17765762d0,-264162310d0,   !   3579  1.00000
     .-24561387d0,  48302518d0, 306753816d0,  18462438d0,-262845541d0,   !   3580  1.00000
     .-24522499d0,  47921365d0, 305842991d0,  19160425d0,-261522508d0,   !   3581  1.00000
     .-24483120d0,  47539424d0, 304926928d0,  19859723d0,-260193200d0,   !   3582  1.00000
     .-24443249d0,  47156695d0, 304005618d0,  20560335d0,-258857607d0,   !   3583  1.00000
     .-24402886d0,  46773178d0, 303079051d0,  21262259d0,-257515718d0,   !   3584  1.00000
     .-24362029d0,  46388872d0, 302147219d0,  21965498d0,-256167523d0,   !   3585  1.00000
     .-24320678d0,  46003777d0, 301210113d0,  22670052d0,-254813011d0,   !   3586  1.00000
     .-24278832d0,  45617892d0, 300267724d0,  23375921d0,-253452173d0,   !   3587  1.00000
     .-24236489d0,  45231216d0, 299320042d0,  24083107d0,-252084998d0,   !   3588  1.00000
     .-24193650d0,  44843750d0, 298367060d0,  24791611d0,-250711475d0,   !   3589  1.00000
     .-24150312d0,  44455492d0, 297408768d0,  25501433d0,-249331594d0,   !   3590  1.00000
     .-24106475d0,  44066442d0, 296445157d0,  26212574d0,-247945344d0,   !   3591  1.00000
     .-24062139d0,  43676600d0, 295476218d0,  26925036d0,-246552716d0,   !   3592  1.00000
     .-24017302d0,  43285964d0, 294501942d0,  27638818d0,-245153698d0,   !   3593  1.00000
     .-23971964d0,  42894536d0, 293522320d0,  28353922d0,-243748280d0,   !   3594  1.00000
     .-23926123d0,  42502313d0, 292537344d0,  29070348d0,-242336453d0,   !   3595  1.00000
     .-23879778d0,  42109296d0, 291547004d0,  29788098d0,-240918204d0,   !   3596  1.00000
     .-23832930d0,  41715485d0, 290551292d0,  30507172d0,-239493524d0,   !   3597  1.00000
     .-23785576d0,  41320877d0, 289550197d0,  31227571d0,-238062402d0,   !   3598  1.00000
     .-23737716d0,  40925474d0, 288543712d0,  31949295d0,-236624828d0,   !   3599  1.00000
     .-23689349d0,  40529275d0, 287531828d0,  32672347d0,-235180792d0,   !   3600  1.00000
     .-23640475d0,  40132278d0, 286514535d0,  33396726d0,-233730282d0,   !   3601  1.00000
     .-23591091d0,  39734484d0, 285491824d0,  34122433d0,-232273289d0,   !   3602  1.00000
     .-23541198d0,  39335892d0, 284463686d0,  34849470d0,-230809801d0,   !   3603  1.00000
     .-23490795d0,  38936502d0, 283430113d0,  35577837d0,-229339809d0,   !   3604  1.00000
     .-23439880d0,  38536313d0, 282391095d0,  36307534d0,-227863302d0,   !   3605  1.00000
     .-23388453d0,  38135324d0, 281346624d0,  37038564d0,-226380268d0,   !   3606  1.00000
     .-23336512d0,  37733536d0, 280296689d0,  37770925d0,-224890699d0,   !   3607  1.00000
     .-23284057d0,  37330947d0, 279241283d0,  38504621d0,-223394582d0,   !   3608  1.00000
     .-23231087d0,  36927557d0, 278180396d0,  39239650d0,-221891909d0,   !   3609  1.00000
     .-23177601d0,  36523365d0, 277114019d0,  39976015d0,-220382667d0,   !   3610  1.00000
     .-23123598d0,  36118372d0, 276042142d0,  40713715d0,-218866847d0,   !   3611  1.00000
     .-23069078d0,  35712576d0, 274964758d0,  41452752d0,-217344437d0,   !   3612  1.00000
     .-23014039d0,  35305977d0, 273881857d0,  42193127d0,-215815428d0,   !   3613  1.00000
     .-22958480d0,  34898575d0, 272793430d0,  42934840d0,-214279809d0,   !   3614  1.00000
     .-22902400d0,  34490369d0, 271699467d0,  43677892d0,-212737569d0,   !   3615  1.00000
     .-22845799d0,  34081359d0, 270599960d0,  44422285d0,-211188697d0,   !   3616  1.00000
     .-22788676d0,  33671543d0, 269494899d0,  45168018d0,-209633183d0,   !   3617  1.00000
     .-22731030d0,  33260922d0, 268384276d0,  45915093d0,-208071017d0,   !   3618  1.00000
     .-22672859d0,  32849495d0, 267268081d0,  46663511d0,-206502187d0,   !   3619  1.00000
     .-22614164d0,  32437262d0, 266146305d0,  47413272d0,-204926684d0,   !   3620  1.00000
     .-22554942d0,  32024221d0, 265018940d0,  48164377d0,-203344496d0,   !   3621  1.00000
     .-22495193d0,  31610374d0, 263885975d0,  48916828d0,-201755613d0,   !   3622  1.00000
     .-22434917d0,  31195718d0, 262747403d0,  49670624d0,-200160024d0,   !   3623  1.00000
     .-22374112d0,  30780254d0, 261603213d0,  50425767d0,-198557718d0,   !   3624  1.00000
     .-22312777d0,  30363981d0, 260453396d0,  51182258d0,-196948685d0,   !   3625  1.00000
     .-22250912d0,  29946898d0, 259297944d0,  51940097d0,-195332915d0,   !   3626  1.00000
     .-22188515d0,  29529005d0, 258136848d0,  52699285d0,-193710396d0,   !   3627  1.00000
     .-22125586d0,  29110302d0, 256970097d0,  53459824d0,-192081118d0,   !   3628  1.00000
     .-22062124d0,  28690788d0, 255797683d0,  54221713d0,-190445070d0,   !   3629  1.00000
     .-21998127d0,  28270463d0, 254619598d0,  54984954d0,-188802242d0,   !   3630  1.00000
     .-21933595d0,  27849326d0, 253435830d0,  55749548d0,-187152623d0,   !   3631  1.00000
     .-21868527d0,  27427376d0, 252246373d0,  56515495d0,-185496201d0,   !   3632  1.00000
     .-21802922d0,  27004613d0, 251051215d0,  57282796d0,-183832967d0,   !   3633  1.00000
     .-21736780d0,  26581036d0, 249850349d0,  58051453d0,-182162910d0,   !   3634  1.00000
     .-21670098d0,  26156646d0, 248643764d0,  58821466d0,-180486019d0,   !   3635  1.00000
     .-21602877d0,  25731441d0, 247431452d0,  59592835d0,-178802282d0,   !   3636  1.00000
     .-21535115d0,  25305421d0, 246213404d0,  60365562d0,-177111691d0,   !   3637  1.00000
     .-21466811d0,  24878586d0, 244989610d0,  61139647d0,-175414233d0,   !   3638  1.00000
     .-21397965d0,  24450935d0, 243760061d0,  61915092d0,-173709898d0,   !   3639  1.00000
     .-21328576d0,  24022467d0, 242524747d0,  62691897d0,-171998675d0,   !   3640  1.00000
     .-21258642d0,  23593182d0, 241283661d0,  63470062d0,-170280554d0,   !   3641  1.00000
     .-21188163d0,  23163080d0, 240036791d0,  64249590d0,-168555523d0,   !   3642  1.00000
     .-21117138d0,  22732160d0, 238784130d0,  65030480d0,-166823573d0,   !   3643  1.00000
     .-21045565d0,  22300421d0, 237525668d0,  65812734d0,-165084691d0,   !   3644  1.00000
     .-20973445d0,  21867864d0, 236261395d0,  66596351d0,-163338868d0,   !   3645  1.00000
     .-20900775d0,  21434486d0, 234991302d0,  67381335d0,-161586092d0,   !   3646  1.00000
     .-20827556d0,  21000289d0, 233715381d0,  68167684d0,-159826353d0,   !   3647  1.00000
     .-20753786d0,  20565272d0, 232433621d0,  68955399d0,-158059641d0,   !   3648  1.00000
     .-20679464d0,  20129433d0, 231146014d0,  69744483d0,-156285943d0,   !   3649  1.00000
     .-20604590d0,  19692773d0, 229852551d0,  70534935d0,-154505249d0,   !   3650  1.00000
     .-20529161d0,  19255291d0, 228553221d0,  71326756d0,-152717549d0,   !   3651  1.00000
     .-20453179d0,  18816986d0, 227248016d0,  72119947d0,-150922832d0,   !   3652  1.00000
     .-20376640d0,  18377858d0, 225936926d0,  72914510d0,-149121087d0,   !   3653  1.00000
     .-20299546d0,  17937907d0, 224619942d0,  73710444d0,-147312302d0,   !   3654  1.00000
     .-20221894d0,  17497132d0, 223297055d0,  74507751d0,-145496468d0,   !   3655  1.00000
     .-20143684d0,  17055532d0, 221968256d0,  75306431d0,-143673572d0,   !   3656  1.00000
     .-20064914d0,  16613108d0, 220633535d0,  76106486d0,-141843605d0,   !   3657  1.00000
     .-19985585d0,  16169857d0, 219292882d0,  76907916d0,-140006556d0,   !   3658  1.00000
     .-19905694d0,  15725781d0, 217946289d0,  77710722d0,-138162413d0,   !   3659  1.00000
     .-19825242d0,  15280878d0, 216593746d0,  78514904d0,-136311166d0,   !   3660  1.00000
     .-19744226d0,  14835149d0, 215235244d0,  79320465d0,-134452804d0,   !   3661  1.00000
     .-19662647d0,  14388591d0, 213870774d0,  80127404d0,-132587316d0,   !   3662  1.00000
     .-19580502d0,  13941206d0, 212500325d0,  80935722d0,-130714690d0,   !   3663  1.00000
     .-19497792d0,  13492992d0, 211123889d0,  81745420d0,-128834917d0,   !   3664  1.00000
     .-19414516d0,  13043949d0, 209741457d0,  82556499d0,-126947985d0,   !   3665  1.00000
     .-19330672d0,  12594077d0, 208353018d0,  83368961d0,-125053884d0,   !   3666  1.00000
     .-19246259d0,  12143375d0, 206958564d0,  84182805d0,-123152601d0,   !   3667  1.00000
     .-19161276d0,  11691842d0, 205558086d0,  84998032d0,-121244127d0,   !   3668  1.00000
     .-19075724d0,  11239478d0, 204151573d0,  85814644d0,-119328451d0,   !   3669  1.00000
     .-18989600d0,  10786282d0, 202739017d0,  86632641d0,-117405560d0,   !   3670  1.00000
     .-18902903d0,  10332255d0, 201320407d0,  87452024d0,-115475446d0,   !   3671  1.00000
     .-18815633d0,   9877395d0, 199895736d0,  88272793d0,-113538096d0,   !   3672  1.00000
     .-18727790d0,   9421702d0, 198464992d0,  89094951d0,-111593499d0,   !   3673  1.00000
     .-18639370d0,   8965175d0, 197028168d0,  89918497d0,-109641645d0,   !   3674  1.00000
     .-18550375d0,   8507815d0, 195585253d0,  90743432d0,-107682523d0,   !   3675  1.00000
     .-18460803d0,   8049620d0, 194136237d0,  91569758d0,-105716121d0,   !   3676  1.00000
     .-18370653d0,   7590589d0, 192681113d0,  92397475d0,-103742429d0,   !   3677  1.00000
     .-18279924d0,   7130724d0, 191219869d0,  93226583d0,-101761436d0,   !   3678  1.00000
     .-18188615d0,   6670022d0, 189752497d0,  94057084d0, -99773130d0,   !   3679  1.00000
     .-18096726d0,   6208483d0, 188278987d0,  94888979d0, -97777501d0,   !   3680  1.00000
     .-18004254d0,   5746108d0, 186799330d0,  95722268d0, -95774537d0,   !   3681  1.00000
     .-17911200d0,   5282895d0, 185313517d0,  96556952d0, -93764228d0,   !   3682  1.00000
     .-17817562d0,   4818844d0, 183821537d0,  97393033d0, -91746562d0,   !   3683  1.00000
     .-17723340d0,   4353955d0, 182323381d0,  98230510d0, -89721529d0,   !   3684  1.00000
     .-17628532d0,   3888226d0, 180819041d0,  99069385d0, -87689117d0,   !   3685  1.00000
     .-17533138d0,   3421658d0, 179308505d0,  99909659d0, -85649315d0,   !   3686  1.00000
     .-17437156d0,   2954250d0, 177791766d0, 100751332d0, -83602113d0,   !   3687  1.00000
     .-17340586d0,   2486001d0, 176268813d0, 101594405d0, -81547498d0,   !   3688  1.00000
     .-17243426d0,   2016911d0, 174739637d0, 102438879d0, -79485461d0,   !   3689  1.00000
     .-17145676d0,   1546980d0, 173204228d0, 103284755d0, -77415990d0,   !   3690  1.00000
     .-17047335d0,   1076206d0, 171662577d0, 104132034d0, -75339074d0,   !   3691  1.00000
     .-16948402d0,    604590d0, 170114674d0, 104980716d0, -73254701d0,   !   3692  1.00000
     .-16848875d0,    132131d0, 168560510d0, 105830803d0, -71162861d0,   !   3693  1.00000
     .-16748755d0,   -341170d0, 167000076d0, 106682295d0, -69063543d0,   !   3694  1.00000
     .-16648039d0,   -815317d0, 165433361d0, 107535193d0, -66956736d0,   !   3695  1.00000
     .-16546728d0,  -1290308d0, 163860357d0, 108389498d0, -64842427d0,   !   3696  1.00000
     .-16444819d0,  -1766144d0, 162281053d0, 109245211d0, -62720607d0,   !   3697  1.00000
     .-16342313d0,  -2242826d0, 160695440d0, 110102332d0, -60591264d0,   !   3698  1.00000
     .-16239207d0,  -2720354d0, 159103509d0, 110960863d0, -58454387d0,   !   3699  1.00000
     .-16135502d0,  -3198729d0, 157505250d0, 111820804d0, -56309964d0,   !   3700  1.00000
     .-16031196d0,  -3677951d0, 155900653d0, 112682156d0, -54157986d0,   !   3701  1.00000
     .-15926288d0,  -4158020d0, 154289709d0, 113544920d0, -51998439d0,   !   3702  1.00000
     .-15820778d0,  -4638937d0, 152672408d0, 114409097d0, -49831314d0,   !   3703  1.00000
     .-15714664d0,  -5120703d0, 151048741d0, 115274687d0, -47656599d0,   !   3704  1.00000
     .-15607945d0,  -5603318d0, 149418698d0, 116141692d0, -45474283d0,   !   3705  1.00000
     .-15500621d0,  -6086782d0, 147782270d0, 117010112d0, -43284355d0,   !   3706  1.00000
     .-15392690d0,  -6571097d0, 146139446d0, 117879948d0, -41086803d0,   !   3707  1.00000
     .-15284152d0,  -7056262d0, 144490218d0, 118751202d0, -38881617d0,   !   3708  1.00000
     .-15175005d0,  -7542278d0, 142834576d0, 119623872d0, -36668784d0,   !   3709  1.00000
     .-15065248d0,  -8029146d0, 141172509d0, 120497962d0, -34448295d0,   !   3710  1.00000
     .-14954882d0,  -8516866d0, 139504009d0, 121373471d0, -32220137d0,   !   3711  1.00000
     .-14843904d0,  -9005439d0, 137829065d0, 122250400d0, -29984300d0,   !   3712  1.00000
     .-14732313d0,  -9494865d0, 136147669d0, 123128750d0, -27740771d0,   !   3713  1.00000
     .-14620110d0,  -9985144d0, 134459810d0, 124008522d0, -25489541d0,   !   3714  1.00000
     .-14507292d0, -10476277d0, 132765480d0, 124889718d0, -23230598d0,   !   3715  1.00000
     .-14393859d0, -10968265d0, 131064667d0, 125772336d0, -20963929d0,   !   3716  1.00000
     .-14279809d0, -11461108d0, 129357363d0, 126656379d0, -18689525d0,   !   3717  1.00000
     .-14165143d0, -11954807d0, 127643557d0, 127541848d0, -16407374d0,   !   3718  1.00000
     .-14049858d0, -12449362d0, 125923241d0, 128428742d0, -14117465d0,   !   3719  1.00000
     .-13933954d0, -12944773d0, 124196405d0, 129317064d0, -11819786d0,   !   3720  1.00000
     .-13817431d0, -13441041d0, 122463038d0, 130206813d0,  -9514326d0,   !   3721  1.00000
     .-13700286d0, -13938168d0, 120723131d0, 131097991d0,  -7201073d0,   !   3722  1.00000
     .-13582519d0, -14436152d0, 118976675d0, 131990598d0,  -4880017d0,   !   3723  1.00000
     .-13464129d0, -14934995d0, 117223660d0, 132884636d0,  -2551146d0,   !   3724  1.00000
     .-13345116d0, -15434696d0, 115464075d0, 133780104d0,   -214449d0,   !   3725  1.00000
     .-13225477d0, -15935258d0, 113697912d0, 134677005d0,   2130084d0,   !   3726  1.00000
     .-13105212d0, -16436680d0, 111925161d0, 135575339d0,   4482467d0,   !   3727  1.00000
     .-12984321d0, -16938962d0, 110145811d0, 136475106d0,   6842711d0,   !   3728  1.00000
     .-12862802d0, -17442105d0, 108359854d0, 137376308d0,   9210826d0,   !   3729  1.00000
     .-12740654d0, -17946110d0, 106567279d0, 138278945d0,  11586824d0,   !   3730  1.00000
     .-12617876d0, -18450977d0, 104768077d0, 139183018d0,  13970716d0,   !   3731  1.00000
     .-12494468d0, -18956707d0, 102962237d0, 140088528d0,  16362514d0,   !   3732  1.00000
     .-12370427d0, -19463300d0, 101149751d0, 140995476d0,  18762229d0,   !   3733  1.00000
     .-12245754d0, -19970757d0,  99330609d0, 141903863d0,  21169872d0,   !   3734  1.00000
     .-12120448d0, -20479077d0,  97504800d0, 142813689d0,  23585455d0,   !   3735  1.00000
     .-11994507d0, -20988262d0,  95672314d0, 143724955d0,  26008989d0,   !   3736  1.00000
     .-11867930d0, -21498313d0,  93833143d0, 144637663d0,  28440485d0,   !   3737  1.00000
     .-11740716d0, -22009229d0,  91987277d0, 145551812d0,  30879955d0,   !   3738  1.00000
     .-11612865d0, -22521011d0,  90134705d0, 146467405d0,  33327411d0,   !   3739  1.00000
     .-11484375d0, -23033660d0,  88275417d0, 147384441d0,  35782864d0,   !   3740  1.00000
     .-11355246d0, -23547176d0,  86409405d0, 148302921d0,  38246324d0,   !   3741  1.00000
     .-11225476d0, -24061559d0,  84536658d0, 149222847d0,  40717805d0,   !   3742  1.00000
     .-11095065d0, -24576811d0,  82657166d0, 150144219d0,  43197316d0,   !   3743  1.00000
     .-10964011d0, -25092931d0,  80770920d0, 151067038d0,  45684870d0,   !   3744  1.00000
     .-10832314d0, -25609921d0,  78877910d0, 151991305d0,  48180477d0,   !   3745  1.00000
     .-10699972d0, -26127780d0,  76978125d0, 152917021d0,  50684150d0,   !   3746  1.00000
     .-10566985d0, -26646509d0,  75071556d0, 153844186d0,  53195900d0,   !   3747  1.00000
     .-10433351d0, -27166109d0,  73158194d0, 154772802d0,  55715739d0,   !   3748  1.00000
     .-10299070d0, -27686580d0,  71238028d0, 155702868d0,  58243677d0,   !   3749  1.00000
     .-10164140d0, -28207923d0,  69311049d0, 156634387d0,  60779726d0,   !   3750  1.00000
     .-10028561d0, -28730138d0,  67377247d0, 157567359d0,  63323898d0,   !   3751  1.00000
     . -9892332d0, -29253226d0,  65436611d0, 158501784d0,  65876205d0,   !   3752  1.00000
     . -9755452d0, -29777186d0,  63489132d0, 159437664d0,  68436657d0,   !   3753  1.00000
     . -9617919d0, -30302021d0,  61534800d0, 160374999d0,  71005267d0,   !   3754  1.00000
     . -9479732d0, -30827729d0,  59573606d0, 161313790d0,  73582045d0,   !   3755  1.00000
     . -9340892d0, -31354312d0,  57605539d0, 162254039d0,  76167003d0,   !   3756  1.00000
     . -9201396d0, -31881771d0,  55630589d0, 163195745d0,  78760154d0,   !   3757  1.00000
     . -9061244d0, -32410105d0,  53648747d0, 164138910d0,  81361508d0,   !   3758  1.00000
     . -8920434d0, -32939315d0,  51660003d0, 165083534d0,  83971077d0,   !   3759  1.00000
     . -8778966d0, -33469402d0,  49664346d0, 166029619d0,  86588872d0,   !   3760  1.00000
     . -8636839d0, -34000366d0,  47661768d0, 166977165d0,  89214905d0,   !   3761  1.00000
     . -8494052d0, -34532208d0,  45652257d0, 167926174d0,  91849188d0,   !   3762  1.00000
     . -8350603d0, -35064927d0,  43635804d0, 168876645d0,  94491732d0,   !   3763  1.00000
     . -8206492d0, -35598526d0,  41612399d0, 169828579d0,  97142549d0,   !   3764  1.00000
     . -8061718d0, -36133004d0,  39582033d0, 170781979d0,  99801650d0,   !   3765  1.00000
     . -7916280d0, -36668361d0,  37544694d0, 171736843d0, 102469047d0,   !   3766  1.00000
     . -7770176d0, -37204599d0,  35500374d0, 172693174d0, 105144752d0,   !   3767  1.00000
     . -7623406d0, -37741717d0,  33449063d0, 173650972d0, 107828775d0,   !   3768  1.00000
     . -7475969d0, -38279716d0,  31390749d0, 174610238d0, 110521130d0,   !   3769  1.00000
     . -7327864d0, -38818597d0,  29325424d0, 175570973d0, 113221826d0,   !   3770  1.00000
     . -7179089d0, -39358361d0,  27253078d0, 176533177d0, 115930877d0,   !   3771  1.00000
     . -7029645d0, -39899006d0,  25173700d0, 177496852d0, 118648293d0,   !   3772  1.00000
     . -6879529d0, -40440536d0,  23087280d0, 178461997d0, 121374087d0,   !   3773  1.00000
     . -6728740d0, -40982948d0,  20993809d0, 179428615d0, 124108269d0,   !   3774  1.00000
     . -6577279d0, -41526246d0,  18893277d0, 180396706d0, 126850852d0,   !   3775  1.00000
     . -6425143d0, -42070427d0,  16785673d0, 181366271d0, 129601847d0,   !   3776  1.00000
     . -6272333d0, -42615494d0,  14670988d0, 182337310d0, 132361266d0,   !   3777  1.00000
     . -6118845d0, -43161447d0,  12549211d0, 183309824d0, 135129120d0,   !   3778  1.00000
     . -5964681d0, -43708286d0,  10420332d0, 184283815d0, 137905422d0,   !   3779  1.00000
     . -5809839d0, -44256011d0,   8284342d0, 185259283d0, 140690182d0,   !   3780  1.00000
     . -5654317d0, -44804624d0,   6141231d0, 186236229d0, 143483413d0,   !   3781  1.00000
     . -5498115d0, -45354124d0,   3990988d0, 187214653d0, 146285126d0,   !   3782  1.00000
     . -5341232d0, -45904513d0,   1833603d0, 188194558d0, 149095333d0,   !   3783  1.00000
     . -5183666d0, -46455790d0,   -330932d0, 189175942d0, 151914046d0,   !   3784  1.00000
     . -5025418d0, -47007957d0,  -2502630d0, 190158808d0, 154741276d0,   !   3785  1.00000
     . -4866485d0, -47561013d0,  -4681499d0, 191143156d0, 157577035d0,   !   3786  1.00000
     . -4706867d0, -48114959d0,  -6867551d0, 192128987d0, 160421335d0,   !   3787  1.00000
     . -4546562d0, -48669796d0,  -9060794d0, 193116302d0, 163274187d0,   !   3788  1.00000
     . -4385570d0, -49225525d0, -11261239d0, 194105101d0, 166135604d0,   !   3789  1.00000
     . -4223890d0, -49782145d0, -13468896d0, 195095386d0, 169005597d0,   !   3790  1.00000
     . -4061521d0, -50339657d0, -15683776d0, 196087156d0, 171884177d0,   !   3791  1.00000
     . -3898462d0, -50898062d0, -17905887d0, 197080415d0, 174771357d0,   !   3792  1.00000
     . -3734711d0, -51457360d0, -20135241d0, 198075161d0, 177667149d0,   !   3793  1.00000
     . -3570268d0, -52017552d0, -22371847d0, 199071395d0, 180571563d0,   !   3794  1.00000
     . -3405131d0, -52578638d0, -24615715d0, 200069120d0, 183484613d0,   !   3795  1.00000
     . -3239300d0, -53140619d0, -26866856d0, 201068335d0, 186406309d0,   !   3796  1.00000
     . -3072774d0, -53703495d0, -29125280d0, 202069041d0, 189336663d0,   !   3797  1.00000
     . -2905551d0, -54267267d0, -31390997d0, 203071239d0, 192275688d0,   !   3798  1.00000
     . -2737632d0, -54831935d0, -33664016d0, 204074931d0, 195223395d0,   !   3799  1.00000
     . -2569013d0, -55397500d0, -35944348d0, 205080116d0, 198179795d0,   !   3800  1.00000
     . -2399695d0, -55963962d0, -38232004d0, 206086795d0, 201144902d0,   !   3801  1.00000
     . -2229677d0, -56531321d0, -40526993d0, 207094970d0, 204118725d0,   !   3802  1.00000
     . -2058958d0, -57099579d0, -42829326d0, 208104642d0, 207101278d0,   !   3803  1.00000
     . -1887536d0, -57668736d0, -45139012d0, 209115811d0, 210092573d0,   !   3804  1.00000
     . -1715410d0, -58238792d0, -47456062d0, 210128477d0, 213092620d0,   !   3805  1.00000
     . -1542580d0, -58809748d0, -49780486d0, 211142642d0, 216101432d0,   !   3806  1.00000
     . -1369045d0, -59381604d0, -52112293d0, 212158307d0, 219119021d0,   !   3807  1.00000
     . -1194803d0, -59954360d0, -54451496d0, 213175473d0, 222145399d0,   !   3808  1.00000
     . -1019853d0, -60528018d0, -56798103d0, 214194140d0, 225180577d0,   !   3809  1.00000
     .  -844195d0, -61102578d0, -59152124d0, 215214309d0, 228224568d0,   !   3810  1.00000
     .  -667828d0, -61678040d0, -61513570d0, 216235981d0, 231277382d0,   !   3811  1.00000
     .  -490749d0, -62254404d0, -63882452d0, 217259156d0, 234339033d0,   !   3812  1.00000
     .  -312960d0, -62831672d0, -66258778d0, 218283837d0, 237409532d0,   !   3813  1.00000
     .  -134457d0, -63409844d0, -68642560d0, 219310023d0, 240488891d0,   !   3814  1.00000
     .    44758d0, -63988919d0, -71033808d0, 220337715d0, 243577122d0,   !   3815  1.00000
     .   224689d0, -64568900d0, -73432532d0, 221366915d0, 246674236d0,   !   3816  1.00000
     .   405336d0, -65149786d0, -75838742d0, 222397622d0, 249780246d0,   !   3817  1.00000
     .   586699d0, -65731578d0, -78252448d0, 223429838d0, 252895164d0,   !   3818  1.00000
     .   768780d0, -66314275d0, -80673661d0, 224463564d0, 256019002d0,   !   3819  1.00000
     .   951580d0, -66897880d0, -83102391d0, 225498800d0, 259151770d0,   !   3820  1.00000
     .  1135099d0, -67482392d0, -85538648d0, 226535548d0, 262293483d0,   !   3821  1.00000
     .  1319340d0, -68067812d0, -87982442d0, 227573808d0, 265444151d0,   !   3822  1.00000
     .  1504302d0, -68654140d0, -90433784d0, 228613581d0, 268603786d0,   !   3823  1.00000
     .  1689987d0, -69241376d0, -92892683d0, 229654868d0, 271772401d0,   !   3824  1.00000
     .  1876396d0, -69829522d0, -95359151d0, 230697669d0, 274950007d0,   !   3825  1.00000
     .  2063530d0, -70418578d0, -97833198d0, 231741986d0, 278136616d0,   !   3826  1.00000
     .  2251390d0, -71008544d0,-100314833d0, 232787820d0, 281332240d0,   !   3827  1.00000
     .  2439977d0, -71599421d0,-102804067d0, 233835170d0, 284536892d0,   !   3828  1.00000
     .  2629292d0, -72191209d0,-105300910d0, 234884039d0, 287750584d0,   !   3829  1.00000
     .  2819336d0, -72783908d0,-107805373d0, 235934426d0, 290973326d0,   !   3830  1.00000
     .  3010111d0, -73377521d0,-110317466d0, 236986333d0, 294205132d0,   !   3831  1.00000
     .  3201616d0, -73972045d0,-112837199d0, 238039761d0, 297446014d0,   !   3832  1.00000
     .  3393853d0, -74567484d0,-115364583d0, 239094710d0, 300695983d0,   !   3833  1.00000
     .  3586824d0, -75163835d0,-117899628d0, 240151181d0, 303955051d0,   !   3834  1.00000
     .  3780529d0, -75761102d0,-120442344d0, 241209175d0, 307223231d0,   !   3835  1.00000
     .  3974969d0, -76359283d0,-122992741d0, 242268694d0, 310500534d0,   !   3836  1.00000
     .  4170145d0, -76958379d0,-125550830d0, 243329736d0, 313786973d0,   !   3837  1.00000
     .  4366059d0, -77558391d0,-128116622d0, 244392305d0, 317082559d0,   !   3838  1.00000
     .  4562711d0, -78159319d0,-130690126d0, 245456400d0, 320387306d0,   !   3839  1.00000
     .  4760103d0, -78761164d0,-133271353d0, 246522021d0, 323701224d0,   !   3840  1.00000
     .  4958235d0, -79363927d0,-135860313d0, 247589171d0, 327024326d0,   !   3841  1.00000
     .  5157108d0, -79967607d0,-138457017d0, 248657850d0, 330356624d0,   !   3842  1.00000
     .  5356724d0, -80572205d0,-141061475d0, 249728059d0, 333698130d0,   !   3843  1.00000
     .  5557084d0, -81177723d0,-143673697d0, 250799798d0, 337048857d0,   !   3844  1.00000
     .  5758188d0, -81784159d0,-146293695d0, 251873069d0, 340408815d0,   !   3845  1.00000
     .  5960037d0, -82391516d0,-148921477d0, 252947871d0, 343778018d0,   !   3846  1.00000
     .  6162634d0, -82999792d0,-151557055d0, 254024207d0, 347156478d0,   !   3847  1.00000
     .  6365978d0, -83608990d0,-154200439d0, 255102077d0, 350544207d0,   !   3848  1.00000
     .  6570070d0, -84219109d0,-156851639d0, 256181482d0, 353941216d0,   !   3849  1.00000
     .  6774913d0, -84830149d0,-159510667d0, 257262422d0, 357347518d0,   !   3850  1.00000
     .  6980506d0, -85442112d0,-162177531d0, 258344898d0, 360763125d0,   !   3851  1.00000
     .  7186852d0, -86054998d0,-164852243d0, 259428912d0, 364188050d0,   !   3852  1.00000
     .  7393950d0, -86668807d0,-167534814d0, 260514464d0, 367622304d0,   !   3853  1.00000
     .  7601802d0, -87283540d0,-170225253d0, 261601554d0, 371065899d0,   !   3854  1.00000
     .  7810409d0, -87899197d0,-172923570d0, 262690185d0, 374518849d0,   !   3855  1.00000
     .  8019773d0, -88515779d0,-175629778d0, 263780356d0, 377981164d0,   !   3856  1.00000
     .  8229893d0, -89133286d0,-178343885d0, 264872068d0, 381452858d0,   !   3857  1.00000
     .  8440771d0, -89751719d0,-181065902d0, 265965323d0, 384933942d0,   !   3858  1.00000
     .  8652409d0, -90371078d0,-183795840d0, 267060120d0, 388424428d0,   !   3859  1.00000
     .  8864807d0, -90991365d0,-186533710d0, 268156462d0, 391924329d0,   !   3860  1.00000
     .  9077966d0, -91612578d0,-189279521d0, 269254349d0, 395433657d0,   !   3861  1.00000
     .  9291888d0, -92234720d0,-192033285d0, 270353781d0, 398952425d0,   !   3862  1.00000
     .  9506573d0, -92857789d0,-194795011d0, 271454759d0, 402480644d0,   !   3863  1.00000
     .  9722022d0, -93481788d0,-197564711d0, 272557285d0, 406018326d0,   !   3864  1.00000
     .  9938237d0, -94106716d0,-200342394d0, 273661359d0, 409565485d0,   !   3865  1.00000
     . 10155218d0, -94732573d0,-203128071d0, 274766981d0, 413122131d0,   !   3866  1.00000
     . 10372967d0, -95359361d0,-205921753d0, 275874154d0, 416688278d0,   !   3867  1.00000
     . 10591485d0, -95987080d0,-208723451d0, 276982877d0, 420263938d0,   !   3868  1.00000
     . 10810772d0, -96615730d0,-211533174d0, 278093152d0, 423849123d0,   !   3869  1.00000
     . 11030831d0, -97245312d0,-214350934d0, 279204978d0, 427443845d0,   !   3870  1.00000
     . 11251660d0, -97875826d0,-217176740d0, 280318358d0, 431048116d0,   !   3871  1.00000
     . 11473263d0, -98507274d0,-220010604d0, 281433292d0, 434661949d0,   !   3872  1.00000
     . 11695640d0, -99139654d0,-222852535d0, 282549781d0, 438285357d0,   !   3873  1.00000
     . 11918792d0, -99772968d0,-225702545d0, 283667825d0, 441918350d0,   !   3874  1.00000
     . 12142719d0,-100407217d0,-228560645d0, 284787426d0, 445560943d0,   !   3875  1.00000
     . 12367424d0,-101042401d0,-231426843d0, 285908584d0, 449213146d0,   !   3876  1.00000
     . 12592907d0,-101678520d0,-234301152d0, 287031300d0, 452874973d0,   !   3877  1.00000
     . 12819169d0,-102315574d0,-237183582d0, 288155576d0, 456546436d0,   !   3878  1.00000
     . 13046211d0,-102953565d0,-240074143d0, 289281411d0, 460227546d0,   !   3879  1.00000
     . 13274035d0,-103592493d0,-242972846d0, 290408806d0, 463918317d0,   !   3880  1.00000
     . 13502641d0,-104232359d0,-245879701d0, 291537763d0, 467618761d0,   !   3881  1.00000
     . 13732030d0,-104873162d0,-248794719d0, 292668283d0, 471328890d0,   !   3882  1.00000
     . 13962204d0,-105514903d0,-251717911d0, 293800365d0, 475048716d0,   !   3883  1.00000
     . 14193164d0,-106157583d0,-254649288d0, 294934012d0, 478778252d0,   !   3884  1.00000
     . 14424910d0,-106801203d0,-257588859d0, 296069223d0, 482517510d0,   !   3885  1.00000
     . 14657444d0,-107445763d0,-260536636d0, 297206000d0, 486266502d0,   !   3886  1.00000
     . 14890766d0,-108091263d0,-263492628d0, 298344343d0, 490025242d0,   !   3887  1.00000
     . 15124879d0,-108737703d0,-266456848d0, 299484254d0, 493793741d0,   !   3888  1.00000
     . 15359782d0,-109385085d0,-269429305d0, 300625733d0, 497572011d0,   !   3889  1.00000
     . 15595477d0,-110033409d0,-272410010d0, 301768781d0, 501360066d0,   !   3890  1.00000
     . 15831965d0,-110682676d0,-275398974d0, 302913398d0, 505157918d0,   !   3891  1.00000
     . 16069247d0,-111332885d0,-278396207d0, 304059587d0, 508965578d0,   !   3892  1.00000
     . 16307325d0,-111984038d0,-281401720d0, 305207346d0, 512783060d0,   !   3893  1.00000
     . 16546198d0,-112636134d0,-284415524d0, 306356679d0, 516610375d0,   !   3894  1.00000
     . 16785869d0,-113289175d0,-287437629d0, 307507584d0, 520447537d0,   !   3895  1.00000
     . 17026338d0,-113943161d0,-290468047d0, 308660063d0, 524294558d0,   !   3896  1.00000
     . 17267607d0,-114598092d0,-293506786d0, 309814117d0, 528151450d0,   !   3897  1.00000
     . 17509676d0,-115253969d0,-296553860d0, 310969747d0, 532018225d0,   !   3898  1.00000
     . 17752546d0,-115910792d0,-299609277d0, 312126954d0, 535894897d0,   !   3899  1.00000
     . 17996219d0,-116568562d0,-302673049d0, 313285737d0, 539781477d0,   !   3900  1.00000
     . 18240696d0,-117227280d0,-305745187d0, 314446099d0, 543677978d0,   !   3901  1.00000
     . 18485978d0,-117886945d0,-308825700d0, 315608040d0, 547584413d0,   !   3902  1.00000
     . 18732065d0,-118547559d0,-311914601d0, 316771561d0, 551500794d0,   !   3903  1.00000
     . 18978960d0,-119209122d0,-315011899d0, 317936662d0, 555427133d0,   !   3904  1.00000
     . 19226662d0,-119871634d0,-318117606d0, 319103345d0, 559363444d0,   !   3905  1.00000
     . 19475173d0,-120535096d0,-321231731d0, 320271610d0, 563309738d0,   !   3906  1.00000
     . 19724495d0,-121199509d0,-324354287d0, 321441459d0, 567266028d0,   !   3907  1.00000
     . 19974627d0,-121864872d0,-327485283d0, 322612891d0, 571232327d0,   !   3908  1.00000
     . 20225572d0,-122531187d0,-330624730d0, 323785908d0, 575208646d0,   !   3909  1.00000
     . 20477331d0,-123198453d0,-333772640d0, 324960511d0, 579195000d0,   !   3910  1.00000
     . 20729903d0,-123866672d0,-336929022d0, 326136701d0, 583191400d0,   !   3911  1.00000
     . 20983291d0,-124535844d0,-340093888d0, 327314478d0, 587197859d0,   !   3912  1.00000
     . 21237496d0,-125205970d0,-343267248d0, 328493843d0, 591214389d0,   !   3913  1.00000
     . 21492519d0,-125877049d0,-346449113d0, 329674797d0, 595241003d0,   !   3914  1.00000
     . 21748360d0,-126549082d0,-349639494d0, 330857341d0, 599277713d0,   !   3915  1.00000
     . 22005021d0,-127222071d0,-352838402d0, 332041476d0, 603324533d0,   !   3916  1.00000
     . 22262503d0,-127896015d0,-356045848d0, 333227202d0, 607381474d0,   !   3917  1.00000
     . 22520807d0,-128570915d0,-359261841d0, 334414521d0, 611448549d0,   !   3918  1.00000
     . 22779934d0,-129246771d0,-362486394d0, 335603433d0, 615525772d0,   !   3919  1.00000
     . 23039885d0,-129923584d0,-365719517d0, 336793939d0, 619613154d0,   !   3920  1.00000
     . 23300661d0,-130601355d0,-368961220d0, 337986040d0, 623710708d0,   !   3921  1.00000
     . 23562263d0,-131280083d0,-372211515d0, 339179736d0, 627818446d0,   !   3922  1.00000
     . 23824693d0,-131959770d0,-375470413d0, 340375029d0, 631936382d0,   !   3923  1.00000
     . 24087952d0,-132640416d0,-378737923d0, 341571920d0, 636064528d0,   !   3924  1.00000
     . 24352040d0,-133322021d0,-382014058d0, 342770409d0, 640202896d0,   !   3925  1.00000
     . 24616958d0,-134004586d0,-385298827d0, 343970496d0, 644351500d0,   !   3926  1.00000
     . 24882709d0,-134688112d0,-388592243d0, 345172184d0, 648510352d0,   !   3927  1.00000
     . 25149292d0,-135372598d0,-391894315d0, 346375472d0, 652679464d0,   !   3928  1.00000
     . 25416709d0,-136058046d0,-395205054d0, 347580362d0, 656858849d0,   !   3929  1.00000
     . 25684961d0,-136744456d0,-398524472d0, 348786855d0, 661048521d0,   !   3930  1.00000
     . 25954049d0,-137431828d0,-401852579d0, 349994950d0, 665248491d0,   !   3931  1.00000
     . 26223974d0,-138120164d0,-405189386d0, 351204650d0, 669458772d0,   !   3932  1.00000
     . 26494738d0,-138809462d0,-408534904d0, 352415955d0, 673679377d0,   !   3933  1.00000
     . 26766341d0,-139499725d0,-411889143d0, 353628865d0, 677910319d0,   !   3934  1.00000
     . 27038784d0,-140190952d0,-415252116d0, 354843382d0, 682151610d0,   !   3935  1.00000
     . 27312069d0,-140883144d0,-418623832d0, 356059507d0, 686403263d0,   !   3936  1.00000
     . 27586196d0,-141576301d0,-422004303d0, 357277239d0, 690665291d0,   !   3937  1.00000
     . 27861167d0,-142270424d0,-425393539d0, 358496581d0, 694937707d0,   !   3938  1.00000
     . 28136983d0,-142965514d0,-428791552d0, 359717532d0, 699220522d0,   !   3939  1.00000
     . 28413645d0,-143661571d0,-432198352d0, 360940095d0, 703513751d0,   !   3940  1.00000
     . 28691154d0,-144358595d0,-435613951d0, 362164269d0, 707817406d0,   !   3941  1.00000
     . 28969511d0,-145056587d0,-439038358d0, 363390055d0, 712131499d0,   !   3942  1.00000
     . 29248717d0,-145755548d0,-442471586d0, 364617455d0, 716456043d0,   !   3943  1.00000
     . 29528773d0,-146455477d0,-445913644d0, 365846469d0, 720791051d0,   !   3944  1.00000
     . 29809681d0,-147156376d0,-449364545d0, 367077097d0, 725136536d0,   !   3945  1.00000
     . 30091441d0,-147858245d0,-452824299d0, 368309342d0, 729492511d0,   !   3946  1.00000
     . 30374055d0,-148561084d0,-456292916d0, 369543203d0, 733858988d0,   !   3947  1.00000
     . 30657523d0,-149264894d0,-459770409d0, 370778682d0, 738235980d0,   !   3948  1.00000
     . 30941847d0,-149969676d0,-463256788d0, 372015778d0, 742623500d0,   !   3949  1.00000
     . 31227028d0,-150675430d0,-466752063d0, 373254494d0, 747021560d0,   !   3950  1.00000
     . 31513067d0,-151382156d0,-470256246d0, 374494830d0, 751430175d0,   !   3951  1.00000
     . 31799965d0,-152089855d0,-473769348d0, 375736787d0, 755849355d0,   !   3952  1.00000
     . 32087723d0,-152798528d0,-477291380d0, 376980365d0, 760279115d0,   !   3953  1.00000
     . 32376342d0,-153508174d0,-480822353d0, 378225566d0, 764719467d0,   !   3954  1.00000
     . 32665824d0,-154218795d0,-484362278d0, 379472390d0, 769170424d0,   !   3955  1.00000
     . 32956169d0,-154930391d0,-487911165d0, 380720838d0, 773631998d0,   !   3956  1.00000
     . 33247379d0,-155642962d0,-491469027d0, 381970911d0, 778104203d0,   !   3957  1.00000
     . 33539454d0,-156356510d0,-495035873d0, 383222610d0, 782587052d0,   !   3958  1.00000
     . 33832396d0,-157071034d0,-498611716d0, 384475936d0, 787080556d0,   !   3959  1.00000
     . 34126206d0,-157786535d0,-502196565d0, 385730889d0, 791584730d0,   !   3960  1.00000
     . 34420885d0,-158503013d0,-505790433d0, 386987470d0, 796099586d0,   !   3961  1.00000
     . 34716434d0,-159220469d0,-509393329d0, 388245681d0, 800625137d0,   !   3962  1.00000
     . 35012854d0,-159938904d0,-513005266d0, 389505521d0, 805161396d0,   !   3963  1.00000
     . 35310147d0,-160658318d0,-516626254d0, 390766992d0, 809708375d0,   !   3964  1.00000
     . 35608313d0,-161378711d0,-520256304d0, 392030095d0, 814266088d0,   !   3965  1.00000
     . 35907353d0,-162100084d0,-523895428d0, 393294831d0, 818834547d0,   !   3966  1.00000
     . 36207269d0,-162822438d0,-527543636d0, 394561199d0, 823413766d0,   !   3967  1.00000
     . 36508062d0,-163545772d0,-531200939d0, 395829202d0, 828003757d0,   !   3968  1.00000
     . 36809733d0,-164270089d0,-534867349d0, 397098840d0, 832604534d0,   !   3969  1.00000
     . 37112282d0,-164995387d0,-538542877d0, 398370114d0, 837216108d0,   !   3970  1.00000
     . 37415712d0,-165721667d0,-542227534d0, 399643024d0, 841838494d0,   !   3971  1.00000
     . 37720023d0,-166448931d0,-545921330d0, 400917572d0, 846471704d0,   !   3972  1.00000
     . 38025216d0,-167177178d0,-549624278d0, 402193758d0, 851115751d0,   !   3973  1.00000
     . 38331292d0,-167906409d0,-553336387d0, 403471583d0, 855770649d0,   !   3974  1.00000
     . 38638253d0,-168636625d0,-557057670d0, 404751049d0, 860436409d0,   !   3975  1.00000
     . 38946100d0,-169367825d0,-560788137d0, 406032155d0, 865113045d0,   !   3976  1.00000
     . 39254833d0,-170100011d0,-564527800d0, 407314903d0, 869800570d0,   !   3977  1.00000
     . 39564454d0,-170833184d0,-568276669d0, 408599293d0, 874498997d0,   !   3978  1.00000
     . 39874964d0,-171567342d0,-572034757d0, 409885327d0, 879208339d0,   !   3979  1.00000
     . 40186364d0,-172302488d0,-575802073d0, 411173004d0, 883928609d0,   !   3980  1.00000
     . 40498656d0,-173038621d0,-579578629d0, 412462327d0, 888659820d0,   !   3981  1.00000
     . 40811839d0,-173775742d0,-583364437d0, 413753296d0, 893401985d0,   !   3982  1.00000
     . 41125916d0,-174513851d0,-587159507d0, 415045911d0, 898155117d0,   !   3983  1.00000
     . 41440888d0,-175252950d0,-590963850d0, 416340174d0, 902919229d0,   !   3984  1.00000
     . 41756755d0,-175993038d0,-594777479d0, 417636085d0, 907694334d0,   !   3985  1.00000
     . 42073519d0,-176734116d0,-598600403d0, 418933646d0, 912480445d0,   !   3986  1.00000
     . 42391181d0,-177476185d0,-602432635d0, 420232856d0, 917277575d0,   !   3987  1.00000
     . 42709742d0,-178219245d0,-606274184d0, 421533718d0, 922085738d0,   !   3988  1.00000
     . 43029203d0,-178963296d0,-610125064d0, 422836231d0, 926904946d0,   !   3989  1.00000
     . 43349566d0,-179708339d0,-613985284d0, 424140396d0, 931735212d0,   !   3990  1.00000
     . 43670830d0,-180454375d0,-617854856d0, 425446215d0, 936576550d0,   !   3991  1.00000
     . 43992999d0,-181201403d0,-621733792d0, 426753688d0, 941428972d0,   !   3992  1.00000
     . 44316072d0,-181949426d0,-625622102d0, 428062817d0, 946292491d0,   !   3993  1.00000
     . 44640050d0,-182698442d0,-629519797d0, 429373601d0, 951167122d0,   !   3994  1.00000
     . 44964936d0,-183448453d0,-633426890d0, 430686041d0, 956052876d0,   !   3995  1.00000
     . 45290729d0,-184199459d0,-637343390d0, 432000140d0, 960949767d0,   !   3996  1.00000
     . 45617432d0,-184951460d0,-641269310d0, 433315896d0, 965857808d0,   !   3997  1.00000
     . 45945045d0,-185704457d0,-645204661d0, 434633312d0, 970777013d0,   !   3998  1.00000
     . 46273569d0,-186458451d0,-649149453d0, 435952388d0, 975707393d0,   !   3999  1.00000
     . 46603006d0,-187213442d0,-653103699d0, 437273125d0, 980648963d0,   !   4000  1.00000
     . 46933356d0,-187969431d0,-657067409d0, 438595524d0, 985601736d0,   !   4001  1.00000
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4002
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4003
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4004
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4005
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4006
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4007
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4008
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4009
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4010
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4011
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4012
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4013
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4014
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4015
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4016
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4017
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4018
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4019
     .  0000000,  00000000,  00000000,  00000000,  00000000,   !   4020
     .  0000000,  00000000,  00000000,  00000000,  00000000/   !   4021
      end
 
 
c-------------------------------------------------------------------------
c
c
c
c
      block data fts0pk
      integer fts0p(41)
      common /fts0p_/fts0p
      data fts0p/
     .  0001, 0101, 0201, 0301, 0401,
     .  0601, 0701, 0801, 0901, 1001,
     .  1001, 1101, 1201, 1301, 1401,
     .  1601, 1701, 1801, 1901, 2001,
     .  2001, 2101, 2201, 2301, 2401,
     .  2601, 2701, 2801, 2901, 3001,
     .  3001, 3101, 3201, 3301, 3401,
     .  3601, 3701, 3801, 3901, 4001,
     .  4100/
      end
 
      block data fts0wk
      real*8 fts0(41,16)
      common /fts0_/fts0
      data fts0/
c ord=01
     .4978.00966d0,4980.35294d0,4982.69231d0,4985.02749d0,4987.35831d0,
     .4989.68449d0,4992.00561d0,4994.32113d0,4996.63057d0,4998.93345d0,
     .5001.22934d0,5003.51791d0,5005.79894d0,5008.07220d0,5010.33747d0,
     .5012.59456d0,5014.84328d0,5017.08346d0,5019.31492d0,5021.53749d0,
     .5023.75099d0,5025.95519d0,5028.14983d0,5030.33464d0,5032.50937d0,
     .5034.67374d0,5036.82756d0,5038.97067d0,5041.10294d0,5043.22422d0,
     .5045.33435d0,5047.43319d0,5049.52054d0,5051.59623d0,5053.66013d0,
     .5055.71210d0,5057.75199d0,5059.77965d0,5061.79496d0,5063.79780d0,
     .5065.78977d0,
c ord=02
     .5049.10901d0,5051.48647d0,5053.85995d0,5056.22917d0,5058.59397d0,
     .5060.95406d0,5063.30900d0,5065.65826d0,5068.00134d0,5070.33776d0,
     .5072.66708d0,5074.98898d0,5077.30321d0,5079.60956d0,5081.90781d0,
     .5084.19774d0,5086.47918d0,5088.75195d0,5091.01587d0,5093.27077d0,
     .5095.51646d0,5097.75271d0,5099.97926d0,5102.19584d0,5104.40218d0,
     .5106.59800d0,5108.78312d0,5110.95739d0,5113.12064d0,5115.27274d0,
     .5117.41353d0,5119.54285d0,5121.66052d0,5123.76637d0,5125.86024d0,
     .5127.94201d0,5130.01153d0,5132.06863d0,5134.11320d0,5136.14513d0,
     .5138.16602d0,
c ord=03
     .5122.26903d0,5124.68161d0,5127.09014d0,5129.49435d0,5131.89406d0,
     .5134.28898d0,5136.67868d0,5139.06261d0,5141.44025d0,5143.81114d0,
     .5146.17482d0,5148.53096d0,5150.87932d0,5153.21967d0,5155.55179d0,
     .5157.87548d0,5160.19055d0,5162.49681d0,5164.79409d0,5167.08220d0,
     .5169.36097d0,5171.63016d0,5173.88950d0,5176.13872d0,5178.37755d0,
     .5180.60570d0,5182.82299d0,5185.02926d0,5187.22436d0,5189.40813d0,
     .5191.58043d0,5193.74108d0,5195.88991d0,5198.02674d0,5200.15142d0,
     .5202.26382d0,5204.36378d0,5206.45114d0,5208.52579d0,5210.58760d0,
     .5212.63821d0,
c ord=04
     .5197.58204d0,5200.03070d0,5202.47526d0,5204.91542d0,5207.35101d0,
     .5209.78173d0,5212.20715d0,5214.62670d0,5217.03987d0,5219.44618d0,
     .5221.84517d0,5224.23650d0,5226.61993d0,5228.99522d0,5231.36217d0,
     .5233.72054d0,5236.07016d0,5238.41084d0,5240.74240d0,5243.06466d0,
     .5245.37742d0,5247.68046d0,5249.97350d0,5252.25627d0,5254.52848d0,
     .5256.78986d0,5259.04021d0,5261.27937d0,5263.50719d0,5265.72351d0,
     .5267.92818d0,5270.12104d0,5272.30189d0,5274.47056d0,5276.62690d0,
     .5278.77077d0,5280.90202d0,5283.02048d0,5285.12602d0,5287.21855d0,
     .5289.29969d0,
c ord=05
     .5275.14438d0,5277.63014d0,5280.11173d0,5282.58884d0,5285.06131d0,
     .5287.52883d0,5289.99096d0,5292.44714d0,5294.89683d0,5297.33955d0,
     .5299.77483d0,5302.20233d0,5304.62181d0,5307.03303d0,5309.43576d0,
     .5311.82979d0,5314.21493d0,5316.59099d0,5318.95779d0,5321.31514d0,
     .5323.66285d0,5326.00069d0,5328.32837d0,5330.64561d0,5332.95214d0,
     .5335.24767d0,5337.53201d0,5339.80498d0,5342.06643d0,5344.31621d0,
     .5346.55416d0,5348.78012d0,5350.99389d0,5353.19529d0,5355.38417d0,
     .5357.56039d0,5359.72379d0,5361.87421d0,5364.01153d0,5366.13562d0,
     .5368.24815d0,
c ord=06
     .5355.05856d0,5357.58247d0,5360.10213d0,5362.61723d0,5365.12762d0,
     .5367.63299d0,5370.13287d0,5372.62669d0,5375.11394d0,5377.59409d0,
     .5380.06669d0,5382.53139d0,5384.98793d0,5387.43608d0,5389.87561d0,
     .5392.30631d0,5394.72796d0,5397.14040d0,5399.54343d0,5401.93686d0,
     .5404.32050d0,5406.69411d0,5409.05741d0,5411.41010d0,5413.75192d0,
     .5416.08256d0,5418.40183d0,5420.70957d0,5423.00560d0,5425.28978d0,
     .5427.56195d0,5429.82194d0,5432.06954d0,5434.30459d0,5436.52693d0,
     .5438.73641d0,5440.93287d0,5443.11615d0,5445.28612d0,5447.44266d0,
     .5449.58746d0,
c ord=07
     .5437.43368d0,5439.99682d0,5442.55563d0,5445.10982d0,5447.65921d0,
     .5450.20349d0,5452.74219d0,5455.27474d0,5457.80060d0,5460.31925d0,
     .5462.83023d0,5465.33318d0,5467.82784d0,5470.31398d0,5472.79135d0,
     .5475.25975d0,5477.71896d0,5480.16881d0,5482.60910d0,5485.03964d0,
     .5487.46022d0,5489.87062d0,5492.27055d0,5494.65970d0,5497.03780d0,
     .5499.40455d0,5501.75975d0,5504.10323d0,5506.43482d0,5508.75438d0,
     .5511.06173d0,5513.35671d0,5515.63912d0,5517.90876d0,5520.16550d0,
     .5522.40918d0,5524.63963d0,5526.85670d0,5529.06025d0,5531.25016d0,
     .5533.42814d0,
c ord=08
     .5522.38446d0,5524.98795d0,5527.58705d0,5530.18144d0,5532.77095d0,
     .5535.35526d0,5537.93390d0,5540.50629d0,5543.07188d0,5545.63014d0,
     .5548.18061d0,5550.72291d0,5553.25678d0,5555.78199d0,5558.29831d0,
     .5560.80549d0,5563.30334d0,5565.79167d0,5568.27029d0,5570.73900d0,
     .5573.19760d0,5575.64585d0,5578.08345d0,5580.51011d0,5582.92554d0,
     .5585.32944d0,5587.72159d0,5590.10185d0,5592.47002d0,5594.82596d0,
     .5597.16951d0,5599.50049d0,5601.81869d0,5604.12393d0,5606.41605d0,
     .5608.69491d0,5610.96033d0,5613.21215d0,5615.45024d0,5617.67447d0,
     .5619.88658d0,
c ord=09
     .5610.03456d0,5612.67958d0,5615.32013d0,5617.95590d0,5620.58670d0,
     .5623.21221d0,5625.83196d0,5628.44534d0,5631.05181d0,5633.65084d0,
     .5636.24193d0,5638.82473d0,5641.39896d0,5643.96438d0,5646.52075d0,
     .5649.06785d0,5651.60546d0,5654.13339d0,5656.65146d0,5659.15945d0,
     .5661.65716d0,5664.14436d0,5666.62072d0,5669.08598d0,5671.53981d0,
     .5673.98193d0,5676.41212d0,5678.83021d0,5681.23602d0,5683.62940d0,
     .5686.01019d0,5688.37820d0,5690.73323d0,5693.07509d0,5695.40362d0,
     .5697.71867d0,5700.02007d0,5702.30764d0,5704.58127d0,5706.84081d0,
     .5709.08803d0,
c ord=10
     .5700.51378d0,5703.20156d0,5705.88479d0,5708.56315d0,5711.23645d0,
     .5713.90438d0,5716.56644d0,5719.22203d0,5721.87059d0,5724.51157d0,
     .5727.14449d0,5729.76898d0,5732.38475d0,5734.99157d0,5737.58919d0,
     .5740.17737d0,5742.75591d0,5745.32462d0,5747.88328d0,5750.43171d0,
     .5752.96969d0,5755.49697d0,5758.01325d0,5760.51823d0,5763.01160d0,
     .5765.49307d0,5767.96240d0,5770.41943d0,5772.86399d0,5775.29591d0,
     .5777.71504d0,5780.12117d0,5782.51411d0,5784.89366d0,5787.25966d0,
     .5789.61197d0,5791.95039d0,5794.27477d0,5796.58497d0,5798.88086d0,
     .5801.16421d0,
c ord=11
     .5793.96241d0,5796.69422d0,5799.42140d0,5802.14361d0,5804.86069d0,
     .5807.57230d0,5810.27793d0,5812.97699d0,5815.66889d0,5818.35308d0,
     .5821.02908d0,5823.69649d0,5826.35504d0,5829.00449d0,5831.64458d0,
     .5834.27507d0,5836.89576d0,5839.50644d0,5842.10692d0,5844.69698d0,
     .5847.27642d0,5849.84498d0,5852.40236d0,5854.94824d0,5857.48232d0,
     .5860.00430d0,5862.51394d0,5865.01107d0,5867.49552d0,5869.96713d0,
     .5872.42572d0,5874.87110d0,5877.30307d0,5879.72142d0,5882.12600d0,
     .5884.51666d0,5886.89321d0,5889.25547d0,5891.60332d0,5893.93662d0,
     .5896.25718d0,
c ord=12
     .5890.52843d0,5893.30560d0,5896.07806d0,5898.84546d0,5901.60764d0,
     .5904.36424d0,5907.11477d0,5909.85860d0,5912.59515d0,5915.32386d0,
     .5918.04423d0,5920.75587d0,5923.45849d0,5926.15185d0,5928.83569d0,
     .5931.50977d0,5934.17388d0,5936.82781d0,5939.47135d0,5942.10431d0,
     .5944.72646d0,5947.33755d0,5949.93726d0,5952.52528d0,5955.10129d0,
     .5957.66500d0,5960.21616d0,5962.75460d0,5965.28014d0,5967.79261d0,
     .5970.29185d0,5972.77766d0,5975.24983d0,5977.70816d0,5980.15248d0,
     .5982.58263d0,5984.99844d0,5987.39973d0,5989.78635d0,5992.15819d0,
     .5994.51706d0,
c ord=13
     .5990.37064d0,5993.19456d0,5996.01368d0,5998.82766d0,6001.63631d0,
     .6004.43929d0,6007.23609d0,6010.02607d0,6012.80864d0,6015.58323d0,
     .6018.34933d0,6021.10655d0,6023.85459d0,6026.59321d0,6029.32214d0,
     .6032.04114d0,6034.75000d0,6037.44849d0,6040.13642d0,6042.81358d0,
     .6045.47975d0,6048.13466d0,6050.77799d0,6053.40943d0,6056.02866d0,
     .6058.63537d0,6061.22931d0,6063.81032d0,6066.37820d0,6068.93279d0,
     .6071.47392d0,6074.00139d0,6076.51498d0,6079.01449d0,6081.49975d0,
     .6083.97060d0,6086.42687d0,6088.86836d0,6091.29494d0,6093.70647d0,
     .6096.10481d0,
c ord=14
     .6093.65896d0,6096.53109d0,6099.39833d0,6102.26032d0,6105.11690d0,
     .6107.96769d0,6110.81219d0,6113.64976d0,6116.47977d0,6119.30167d0,
     .6122.11492d0,6124.91913d0,6127.71401d0,6130.49928d0,6133.27471d0,
     .6136.04002d0,6138.79501d0,6141.53945d0,6144.27314d0,6146.99587d0,
     .6149.70742d0,6152.40750d0,6155.09581d0,6157.77202d0,6160.43580d0,
     .6163.08683d0,6165.72488d0,6168.34976d0,6170.96130d0,6173.55930d0,
     .6176.14361d0,6178.71402d0,6181.27031d0,6183.81227d0,6186.33973d0,
     .6188.85254d0,6191.35050d0,6193.83343d0,6196.30118d0,6198.75363d0,
     .6201.19265d0,
c ord=15
     .6200.57486d0,6203.49671d0,6206.41358d0,6209.32511d0,6212.23112d0,
     .6215.13124d0,6218.02495d0,6220.91160d0,6223.79056d0,6226.66125d0,
     .6229.52314d0,6232.37582d0,6235.21900d0,6238.05241d0,6240.87579d0,
     .6243.68887d0,6246.49144d0,6249.28328d0,6252.06418d0,6254.83391d0,
     .6257.59225d0,6260.33894d0,6263.07363d0,6265.79601d0,6268.50574d0,
     .6271.20249d0,6273.88603d0,6276.55616d0,6279.21271d0,6281.85549d0,
     .6284.48433d0,6287.09902d0,6289.69934d0,6292.28508d0,6294.85606d0,
     .6297.41213d0,6299.95308d0,6302.47874d0,6304.98896d0,6307.48360d0,
     .6309.96457d0,
c ord=16
     .6311.31427d0,6314.28744d0,6317.25553d0,6320.21818d0,6323.17521d0,
     .6326.12623d0,6329.07072d0,6332.00802d0,6334.93749d0,6337.85853d0,
     .6340.77061d0,6343.67332d0,6346.56634d0,6349.44942d0,6352.32227d0,
     .6355.18465d0,6358.03632d0,6360.87706d0,6363.70665d0,6366.52489d0,
     .6369.33152d0,6372.12629d0,6374.90884d0,6377.67885d0,6380.43598d0,
     .6383.17991d0,6385.91037d0,6388.62720d0,6391.33018d0,6394.01916d0,
     .6396.69393d0,6399.35431d0,6399.99881d0,6399.99881d0,6399.99881d0,
     .6399.99881d0,6399.99881d0,6399.99881d0,6399.99881d0,6399.99881d0,
     .6399.99881d0/
      end
 
 

      block data le_jq_bdata              
      integer*2  le_jq_x10000(13,201)     
      common /le_jq_x10000_/le_jq_x10000
      data le_jq_x10000/                  
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 001  x=-10.0
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 002  x= -9.9
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 003  x= -9.8
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 004  x= -9.7
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 005  x= -9.6
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 006  x= -9.5
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 007  x= -9.4
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 008  x= -9.3
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 009  x= -9.2
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 010  x= -9.1
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 011  x= -9.0
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 012  x= -8.9
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 013  x= -8.8
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 014  x= -8.7
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 015  x= -8.6
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 016  x= -8.5
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 017  x= -8.4
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 018  x= -8.3
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 019  x= -8.2
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 020  x= -8.1
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 021  x= -8.0
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 022  x= -7.9
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 023  x= -7.8
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 024  x= -7.7
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 025  x= -7.6
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 026  x= -7.5
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 027  x= -7.4
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 028  x= -7.3
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 029  x= -7.2
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 030  x= -7.1
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 031  x= -7.0
     . 1031,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 032  x= -6.9
     . 2063,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 033  x= -6.8
     . 3094,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 034  x= -6.7
     . 4126,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 035  x= -6.6
     . 5158,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 036  x= -6.5
     . 6189,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 037  x= -6.4
     . 7221,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 038  x= -6.3
     . 8252,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 039  x= -6.2
     . 9284,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 040  x= -6.1
     .10315,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 041  x= -6.0
     . 9284, 1002,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 042  x= -5.9
     . 8252, 2004,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 043  x= -5.8
     . 7221, 3006,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 044  x= -5.7
     . 6189, 4008,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 045  x= -5.6
     . 5158, 5010,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 046  x= -5.5
     . 4126, 6012,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 047  x= -5.4
     . 3095, 7014,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 048  x= -5.3
     . 2063, 8016,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 049  x= -5.2
     . 1032, 9018,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 050  x= -5.1
     .    0,10020,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 051  x= -5.0
     .  -83, 9175,  919,   -3,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 052  x= -4.9
     . -263, 8497, 1790,  -30,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 053  x= -4.8
     . -454, 7803, 2720,  -91,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 054  x= -4.7
     . -593, 6974, 3769, -187,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 055  x= -4.6
     . -644, 5949, 4949, -305,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 056  x= -4.5
     . -593, 4729, 6224, -422,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 057  x= -4.4
     . -454, 3374, 7511, -502,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 058  x= -4.3
     . -263, 2004, 8679, -500,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 059  x= -4.2
     .  -83,  799, 9550, -356,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 060  x= -4.1
     .    0,    0, 9898,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 061  x= -4.0
     .    7, -419, 9667,  683,  -52,    1,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 062  x= -3.9
     .   22, -688, 9038, 1685, -190,    5,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 063  x= -3.8
     .   38, -825, 8103, 2904, -381,   16,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 064  x= -3.7
     .   50, -849, 6953, 4242, -589,   32,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 065  x= -3.6
     .   54, -782, 5671, 5605, -778,   52,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 066  x= -3.5
     .   50, -648, 4339, 6903, -908,   72,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 067  x= -3.4
     .   38, -474, 3035, 8050, -940,   86,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 068  x= -3.3
     .   22, -288, 1832, 8964, -829,   86,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 069  x= -3.2
     .    7, -118,  798, 9566, -531,   61,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 070  x= -3.1
     .    0,    0,    0, 9783,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 071  x= -3.0
     .    0,   68, -541, 9567,  818, -127,
     .                                  7,
     .                                  0,    0,    0,    0,    0,    0, ! 072  x= -2.9
     .    0,  107, -870, 8974, 1888, -320,
     .                                 27,
     .                                  0,    0,    0,    0,    0,    0, ! 073  x= -2.8
     .    0,  123,-1021, 8080, 3131, -543,
     .                                 52,
     .                                  0,    0,    0,    0,    0,    0, ! 074  x= -2.7
     .    0,  120,-1028, 6965, 4468, -761,
     .                                 80,
     .                                  0,    0,    0,    0,    0,    0, ! 075  x= -2.6
     .    0,  104, -927, 5707, 5817, -939,
     .                                104,
     .                                  0,    0,    0,    0,    0,    0, ! 076  x= -2.5
     .    0,   80, -751, 4383, 7100,-1042,
     .                                120,
     .                                  0,    0,    0,    0,    0,    0, ! 077  x= -2.4
     .    0,   53, -536, 3072, 8237,-1034,
     .                                122,
     .                                  0,    0,    0,    0,    0,    0, ! 078  x= -2.3
     .    0,   27, -316, 1852, 9148, -881,
     .                                107,
     .                                  0,    0,    0,    0,    0,    0, ! 079  x= -2.2
     .    0,    8, -125,  802, 9753, -548,
     .                                 67,
     .                                  0,    0,    0,    0,    0,    0, ! 080  x= -2.1
     .    0,    0,    0,    0, 9972,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 081  x= -2.0
     .    0,    0,   67, -535, 9753,  822,
     .                               -126,
     .                                  8,    0,    0,    0,    0,    0, ! 082  x= -1.9
     .    0,    0,  106, -860, 9148, 1899,
     .                               -319,
     .                                 27,    0,    0,    0,    0,    0, ! 083  x= -1.8
     .    0,    0,  121,-1009, 8237, 3149,
     .                               -541,
     .                                 53,    0,    0,    0,    0,    0, ! 084  x= -1.7
     .    0,    0,  119,-1016, 7100, 4492,
     .                               -759,
     .                                 80,    0,    0,    0,    0,    0, ! 085  x= -1.6
     .    0,    0,  103, -916, 5817, 5850,
     .                               -937,
     .                                104,    0,    0,    0,    0,    0, ! 086  x= -1.5
     .    0,    0,   79, -742, 4468, 7140,
     .                              -1039,
     .                                120,    0,    0,    0,    0,    0, ! 087  x= -1.4
     .    0,    0,   52, -530, 3131, 8283,
     .                              -1032,
     .                                123,    0,    0,    0,    0,    0, ! 088  x= -1.3
     .    0,    0,   26, -312, 1888, 9199,
     .                               -879,
     .                                107,    0,    0,    0,    0,    0, ! 089  x= -1.2
     .    0,    0,    7, -124,  818, 9807,
     .                               -546,
     .                                 68,    0,    0,    0,    0,    0, ! 090  x= -1.1
     .    0,    0,    0,    0,    0,10028,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 091  x= -1.0
     .    0,    0,    0,   66, -545, 9807,
     .                                820,
     .                               -127,    7,    0,    0,    0,    0, ! 092  x= -0.9
     .    0,    0,    0,  104, -877, 9199,
     .                               1893,
     .                               -320,   27,    0,    0,    0,    0, ! 093  x= -0.8
     .    0,    0,    0,  120,-1029, 8283,
     .                               3140,
     .                               -543,   52,    0,    0,    0,    0, ! 094  x= -0.7
     .    0,    0,    0,  117,-1036, 7140,
     .                               4480,
     .                               -761,   80,    0,    0,    0,    0, ! 095  x= -0.6
     .    0,    0,    0,  102, -934, 5850,
     .                               5833,
     .                               -939,  104,    0,    0,    0,    0, ! 096  x= -0.5
     .    0,    0,    0,   78, -757, 4492,
     .                               7120,
     .                              -1042,  120,    0,    0,    0,    0, ! 097  x= -0.4
     .    0,    0,    0,   51, -540, 3149,
     .                               8260,
     .                              -1034,  122,    0,    0,    0,    0, ! 098  x= -0.3
     .    0,    0,    0,   26, -318, 1899,
     .                               9173,
     .                               -881,  106,    0,    0,    0,    0, ! 099  x= -0.2
     .    0,    0,    0,    7, -126,  822,
     .                               9780,
     .                               -548,   67,    0,    0,    0,    0, ! 100  x= -0.1
     .    0,    0,    0,    0,    0,    0,
     .                              10000,
     .                                  0,    0,    0,    0,    0,    0, ! 101  x=  0.0
     .    0,    0,    0,    0,   67, -548,
     .                               9780,
     .                                822, -126,    7,    0,    0,    0, ! 102  x=  0.1
     .    0,    0,    0,    0,  106, -881,
     .                               9173,
     .                               1899, -318,   26,    0,    0,    0, ! 103  x=  0.2
     .    0,    0,    0,    0,  122,-1034,
     .                               8260,
     .                               3149, -540,   51,    0,    0,    0, ! 104  x=  0.3
     .    0,    0,    0,    0,  120,-1042,
     .                               7120,
     .                               4492, -757,   78,    0,    0,    0, ! 105  x=  0.4
     .    0,    0,    0,    0,  104, -939,
     .                               5833,
     .                               5850, -934,  102,    0,    0,    0, ! 106  x=  0.5
     .    0,    0,    0,    0,   80, -761,
     .                               4480,
     .                               7140,-1036,  117,    0,    0,    0, ! 107  x=  0.6
     .    0,    0,    0,    0,   52, -543,
     .                               3140,
     .                               8283,-1029,  120,    0,    0,    0, ! 108  x=  0.7
     .    0,    0,    0,    0,   27, -320,
     .                               1893,
     .                               9199, -877,  104,    0,    0,    0, ! 109  x=  0.8
     .    0,    0,    0,    0,    7, -127,
     .                                820,
     .                               9807, -545,   66,    0,    0,    0, ! 110  x=  0.9
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                              10028,    0,    0,    0,    0,    0, ! 111  x=  1.0
     .    0,    0,    0,    0,    0,   68,
     .                               -546,
     .                               9807,  818, -124,    7,    0,    0, ! 112  x=  1.1
     .    0,    0,    0,    0,    0,  107,
     .                               -879,
     .                               9199, 1888, -312,   26,    0,    0, ! 113  x=  1.2
     .    0,    0,    0,    0,    0,  123,
     .                              -1032,
     .                               8283, 3131, -530,   52,    0,    0, ! 114  x=  1.3
     .    0,    0,    0,    0,    0,  120,
     .                              -1039,
     .                               7140, 4468, -742,   79,    0,    0, ! 115  x=  1.4
     .    0,    0,    0,    0,    0,  104,
     .                               -936,
     .                               5850, 5817, -916,  103,    0,    0, ! 116  x=  1.5
     .    0,    0,    0,    0,    0,   80,
     .                               -759,
     .                               4492, 7100,-1016,  119,    0,    0, ! 117  x=  1.6
     .    0,    0,    0,    0,    0,   53,
     .                               -541,
     .                               3149, 8237,-1009,  121,    0,    0, ! 118  x=  1.7
     .    0,    0,    0,    0,    0,   27,
     .                               -319,
     .                               1899, 9148, -860,  106,    0,    0, ! 119  x=  1.8
     .    0,    0,    0,    0,    0,    8,
     .                               -126,
     .                                822, 9753, -535,   67,    0,    0, ! 120  x=  1.9
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0, 9972,    0,    0,    0,    0, ! 121  x=  2.0
     .    0,    0,    0,    0,    0,    0,
     .                                 67,
     .                               -548, 9753,  802, -125,    8,    0, ! 122  x=  2.1
     .    0,    0,    0,    0,    0,    0,
     .                                107,
     .                               -881, 9148, 1852, -316,   27,    0, ! 123  x=  2.2
     .    0,    0,    0,    0,    0,    0,
     .                                122,
     .                              -1034, 8237, 3072, -536,   53,    0, ! 124  x=  2.3
     .    0,    0,    0,    0,    0,    0,
     .                                120,
     .                              -1042, 7100, 4383, -751,   80,    0, ! 125  x=  2.4
     .    0,    0,    0,    0,    0,    0,
     .                                104,
     .                               -939, 5817, 5707, -927,  104,    0, ! 126  x=  2.5
     .    0,    0,    0,    0,    0,    0,
     .                                 80,
     .                               -761, 4468, 6965,-1028,  120,    0, ! 127  x=  2.6
     .    0,    0,    0,    0,    0,    0,
     .                                 52,
     .                               -543, 3131, 8080,-1021,  123,    0, ! 128  x=  2.7
     .    0,    0,    0,    0,    0,    0,
     .                                 27,
     .                               -320, 1888, 8974, -870,  107,    0, ! 129  x=  2.8
     .    0,    0,    0,    0,    0,    0,
     .                                  7,
     .                               -127,  818, 9567, -541,   68,    0, ! 130  x=  2.9
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0, 9783,    0,    0,    0, ! 131  x=  3.0
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                 61, -531, 9566,  798, -118,    7, ! 132  x=  3.1
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                 86, -829, 8964, 1832, -288,   22, ! 133  x=  3.2
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                 86, -940, 8050, 3035, -474,   38, ! 134  x=  3.3
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                 72, -908, 6903, 4339, -648,   50, ! 135  x=  3.4
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                 52, -778, 5605, 5671, -782,   54, ! 136  x=  3.5
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                 32, -589, 4242, 6953, -849,   50, ! 137  x=  3.6
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                 16, -381, 2904, 8103, -825,   38, ! 138  x=  3.7
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  5, -190, 1685, 9038, -688,   22, ! 139  x=  3.8
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  1,  -52,  683, 9667, -419,    7, ! 140  x=  3.9
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0, 9898,    0,    0, ! 141  x=  4.0
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0, -356, 9550,  799,  -83, ! 142  x=  4.1
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0, -500, 8679, 2004, -263, ! 143  x=  4.2
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0, -502, 7511, 3374, -454, ! 144  x=  4.3
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0, -422, 6224, 4729, -593, ! 145  x=  4.4
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0, -305, 4949, 5949, -644, ! 146  x=  4.5
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0, -187, 3769, 6974, -593, ! 147  x=  4.6
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,  -91, 2720, 7803, -454, ! 148  x=  4.7
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,  -30, 1790, 8497, -263, ! 149  x=  4.8
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,   -3,  919, 9175,  -83, ! 150  x=  4.9
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,10020,    0, ! 151  x=  5.0
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0, 9018, 1032, ! 152  x=  5.1
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0, 8016, 2063, ! 153  x=  5.2
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0, 7014, 3095, ! 154  x=  5.3
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0, 6012, 4126, ! 155  x=  5.4
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0, 5010, 5158, ! 156  x=  5.5
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0, 4008, 6189, ! 157  x=  5.6
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0, 3006, 7221, ! 158  x=  5.7
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0, 2004, 8252, ! 159  x=  5.8
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0, 1002, 9284, ! 160  x=  5.9
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,10315, ! 161  x=  6.0
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0, 9284, ! 162  x=  6.1
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0, 8252, ! 163  x=  6.2
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0, 7221, ! 164  x=  6.3
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0, 6189, ! 165  x=  6.4
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0, 5158, ! 166  x=  6.5
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0, 4126, ! 167  x=  6.6
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0, 3094, ! 168  x=  6.7
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0, 2063, ! 169  x=  6.8
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0, 1031, ! 170  x=  6.9
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 171  x=  7.0
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 172  x=  7.1
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 173  x=  7.2
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 174  x=  7.3
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 175  x=  7.4
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 176  x=  7.5
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 177  x=  7.6
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 178  x=  7.7
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 179  x=  7.8
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 180  x=  7.9
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 181  x=  8.0
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 182  x=  8.1
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 183  x=  8.2
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 184  x=  8.3
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 185  x=  8.4
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 186  x=  8.5
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 187  x=  8.6
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 188  x=  8.7
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 189  x=  8.8
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 190  x=  8.9
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 191  x=  9.0
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 192  x=  9.1
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 193  x=  9.2
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 194  x=  9.3
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 195  x=  9.4
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 196  x=  9.5
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 197  x=  9.6
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 198  x=  9.7
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 199  x=  9.8
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0, ! 200  x=  9.9
     .    0,    0,    0,    0,    0,    0,
     .                                  0,
     .                                  0,    0,    0,    0,    0,    0/ ! 201  x= 10.0
      end



      subroutine block_meteor(RDSKFILE,uuu,METEORFILE)
      implicit none

      character*80 RDSKFILE
      character*200 METEORFILE
      integer*1 uuu(_IDIM_,_MDIM_)

      character*80 RDSKFILE_RJ

      character*10 obsid
      real y7at2k
      integer meteor(16)
      real zmed
      real mar
  
      integer m
      integer i

      integer A, B, AA, BB
      logical match

    
      write(6,*) 'BLOCK_METEOR:  METEORFILE: ',METEORFILE 
      if (METEORFILE(1:4).eq.'NONE') return

      RDSKFILE_RJ = 'NULL'
      do i = 1, 80
         if (RDSKFILE(i:i).eq.'.') RDSKFILE(i:i) = ' '
         enddo
      do i = 1, 75
         if (RDSKFILE(i:i+1).eq.'rj') RDSKFILE_RJ = RDSKFILE(i:80)
         enddo
      read(RDSKFILE_RJ(3:80),*) A,B 
   
      open(37,file=METEORFILE,status='old')
    1 read(37,*,end=2) obsid, y7at2k, (meteor(m),m=1,16), zmed, mar
      do i = 1, 10
         if (OBSID(i:i).eq.'.') OBSID(i:i) = ' '
         enddo
      read(OBSID(3:10),*) AA,BB 
      if (A.eq.AA.and.B.eq.BB) then
         write(6,*) 'RJFILE: ',RDSKFILE_RJ(1:40)
         write(6,*) 'BKFILE: ',OBSID
         do m = 02, 15
            do i = 0001, _IDIM_
               uuu(i,m) = 0
               enddo
            do i = (meteor(m+1)+meteor(m  ))/2,
     .             (meteor(m  )+meteor(m-1))/2
               uuu(i,m) = 1
               enddo 
            do i = meteor(m)-12, meteor(m)+12
               uuu(i,m) = 2
               enddo
            i = meteor(m)
            uuu(i,m) = 9  !  center of the meteor
            enddo 
         close(37)
         return
         endif
      goto 1
    2 continue
      close(37)
      write(6,*) 'could not find match in meteor file...'
      write(6,*) 'FILEU   : ',RDSKFILE
      write(6,*) 'FILEU_RJ: ',RDSKFILE_RJ
      stop 

      return
      end




c-----------------------------------------------------
c
c
c
      subroutine smoo_tem(tems,NIT)
      implicit none
   
      real    tems(_SDIM_)   ! current template
      integer NIT

      real    tem0(_SDIM_)
      real    temz(_SDIM_) 

      real    temb(_SDIM_)
      real    temd(_SDIM_)
      real    temx(_SDIM_)

      integer s
      real*8  ww_s
      character*80 FILENAME

      do s = 0001, _SDIM_
         temz(s) = tems(s)
         tem0(s) = tems(s)
         enddo

      do s = 0003, _SDIM_-2
         temb(s) =  0.333*tem0(s-1)
     .           +  0.334*tem0(s  )
     .           +  0.333*tem0(s+1)
         temd(s) = tem0(s)-temb(s)
         enddo

c 
c sextic 9-point (7 terms, up to x**6)
c        
      do s = 0005, _SDIM_-4
         temx(s) = 
     .           -  0.005439*temd(s-4)  
     .           +  0.043512*temd(s-3)
     .           -  0.152292*temd(s-2)
     .           +  0.304584*temd(s-1)
     .           +  0.619270*temd(s  )
     .           +  0.304584*temd(s+1)
     .           -  0.152292*temd(s+2)
     .           +  0.043512*temd(s+3)
     .           -  0.005439*temd(s+4)
         enddo

      do s = 0001, _SDIM_
         tems(s) = temb(s) + temx(s)
         enddo

      write(FILENAME,'(''LOG.SMOO_TEM_NIT'',i2.2)') NIT 
      open(92,file=FILENAME,status='unknown')
      do s = 0001, _SDIM_
         write(92,192) s,ww_s(s),tem0(s),tems(s),
     .                           temb(s),temd(s), temx(s)
  192    format(i6.6,1x,f9.4,1x,8(1x,f8.5))
         enddo
      close(92)

      return
      end






c----------------------------------------------------------------------
c
c this is a new version of the original wavelength solution function...
c
c
      real*8 function wav_rim_frompoly(ri,m)
      implicit none

      real    ri
      integer m

      real*8 wav5m(5,16)
      common /wav5m_/wav5m
      data wav5m/
c        
c      i=0001       i=1001       i=2001       i=3001       i=4001
c   
c      x=-1.0       x=-0.5       x= 0.0       x= 0.5       x= 1.0
c   
     . 4978.00703d0,5001.22774d0,5023.75220d0,5045.33407d0,5065.78966d0, ! m=01
     . 5049.10635d0,5072.66547d0,5095.51769d0,5117.41324d0,5138.16591d0, ! m=02
     . 5122.26633d0,5146.17318d0,5169.36222d0,5191.58013d0,5212.63810d0, ! m=03
     . 5197.57929d0,5221.84351d0,5245.37869d0,5267.92789d0,5289.29958d0, ! m=04
     . 5275.14159d0,5299.77314d0,5323.66413d0,5346.55386d0,5368.24805d0, ! m=05
     . 5355.05573d0,5380.06498d0,5404.32180d0,5427.56165d0,5449.58736d0, ! m=06
     . 5437.43081d0,5462.82849d0,5487.46155d0,5511.06142d0,5533.42804d0, ! m=07
     . 5522.38154d0,5548.17884d0,5573.19894d0,5597.16920d0,5619.88648d0, ! m=08
     . 5610.03159d0,5636.24014d0,5661.65853d0,5686.00987d0,5709.08793d0, ! m=09
     . 5700.51077d0,5727.14267d0,5752.97107d0,5777.71471d0,5801.16412d0, ! m=10
     . 5793.95935d0,5821.02722d0,5847.27783d0,5872.42539d0,5896.25708d0, ! m=11
     . 5890.52532d0,5918.04234d0,5944.72789d0,5970.29152d0,5994.51696d0, ! m=12
     . 5990.36747d0,6018.34741d0,6045.48120d0,6071.47358d0,6096.10472d0, ! m=13
     . 6093.65574d0,6122.11297d0,6149.70890d0,6176.14327d0,6201.19256d0, ! m=14
     . 6200.57158d0,6229.52116d0,6257.59376d0,6284.48398d0,6309.96448d0, ! m=15
     . 6311.08441d0,6340.56336d0,6369.67956d0,6395.22480d0,6398.17334d0/ ! m=16
c        
c      i=0001       i=1001       i=2001       i=3001       i=4001
c   
c      x=-1.0       x=-0.5       x= 0.0       x= 0.5       x= 1.0
c   


      real*8 matrix55(5,5)
      common /matrix55_/matrix55  
      data matrix55/  
     . 0.0000000d0, 0.0000000d0, 1.0000000d0, 0.0000000d0, 0.0000000d0,
     . 0.1666667d0,-1.3333333d0, 0.0000000d0, 1.3333333d0,-0.1666667d0,
     .-0.1666667d0, 2.6666667d0,-5.0000000d0, 2.6666667d0,-0.1666667d0,
     .-0.6666667d0, 1.3333333d0, 0.0000000d0,-1.3333333d0, 0.6666667d0,
     . 0.6666667d0,-2.6666667d0, 4.0000000d0,-2.6666667d0, 0.6666667d0/

      real*8 APOLY, BPOLY, CPOLY, DPOLY, EPOLY
      real*8 f(5)
      real*8 A, B, C, D, E

      integer i
      real*8 x

      f(1) = wav5m(1,m)
      f(2) = wav5m(2,m)
      f(3) = wav5m(3,m)
      f(4) = wav5m(4,m)
      f(5) = wav5m(5,m)

      APOLY = 0.
      BPOLY = 0.
      CPOLY = 0.
      DPOLY = 0.
      EPOLY = 0.
      do i = 1, 5
         APOLY = APOLY + f(i)*matrix55(i,1)
         BPOLY = BPOLY + f(i)*matrix55(i,2)
         CPOLY = CPOLY + f(i)*matrix55(i,3)
         DPOLY = DPOLY + f(i)*matrix55(i,4)
         EPOLY = EPOLY + f(i)*matrix55(i,5)
         enddo

      x = (ri-2001.0)/2000.0
 
      wav_rim_frompoly = APOLY
     .                 + BPOLY*x
     .                 + CPOLY*x**2
     .                 + DPOLY*x**3
     .                 + EPOLY*x**4

      return
      end 



     
c----------------------------------------------
c
c 
      subroutine cnorm_tem_bar(tems,temo,NIT,LABEL59,w1,w2)
      implicit none
 
      real   tems(_SDIM_)
      real   temo(_SDIM_)
      integer     NIT
      character*59 LABEL59
      real*8 w1, w2

      real   tem2(_SDIM_)
      real   tem3(_SDIM_)
      real   tem95
 
      real   temh(_SDIM_)
      real   teml(_SDIM_)
 
      integer s, ss ,s1, s2
      character*80 FILENAME

      integer NITu
      integer OUTPUT_NRM
      common /OUTPUT_NRM_/OUTPUT_NRM

      integer ntot
      real*8  ptot, pnorm
      real*8  ww_s
 
      write(6,*) 'enter cnorm_tem_bar...'

      ntot = 0
      ptot = 0. 
      do s = 0001, _SDIM_
         if (ww_s(s).ge.w1.and.ww_s(s).le.w2) then
            ntot = ntot + 1
            ptot = ptot + tems(s)
            endif
         enddo

      pnorm = ptot/ntot

      do s = 0001, _SDIM_
         tems(s) = tems(s)/pnorm
         enddo

      return
      end




c------------------------------------------------------
c
c
c
      subroutine find_twktem(tems,
     .                       pix_imn,mod_imn,
     .                       rat_imn,bar_imn,uuu_imn,
     .                       wav_imn,zn,unim,
     .                       NIMs,MMMM)
      implicit none
 
      real    tems(_SDIM_)
      real       pix_imn(_IDIM_,_MDIM_,NIMs)
      real       mod_imn(_IDIM_,_MDIM_,NIMs)
      real       rat_imn(_IDIM_,_MDIM_,NIMs)
      real       bar_imn(_IDIM_,_MDIM_,NIMs)
      integer*1  uuu_imn(_IDIM_,_MDIM_,NIMs)
      real*8     wav_imn(_IDIM_,_MDIM_,NIMs)
      real*8     zn(NIMs)
      logical    unim(NIMs)
      integer    NIMs
      integer    MMMM
 
      integer s
      integer s1, s2, su
      real    teml(_SDIM_)
      real    temh(_SDIM_)
      real*8  ww_s
      integer u
      integer sau, sbu
      real*8  wmin, wmax
      integer i, m, n
      real*8  r8s_ww
      integer ss
      real*8  rl(9999), bar, sig
      integer Ls, Lu
      real*8  wwa, wwb, www, wwn 
      integer ic, i_w

      wmin = 9e9
      wmax = 0e0
      do i = 0001, 4020
      do m = MMMM, MMMM
      do n = 0001, NIMs
         if (unim(n)) then
            if (wav_imn(i,m,n).gt.wmax) wmax = wav_imn(i,m,n)
            if (wav_imn(i,m,n).lt.wmin) wmin = wav_imn(i,m,n)
            endif
         enddo
         enddo
         enddo

      sau = int(r8s_ww(wmin)+0.5)+1
      sbu = int(r8s_ww(wmax)+0.5)-1

      do s = sau, sbu, 1
         wwa = ww_s(s-1)
         www = ww_s(s)
         wwb = ww_s(s+1)
         Ls = 0
         do n = 0001, NIMs
            do M = MMMM, MMMM
               ic = i_w(www*(1+zn(n)),wav_imn(1,m,n))
               wwn = wav_imn(ic,m,n)/(1+zn(n))
               do i = max(0001,ic-2), min(4020,ic+2)
                  wwn = wav_imn(i,m,n)/(1+zn(n))
                  if (wwn.ge.wwa.and.wwn.le.wwb) then
                      Ls = Ls + 1
                      rl(Ls) = (pix_imn(i,m,n)
     .                             -mod_imn(i,m,n)*rat_imn(i,m,n))/
     .                                     bar_imn(i,m,n)
                      endif
                  enddo
               enddo
            enddo
         call barsig(rl,Ls,bar,sig,Lu)
         write(47,147) s,www,wwa,wwb,bar,sig,Ls,Lu,tems(s)
  147    format(1x,i6.6,3x,f12.6,1x,f12.6,1x,f12.6,
     .                  3x,f9.5,1x,f9.5,1x,i3,1x,i3,
     .                  3x,f9.5)
         tems(s) = tems(s) + bar
         enddo
      do s = 00001, sau-1
         tems(s) = tems(sau)
         enddo
      do s = sbu, _SDIM_
         tems(s) = tems(sbu)
         enddo
      close(47)

      return
      end
 
 
