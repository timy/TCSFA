#include '../include/inc_misc.h'
#ifdef MISC_PLOT
#include '../include/inc_plot_crf_map.h'
#endif
#ifdef  CRF_PLOT_MAP_FIX
#include '../include/inc_ts_guess.h'
#endif
!#define CRF_PLOT_MAP_VAR
!#define CRF_PLOT_ARRY_MAX 4000
!#define CRF_PLOT_MAP_X_FACTOR 2d0
!#define CRF_PLOT_MAP_Y_FACTOR 2d0
!#define CRF_PLOT_MAP_X_OFFSET 0d0
!#define CRF_PLOT_MAP_Y_OFFSET 0d0


!#define CRF_PLOT_FILE_NAME1 'dat/crf_track.dat'
!#define CRF_PRNT_MAP_INFO

!#define CRF_PRNT_DATA

!#define CRF_NORM_ADD
    
    
    ! --------------------------------------------------------------------------------    
    
    double complex function find_root( zs, func, hs, ierr )

        implicit none
        double complex:: zs;
        double precision, intent(in):: hs;
        integer:: ierr;

        double complex:: ze;
        double precision:: hm, dm, ds, he, de
        integer:: n;
        double complex, external:: func;
        !integer, parameter:: n_arr_max = CRF_PLOT_ARRY_MAX
        !double precision:: zh_re(n_arr_max), zh_im(n_arr_max), wh(n_arr_max);
#ifdef CRF_PLOT_MAP_VAR
        double precision:: re_max, re_min, im_max, im_min;
#endif
#ifdef CRF_PLOT_MAP
        double precision:: re_hi, re_lo, im_hi, im_lo, d_re, d_im
        integer, parameter:: n_grid_x = CRF_PLOT_N_GRID_X
        integer, parameter:: n_grid_y = CRF_PLOT_N_GRID_Y
        double complex:: grid(n_grid_x, n_grid_y), cw
        double precision:: grid_re(n_grid_x), grid_im(n_grid_y);
#endif
#ifdef CRF_PLOT_MAP_VAR
        character(len=200):: line;
#endif

        double precision, external:: get_tp

        hm = 1d-9;
        dm = 1d-7;

#ifdef CRF_PLOT_MAP_FIX
        re_lo = LMS_RE_LOWER; re_hi = LMS_RE_UPPER; 
        im_lo = LMS_IM_LOWER; im_hi = LMS_IM_UPPER;
        call crf_plot_map( n_grid_x, n_grid_y, re_lo, re_hi, im_lo, im_hi, func )
#endif


        call crf ( zs, hs, hm, dm, func, ds, ze, he, de, n );

        if( de > dm ) then
#ifdef CRF_PRNT_INFO
            print*, 'searching failed!'
#endif
            ierr = 1;
        end if


#ifdef CRF_PRNT_DATA
            print*, 'zs:', zs;
            print*, 'hs:', hs;
            print*, 'hm:', hm;
            print*, 'dm:', dm;
            print*, 'ds:', ds;
            print*, 'ze:', ze;
            print*, 'he:', he;
            print*, 'de:', de;
            print*, 'n:', n;
#endif


        find_root = ze;
        ierr = 0;

#ifdef CRF_PLOT_MAP_VAR

        ! open the history file

        open(CRF_PLOT_FILE_ID, file=CRF_PLOT_FILE_NAME1);
        n = 0;
        do
            read (CRF_PLOT_FILE_ID, '(a)', iostat = io) line;
            if (io < 0) exit;                 ! there is no data
            if (len_trim(line) == 0) cycle;   ! if it is a blank line
            n = n + 1;                        ! it is a new line
            if( n > n_arr_max ) stop 'exceed n_arr_max when processing plot file!'
            read (line, '(3(e15.8,2x))'), zh_re(n), zh_im(n), wh(n);
        end do
        close(CRF_PLOT_FILE_ID);

        if( n == 0 ) stop 'Forget to generate crf_track.dat first?'

        ! get the boundary
        re_max = maxval(zh_re(1:n)); re_min = minval(zh_re(1:n));
        im_max = maxval(zh_im(1:n)); im_min = minval(zh_im(1:n));
        re_lo = re_min - (re_max - re_min) * CRF_PLOT_MAP_X_FACTOR - CRF_PLOT_MAP_X_OFFSET;
        re_hi = re_max + (re_max - re_min) * CRF_PLOT_MAP_X_FACTOR - CRF_PLOT_MAP_X_OFFSET;
        im_lo = im_min - (im_max - im_min) * CRF_PLOT_MAP_Y_FACTOR - CRF_PLOT_MAP_Y_OFFSET;
        im_hi = im_max + (im_max - im_min) * CRF_PLOT_MAP_Y_FACTOR - CRF_PLOT_MAP_Y_OFFSET;


#ifdef CRF_PRNT_MAP_INFO
        print*, 're_max', re_max;
        print*, 're_min', re_min;
        print*, 'im_max', im_max;
        print*, 'im_min', im_min;
#endif ! CRF_PRNT_MAP_INFO

        call crf_plot_map( n_grid_x, n_grid_y, re_lo, re_hi, im_lo, im_hi, func )

#endif ! CRF_PLOT_MAP_VAR

        return;
        
    end function find_root


#ifdef CRF_PLOT_MAP
    subroutine crf_plot_map( n_grid_x, n_grid_y, re_lo, re_hi, im_lo, im_hi, func )

        implicit none;
        integer, intent(in):: n_grid_x, n_grid_y;
        double precision, intent(in):: re_lo, re_hi, im_lo, im_hi;
        double precision:: d_re, d_im
        double complex, external:: func;
        double complex:: grid(n_grid_x,n_grid_y), cw
        double precision:: grid_re(n_grid_x), grid_im(n_grid_y), w
        integer:: i, j

        ! create the grid
        d_re = (re_hi - re_lo) / (n_grid_x - 1);
        d_im = (im_hi - im_lo) / (n_grid_y - 1);
        forall(i=1:n_grid_x) grid_re(i) = re_lo + (i-1) * d_re;
        forall(i=1:n_grid_y) grid_im(i) = im_lo + (i-1) * d_im;
        forall(i=1:n_grid_x,j=1:n_grid_y) grid(i,j) = dcmplx( grid_re(i), grid_im(j) );



#ifdef CRF_PRNT_MAP_INFO
        print*, 're_a', minval(grid_re); ! left bounary of the rectangle
        print*, 're_b', maxval(grid_re); ! right bounary of the rectangle
        print*, 'im_a', minval(grid_im); ! bottom bounary of the rectangle
        print*, 'im_b', maxval(grid_im); ! top bounary of the rectangle
#endif

        ! evaluate the function and record the elevations
        open( CRF_PLOT_FILE_ID, file=CRF_PLOT_MAP_FILE_NAME );
        do j = 1, n_grid_y
            do i = 1, n_grid_x
                cw = func(grid(i,j));


#ifdef CRF_NORM_ADD
                w = dabs ( dreal( cw ) ) + dabs ( dimag ( cw ) ); 
#else
                w = cdabs(cw);
#endif
                write(CRF_PLOT_FILE_ID, '(3(e15.8,2x))'), grid(i,j), w;
            end do
        end do
        close(CRF_PLOT_FILE_ID);

        return;

    end subroutine crf_plot_map
#endif !CRF_PLOT_MAP
