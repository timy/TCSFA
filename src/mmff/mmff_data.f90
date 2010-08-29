module mod_mmff_mpi_info

    implicit none;

    type T_mmff_mpi_info
        integer:: i_rank;
        integer:: n_proc;
    end type T_mmff_mpi_info

end module mod_mmff_mpi_info



module mod_mmff_core_para

    parameter ( MAX_COUNT = 40 )
    parameter ( MAX_STR_LEN = 20 )
    parameter ( MPI_DOUBLE_PRECISION = 1275070495 )
    parameter ( MPI_DOUBLE_COMPLEX = 1275072546 )
    parameter ( MPI_INTEGER = 1275069467 )


    type T_ptr_arr_dbl
        double precision, pointer:: ptr(:);
    end type T_ptr_arr_dbl

    type T_ptr_arr_dcp
        double complex, pointer:: ptr(:);
    end type T_ptr_arr_dcp

    type T_ptr_arr_int
        integer, pointer:: ptr(:);
    end type T_ptr_arr_int


    type T_arry_info
        integer:: id;                                      ! id of the data array
        integer:: type;                                    ! data type of the array
        integer:: queu;                                    ! the queue where the array is
        character(len=MAX_STR_LEN):: name;                 ! the alias
        type( T_ptr_arr_dbl ):: arry_dbl;                  ! pointer to double
        type( T_ptr_arr_dcp ):: arry_dcp;                  ! pointer to double complex
        type( T_ptr_arr_int ):: arry_int;                  ! pointer to integer
        integer:: b_stat;                                  ! if the array is allocated
    end type T_arry_info


    type T_queu_info
        integer:: id;                                      ! id
        integer:: n_data;                                  ! the array size of the this queue
        integer:: i_arry;                                  ! the id of the first array in this queue
        integer:: b_arry;                                  ! if this queue contains the data array
        integer, allocatable:: a_task(:);                  ! number of workload for each process
    end type T_queu_info


    integer:: arry_count;                                  ! number of all arrays from start on
    integer:: queu_count;                                  ! number of all queues from start on
    integer:: n_resc;                                      ! total number of process resources
    integer:: i_rank;                                      ! current process id
    type(T_arry_info):: arry_info(MAX_COUNT);              ! information of the arrays in queues
    type(T_queu_info):: queu_info(MAX_COUNT);              ! information of the queues

end module mod_mmff_core_para
