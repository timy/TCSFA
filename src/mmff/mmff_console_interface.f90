
subroutine mmff_init( mpi_info )
    use mod_mmff_mpi_info, only: T_mmff_mpi_info
    type(T_mmff_mpi_info):: mpi_info;
end subroutine mmff_init



subroutine mmff_reset( n_data, flag, h_queu )
    integer, intent(in)::  n_data, flag, h_queu
end subroutine mmff_reset



subroutine mmff_get_task( h_queu, a_task )
    integer, intent(in):: h_queu
    integer, pointer:: a_task(:)
end subroutine mmff_get_task



subroutine mmff_set_task( h_queu, a_task )
    integer, intent(in):: h_queu
    integer, pointer:: a_task(:)
end subroutine mmff_set_task


subroutine mmff_broadcast_master_int( na, a )
    integer, intent(in):: na, a(na)
end subroutine mmff_broadcast_master_int


subroutine mmff_broadcast_slave_int( a )
    integer, pointer:: a(:)
end subroutine mmff_broadcast_slave_int









subroutine mmff_send_data( h_queu )
    integer, intent(in):: h_queu;
end subroutine mmff_send_data



subroutine mmff_recv_data( h_queu )
    integer, intent(in):: h_queu;
end subroutine mmff_recv_data




subroutine mmff_create_data_dbl( a, f, str, h_arry )
    double precision, pointer:: a(:)
    external:: f;
    character(len=*):: str;
    integer, intent(in):: h_arry;
end subroutine mmff_create_data_dbl



subroutine mmff_delete_data_dbl( a, h_arry )
    double precision, pointer:: a(:)
    integer, intent(in):: h_arry
end subroutine mmff_delete_data_dbl




subroutine mmff_create_data_dcp( a, f, str, h_arry )
    double complex, pointer:: a(:)
    external:: f;
    character(len=*):: str;
    integer, intent(in):: h_arry;
end subroutine mmff_create_data_dcp



subroutine mmff_delete_data_dcp( a, h_arry )
    double complex, pointer:: a(:)
    integer, intent(in):: h_arry
end subroutine mmff_delete_data_dcp



subroutine mmff_create_data_int( a, f, str, h_arry )
    integer, pointer:: a(:)
    external:: f;
    character(len=*):: str;
    integer, intent(in):: h_arry;
end subroutine mmff_create_data_int



subroutine mmff_delete_data_int( a, h_arry )
    integer, pointer:: a(:)
    integer, intent(in):: h_arry
end subroutine mmff_delete_data_int
