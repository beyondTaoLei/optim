#include "fd.h"
double time_count(double time1, double *time2){

    double dt_time;
    *time2=MPI_Wtime();
    dt_time=*time2-time1;
    return dt_time;

}
