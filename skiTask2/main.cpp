#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "mpi.h"

#define TARGET log(2.0)/2.0 - 5.0/16.0

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Request request_data, request_ind;
    MPI_Status status_data, status_ind;
    int index = 1, block_size = 50000/(size-1), iter = 0;
    double x, y, z, *xyzBuff, *recv_buff;
    double int_sum = 0.0, master_sum = 0.0, result = 0.0;
    double eps = atof(argv[1]);
    double delta = 1.0/(size-1);
    if (!rank) {
        srand(10);
        xyzBuff = new double[(size - 1) * block_size * 3];
        for (int rank_slave = 1; rank_slave < size; ++rank_slave)
            MPI_Isend(&index, 1, MPI_INT, rank_slave, 2, MPI_COMM_WORLD, &request_ind);
    } else {
        recv_buff = new double[block_size * 3];
    }
    double start_time = MPI_Wtime();
    while (index) {
        if (!rank) //Master
        {
            for (int rank_slave = 1; rank_slave < size; ++rank_slave) {
                // generating points in interval (0,1)
                for (int rand_point_i = 0; rand_point_i < block_size; ++rand_point_i) {
                    x = delta*(rank_slave-1) + delta*(double) rand() / (RAND_MAX);
                    y = (double) rand() / (RAND_MAX);
                    z = (double) rand() / (RAND_MAX);
                    xyzBuff[block_size * (rank_slave - 1) * 3 + 3 * rand_point_i] = x;
                    xyzBuff[block_size * (rank_slave - 1) * 3 + 3 * rand_point_i + 1] = y;
                    xyzBuff[block_size * (rank_slave - 1) * 3 + 3 * rand_point_i + 2] = z;
                }
                MPI_Isend(xyzBuff + (rank_slave - 1) * block_size * 3, 3 * block_size, MPI_DOUBLE,
                          rank_slave, 1, MPI_COMM_WORLD, &request_data);
            }
        }
        if (rank) // Slave
        {
            MPI_Recv(&index, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status_ind);
            if (!index) {
                break;
            }
            MPI_Recv(recv_buff, block_size * 3, MPI_DOUBLE, 0, 1,
                     MPI_COMM_WORLD, &status_data);
            for (int rand_point_i = 0; rand_point_i < block_size; ++rand_point_i) {
                double sum_xyz =
                        recv_buff[3 * rand_point_i] + recv_buff[3 * rand_point_i + 1] + recv_buff[3 * rand_point_i + 2];
                if (sum_xyz <= 1) {
                    sum_xyz += 1;
                    int_sum += 1.0 / (sum_xyz * sum_xyz * sum_xyz);
                }
            }
        }
        MPI_Reduce(&int_sum, &master_sum, 1, MPI_DOUBLE,
                   MPI_SUM, 0, MPI_COMM_WORLD);
        if (!rank) {
            ++iter;
            result = master_sum / (iter * block_size * (size - 1));
            if (std::fabs(TARGET - result) < eps) {
                index = 0;
            }
            for (int i = 1; i < size; ++i) {
                MPI_Isend(&index, 1, MPI_INT, i, 2, MPI_COMM_WORLD,
                          &request_ind);
            }
        }
    }
    double end_time = MPI_Wtime();
    double result_time, time = end_time - start_time;
    MPI_Reduce(&time, &result_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "True: " << TARGET << std::endl;
        std::cout << "Integral: " << result << std::endl;
        std::cout << "Eps: " << std::abs(TARGET - result) << std::endl;
        std::cout << "N points: " << (size - 1) * block_size * iter << std::endl;
        std::cout << "Time: " << result_time << std::endl;
    }
    MPI_Finalize();
    return 0;
}
