/* Copyright (c) 2019 CSC Training */
/* Copyright (c) 2021 ENCCS */
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void test_on_gpu(void) {
    int on_device = 0;
    #pragma omp target teams map(from:on_device)
    {
        #pragma omp parallel
        {
            #pragma omp master
            {
                if (0 == omp_get_team_num()) {
                    on_device = !omp_is_initial_device();
                }
            }
        }
    }
    printf("on GPU: %s\n", on_device ? "yes" : "no");
}

int main() 
{
  test_on_gpu();
}
