for dim in 64 128 256 512 1024
do
    echo ======= ${dim}x${dim}
    for p in 1 2 4
    do
        for t in 1 2 4 8
        do
            printf "%3d nodes, %3d ranks,  %3d threads : "  $p $p $t
            #OMP_NUM_THREADS=$t srun -n$p -N$p -c10 --hint=nomultithread ./main $dim $dim 100 0.01 | grep took;
            OMP_NUM_THREADS=$t mpirun -n $p -N 1 ./main $dim $dim 100 0.01 | grep took;
        done
    done
done
