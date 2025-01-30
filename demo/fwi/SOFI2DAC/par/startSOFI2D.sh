#----excecute with LAMMPI
#lamboot -v mpihosts
#lamboot -v

mpirun -np 4 nice -19 ../bin/sofi2D ./in_and_out/sofi2D.json | tee ./in_and_out/sofi2D.jout

#----execute with OPENMPI2
#mpirun  --hostfile mpihosts -np 4 nice -19  ../bin/sofi2D ./in_and_out/sofi2D.json | tee ./in_and_out/sofi2D.jout


#merge snapshots
#../bin/snapmerge ./in_and_out/sofi2D.json

