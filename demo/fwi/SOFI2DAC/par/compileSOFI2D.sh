#---- general compilation

cd ../src
#make all
make clean
make sofi2D
cd ../par

#---- compilation only used for speed-up or scale-up benchmark

#cd ../src
#cp hh_elastic.c benchmod_el.c
#cp hh_visco.c benchmod.c
#make sofi2D_bench
#cd ../par
