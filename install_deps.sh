source setup.sh

mkdir dependencies
cd dependencies
TOP_DIR=$(pwd)

#MPICH
mkdir mpich
cd mpich
MPICH_DIR=$(pwd)
wget https://www.mpich.org/static/downloads/4.0.3/mpich-4.0.3.tar.gz
tar -xvkf mpich-4.0.3.tar.gz
cd mpich-4.0.3/
./configure --prefix=$MPICH_DIR/mpich-install --disable-f08 --disable-collalgo-tests  &> con.log
make &> make&>.log
make install &> install.log
cd $TOP_DIR

#HDF5
mkdir hdf5
cd hdf5
HDF5_DIR=$(pwd)
cp /pnfs/uboone/resilient/users/markross/tars/hdf5-1.12.2.tar.gz .
tar -xvkf hdf5-1.12.2.tar.gz
cd hdf5-1.12.2/
CC=$MPICH_DIR/mpich-install/bin/mpicc ./configure  --enable-parallel --prefix=$HDF5_DIR-install/
make &> make.log
make check &> make.check.log
make install &> install.log
make check-install &> install.check.log
cd $TOP_DIR

#HighFive
mkdir highfive
cd highfive
HIGHFIVE_DIR=$(pwd)
git clone https://github.com/BlueBrain/HighFive.git
cd HighFive/
git checkout -b feature/sbnfit_compatable_v2.4 v2.3
cd $TOP_DIR

#DIY
mkdir di
cd diy
DIY_DIR=$(pwd)
git clone https://github.com/diatomic/diy.git
cd diy
git checkout -b feature/confirmed_sbnfit_compat ef460f828c36e8f533ccec1f551024f4da896165
cd $TOP_DIR

