# Instructions for using the free-space solver in MPI-AMRVAC
#
# 1. Clone the octree-mg repository somewhere. Below, it is assumed to be in
# ~/git/octree-mg
#
# 2. Type `make` in `the octree-mg/poisson_3d_fft` folder
#
# 3. Copy the file `octree-mg/single_module/m_free_space.f90` to the folder
# with your MPI-AMRVAC user code (with the mod_usr.t in it), and rename it to `m_free_space.t`
#
# 3. Use/Create this file `local.make` with the following content

# Adjust this to where you have installed octree-mg
OCTREE_MG_DIR:=$${HOME}/git/octree-mg

INC_DIRS+=$(OCTREE_MG_DIR)/poisson_3d_fft
LIB_DIRS+=$(OCTREE_MG_DIR)/poisson_3d_fft
LIBS+=pois3dfft

mod_usr.o: m_free_space.o
amrvac: m_free_space.o
