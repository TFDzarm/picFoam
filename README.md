# picFoam
**picFoam** is a fully kinetic electrostatic Particle-in-Cell(PIC) solver, including Monte Carlo Collisions(MCC), for non-equilibrium plasma research.

The solver is developed at ZARM (Center of Applied Space Technology and Microgravity, University of Bremen).

## Publications

1. [picFoam: An OpenFOAM based electrostatic Particle-in-Cell solver](https://doi.org/10.1016/j.cpc.2021.107853) ([preprint](http://arxiv.org/abs/2012.14724))

# Instructions

picFoam requires [OpenFOAM version 8](https://github.com/OpenFOAM/OpenFOAM-8).

Detailed installation instructions for OpenFOAM can be found [here](https://openfoam.org/download/source/).

## OpenFOAM installation

Summery of the OpenFOAM installtion (link above):

1. Compiler: GCC version > 4.8 or above; or LLVM Clang version > 3.6 or above; or Intel ICC version > 17.0.4 or above
2. Package dependencies (for Ubuntu version 18.04 or above):

       OpenFOAM:

		build-essential
		flex
		bison
		git-core
		cmake
		zlib1g-dev
		libboost-system-dev
		libboost-thread-dev
		libopenmpi-dev
		openmpi-bin
		gnuplot
		libreadline-dev
		libncurses-dev
		libxt-dev

	ParaView:

		libqt5x11extras5-dev
		libxt-dev
		qt5-default
		qttools5-dev curl

3. Clone OpenFOAM repositories (recommended installation directory $HOME/OpenFOAM):
    git clone https://github.com/OpenFOAM/OpenFOAM-8.git
    git clone https://github.com/OpenFOAM/ThirdParty-8.git

4. Set the environment variables:
    source $HOME/OpenFOAM/OpenFOAM-8/etc/bashrc

5. Compile Scotch/PT-Scotch: Run the **Allwmake** script in $HOME/OpenFOAM/ThirdParty-8
6. Compile OpenFOAM: Run the **Allwmake** script in $HOME/OpenFOAM/OpenFOAM-8
7. Compile ParaView: Run the **makeParaView** script in $HOME/OpenFOAM/ThirdParty-8

## picFoam installation

Run the **Allwmake** script inside the picFoam directory.

---
## picFoam tutorials

To run a tutorial execute the respective **run.sh** script.
By default meshes are constructed with *blockMesh*, but [Gmsh](https://gmsh.info/) geometry files are also included.

If you wish to use the Gmsh geometry, installation of Gmsh version [3.0.6](https://gmsh.info/bin/Linux/gmsh-3.0.6-Linux64.tgz) is recommended.
For newer versions run Gmsh with **-format msh22** to output an older mesh format, which OpenFOAM can handle.

---
## picFoam extension

If you wish to extend picFoam with new submodels see the file EXTENSION_INSTRUCTIONS.md.
