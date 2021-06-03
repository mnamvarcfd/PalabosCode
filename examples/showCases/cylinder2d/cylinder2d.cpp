/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 * 
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * The most recent release of Palabos can be downloaded at 
 * <https://palabos.unige.ch/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file
 * Flow around a 2D cylinder inside a channel, with the creation of a von
 * Karman vortex street. This example makes use of bounce-back nodes to
 * describe the shape of the cylinder. The outlet is modeled through a
 * Neumann (zero velocity-gradient) condition.
 */

#include "palabos2D.h"
#include "palabos2D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor

/// Velocity on the parabolic Poiseuille profile
T poiseuilleVelocity(plint iY, IncomprFlowParam<T> const& parameters) {
    T y = (T)iY / parameters.getResolution();
    return 4.*parameters.getLatticeU() * (y-y*y);
}

/// Linearly decreasing pressure profile
T poiseuillePressure(plint iX, IncomprFlowParam<T> const& parameters) {
    T Lx = parameters.getNx()-1;
    T Ly = parameters.getNy()-1;
    return 8.*parameters.getLatticeNu()*parameters.getLatticeU() / (Ly*Ly) * (Lx/(T)2-(T)iX);
}

/// Convert pressure to density according to ideal gas law
T poiseuilleDensity(plint iX, IncomprFlowParam<T> const& parameters) {
    return poiseuillePressure(iX,parameters)*DESCRIPTOR<T>::invCs2 + (T)1;
}

/// A functional, used to initialize the velocity for the boundary conditions
template<typename T>
class PoiseuilleVelocity {
public:
    PoiseuilleVelocity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, Array<T,2>& u) const {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }
private:
    IncomprFlowParam<T> parameters;
};

/// A functional, used to initialize a pressure boundary to constant density
template<typename T>
class ConstantDensity {
public:
    ConstantDensity(T density_)
        : density(density_)
    { }
    T operator()(plint iX, plint iY) const {
        return density;
    }
private:
    T density;
};

/// A functional, used to create an initial condition for the density and velocity
template<typename T>
class PoiseuilleVelocityAndDensity {
public:
    PoiseuilleVelocityAndDensity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, T& rho, Array<T,2>& u) const {
        rho = poiseuilleDensity(iX,parameters);
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }
private:
    IncomprFlowParam<T> parameters;
};

/// A functional, used to instantiate bounce-back nodes at the locations of the cylinder
template<typename T>
class CylinderShapeDomain2D : public plb::DomainFunctional2D {
public:
    CylinderShapeDomain2D(plb::plint cx_, plb::plint cy_, plb::plint radius)
        : cx(cx_),
          cy(cy_),
          radiusSqr(plb::util::sqr(radius))
    { }
    virtual bool operator() (plb::plint iX, plb::plint iY) const {
        return plb::util::sqr(iX-cx) + plb::util::sqr(iY-cy) <= radiusSqr;
    }
    virtual CylinderShapeDomain2D<T>* clone() const {
        return new CylinderShapeDomain2D<T>(*this);
    }
private:
    plb::plint cx;
    plb::plint cy;
    plb::plint radiusSqr;
};


void cylinderSetup( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                    IncomprFlowParam<T> const& parameters,
                    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    Box2D outlet(nx-1,nx-1, 1, ny-2);

    // Create Velocity boundary conditions everywhere
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, 0, 1, ny-2) );
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, nx-1, 0, 0) );
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, nx-1, ny-1, ny-1) );
    // .. except on right boundary, where we prefer an outflow condition
    //    (zero velocity-gradient).
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(nx-1, nx-1, 1, ny-2), boundary::outflow );

    setBoundaryVelocity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocity<T>(parameters) );
    setBoundaryDensity (
            lattice, outlet,
            ConstantDensity<T>(1.) );
    initializeAtEquilibrium (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocityAndDensity<T>(parameters) );


    double dx = 1.0;
    double Xc = ny;
    double Yc = ny;

    double R = 1.285 * ny;

    double offset = .1 * ny;

    int nCircles = 19;
    double* xc = new double[nCircles];
    double* yc = new double[nCircles];
    double* radius = new double[nCircles];

    xc[0] = offset + Xc * 0.74; yc[0] = Yc * 0.26; radius[0] = R * 0.06;
    xc[1] = offset + Xc * 0.26; yc[1] = Yc * 0.24; radius[1] = R * 0.21;
    xc[2] = offset + Xc * 0.15; yc[2] = Yc * 0.54; radius[2] = R * 0.09;
    xc[3] = offset + Xc * 0.23; yc[3] = Yc * 0.81; radius[3] = R * 0.17;
    xc[4] = offset + Xc * 0.60; yc[4] = Yc * 0.13; radius[4] = R * 0.11;
    xc[5] = offset + Xc * 0.62; yc[5] = Yc * 0.58; radius[5] = R * 0.26;
    xc[6] = offset + Xc * 0.48; yc[6] = Yc * 0.91; radius[6] = R * 0.08;
    xc[7] = offset + Xc * 0.81; yc[7] = Yc * 0.10; radius[7] = R * 0.09;
    xc[8] = offset + Xc * 1.10; yc[8] = Yc * 0.34; radius[8] = R * 0.26;
    xc[9] = offset + Xc * 0.75; yc[9] = Yc * 0.91; radius[9] = R * 0.08;
    xc[10] = offset + Xc * 1.03; yc[10] = Yc * 0.8; radius[10] = R * 0.19;
    xc[11] = offset + Xc * 1.36; yc[11] = Yc * 0.09; radius[11] = R * 0.07;
    xc[12] = offset + Xc * 1.39; yc[12] = Yc * 0.68; radius[12] = R * 0.16;
    xc[13] = offset + Xc * 1.29; yc[13] = Yc * 0.91; radius[13] = R * 0.07;
    xc[14] = offset + Xc * 1.66; yc[14] = Yc * 0.31; radius[14] = R * 0.28;
    xc[15] = offset + Xc * 1.68; yc[15] = Yc * 0.72; radius[15] = R * 0.1;
    xc[16] = offset + Xc * 1.56; yc[16] = Yc * 0.89; radius[16] = R * 0.09;
    xc[17] = offset + Xc * 1.86; yc[17] = Yc * 0.63; radius[17] = R * 0.08;
    xc[18] = offset + Xc * 1.83; yc[18] = Yc * 0.88; radius[18] = R * 0.10;


    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[0], yc[0], radius[0]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[1], yc[1], radius[1]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[2], yc[2], radius[2]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[3], yc[3], radius[3]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[4], yc[4], radius[4]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[5], yc[5], radius[5]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[6], yc[6], radius[6]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[7], yc[7], radius[7]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[8], yc[8], radius[8]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[9], yc[9], radius[9]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[10], yc[10], radius[10]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[11], yc[11], radius[11]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[12], yc[12], radius[12]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[13], yc[13], radius[13]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[14], yc[14], radius[14]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[15], yc[15], radius[15]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[16], yc[16], radius[16]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[17], yc[17], radius[17]), new plb::BounceBack<T, DESCRIPTOR>);
    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(xc[18], yc[18], radius[18]), new plb::BounceBack<T, DESCRIPTOR>);

    lattice.initialize();
}

void writeGif(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, plint iter)
{
    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("u", iter, 6),
                               *computeVelocityNorm(lattice) );
}

void writeVTK(MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<2,float>(*computeVelocity(lattice), "velocity", dx/dt);
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");

    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // uMax
            (T) 600.,  // Re
            1000,       // N
            2.,        // lx
            1.         // ly 
    );
    const T logT     = (T)0.02;
    const T imSave   = (T)0.06;
    const T vtkSave  = (T)1.;
    const T maxT     = (T)20.1;

    writeLogFile(parameters, "Poiseuille flow");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
            parameters.getNx(), parameters.getNy(),
            new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    cylinderSetup(lattice, parameters, *boundaryCondition);

    // Main loop over time iterations.
    for (plint iT=0; iT*parameters.getDeltaT()<maxT; ++iT) {
        // At this point, the state of the lattice corresponds to the
        //   discrete time iT. However, the stored averages (getStoredAverageEnergy
        //   and getStoredAverageDensity) correspond to the previous time iT-1.

       if (iT%parameters.nStep(imSave)==0) {
            pcout << "Saving Gif ..." << endl;
            writeGif(lattice, iT);
        }

        if (iT%parameters.nStep(vtkSave)==0 && iT>0) {
            pcout << "Saving VTK file ..." << endl;
            writeVTK(lattice, parameters, iT);
        }

        if (iT%parameters.nStep(logT)==0) {
            pcout << "step " << iT
                  << "; t=" << iT*parameters.getDeltaT();
        }

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();

        // At this point, the state of the lattice corresponds to the
        //   discrete time iT+1, and the stored averages are upgraded to time iT.
        if (iT%parameters.nStep(logT)==0) {
            pcout << "; av energy ="
                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << "; av rho ="
                  << getStoredAverageDensity<T>(lattice) << endl;
        }
    }
    
    delete boundaryCondition;
}
