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
 * Groups all the generic 3D template files in the boundaryConditions directory.
 */
#include "boundaryCondition/boundaryDynamics.hh"
#include "boundaryCondition/regularizedBoundaryDynamics.hh"
#include "boundaryCondition/regularizedBoundaryDynamics3D.hh"
#include "boundaryCondition/equilibriumBoundaryDynamics.hh"
#include "boundaryCondition/inamuroAnalyticalDynamics.hh"
#include "boundaryCondition/zouHeBoundary3D.hh"
#include "boundaryCondition/zouHeDynamics.hh"
#include "boundaryCondition/boundaryCondition3D.hh"
#include "boundaryCondition/finiteDifferenceBoundaryProcessor3D.hh"
#include "boundaryCondition/wrappedLocalBoundaryProcessor3D.hh"
#include "boundaryCondition/neumannCondition3D.hh"
#include "boundaryCondition/bounceBackModels.hh"
#include "boundaryCondition/bounceBackModels3D.hh"
#include "boundaryCondition/NLD_boundaryDynamics3D.hh"
#include "boundaryCondition/NLD_boundaries3D.hh"
#include "boundaryCondition/spongeZones3D.hh"
#include "boundaryCondition/partialBBdynamics.hh"

#ifndef PLB_BGP
#ifdef PLB_USE_EIGEN
#include "boundaryCondition/generalizedBoundaryDynamics.hh"
#include "boundaryCondition/generalizedBoundaryCondition3D.hh"
#endif
#endif

