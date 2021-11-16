/*************************************************************************
 *    CompuCell - A software framework for multimodel simulations of     *
 * biocomplexity problems Copyright (C) 2003 University of Notre Dame,   *
 *                             Indiana                                   *
 *                                                                       *
 * This program is free software; IF YOU AGREE TO CITE USE OF CompuCell  *
 *  IN ALL RELATED RESEARCH PUBLICATIONS according to the terms of the   *
 *  CompuCell GNU General Public License RIDER you can redistribute it   *
 * and/or modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of the   *
 *         License, or (at your option) any later version.               *
 *                                                                       *
 * This program is distributed in the hope that it will be useful, but   *
 *      WITHOUT ANY WARRANTY; without even the implied warranty of       *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    *
 *             General Public License for more details.                  *
 *                                                                       *
 *  You should have received a copy of the GNU General Public License    *
 *     along with this program; if not, write to the Free Software       *
 *      Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.        *
 *************************************************************************/



#ifndef MITOSISINITIALISATIONSTEPPABLE_H
#define MITOSISINITIALISATIONSTEPPABLE_H

#include <CompuCell3D/Plugin.h>
#include <CompuCell3D/Potts3D/Cell.h>
#include <CompuCell3D/Steppable.h>
#include <BasicUtils/BasicClassAccessor.h>
#include <CompuCell3D/Field3D/Dim3D.h>
#include <vector>
#include <iostream>
#include <fstream>




#include "MitosisInitialisationDLLSpecifier.h"

namespace CompuCell3D {


  template <class T> class Field3D;
  template <class T> class WatchableField3D;
  
  class Potts3D;
  class PixelTracker;
  class PixelTrackerPlugin;
  class BoundaryPixelTracker;
  class BoundaryPixelTrackerPlugin;
  class CellG;
  class BoundaryStrategy;
  class PixelTrackerData;
  class BoundaryPixelTrackerData;
  


  class MITOSISINITIALISATION_EXPORT MitosisInitialisation : public Steppable {

    WatchableField3D<CellG *> *cellFieldG;
    Simulator * sim;
    Potts3D *potts;
    int VCAFtype;
    int SCCtype;
    int SCCVolume;		
    int VCAFVolume;		
    float SCCLambdaVolume;
    float VCAFLambdaVolume;
    int SCCSurface;		
    int VCAFSurface;		
    float SCCLambdaSurface;
    float VCAFLambdaSurface;
    int InitialTimeToMitosisSCC;	//The maximum time for mitosis, currently a random numbner from uniform dist. (0,1] is multiplied by this number
    int MeanTimeToMitosisCheckSCC; // the mean time it takes a cell to divide. 
     int TimeToVolumeCheck;	//The maximum time the cell will wait before giving up on trying to grow prior to mitosis.
    float VolumeMultiplierForGrowth;	//The scaling factor we increase the targetVolume of the cell
    float DivisionThresholdMultiplier;	//The scaling factor to set limit to decide if cell qualifies for division (the already increased targetVolume*this number sets the threshold)
    string InitialCellOutputFileName;
    ofstream InitialCellOutputFile;
    string TimeToCellDivisionFileName;
    ofstream TimeToCellDivisionFile;
    int TimeToWriteCells;
    string DivisionTimeDistribution;
    int ErlangShapeParameter;
    string SurfacevolumeRatio; 
    
    BasicClassAccessor<PixelTracker> *pixelTrackerAccessorPtr;
    PixelTrackerPlugin * pixelTrackerPlugin;
    
    
  public:
    MitosisInitialisation();

    virtual ~MitosisInitialisation();
    // SimObject interface
    virtual void init(Simulator *simulator, CC3DXMLElement *_xmlData);
    virtual void extraInit(Simulator *simulator);


    virtual void start();
    virtual void step(const unsigned int currentStep);
    //virtual void finish() {}
    virtual void finish();
    //virtual void InitialCellOutput();

    //SteerableObject interface
    virtual void update(CC3DXMLElement *_xmlData, bool _fullInitFlag=false);
    virtual std::string steerableName();
	 virtual std::string toString();
  };
};
#endif
