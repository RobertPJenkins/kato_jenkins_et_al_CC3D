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

#include <stdlib.h>
#include <CompuCell3D/Plugin.h>
#include <CompuCell3D/Simulator.h>
#include <CompuCell3D/Potts3D/Potts3D.h>
#include <CompuCell3D/Field3D/Field3D.h>
#include <CompuCell3D/Field3D/WatchableField3D.h>
using namespace CompuCell3D;


#include <CompuCell3D/Potts3D/CellInventory.h>

#include <iostream>
#include <numeric>
using namespace std;


#include "SpheroidCentroid.h"

#include "CompuCell3D/plugins/PixelTracker/PixelTrackerPlugin.h"
#include "CompuCell3D/plugins/PixelTracker/PixelTracker.h"

SpheroidCentroid::SpheroidCentroid() :  cellFieldG(0),sim(0),potts(0) {}

SpheroidCentroid::~SpheroidCentroid() {
}


void SpheroidCentroid::init(Simulator *simulator, CC3DXMLElement *_xmlData) {
  potts = simulator->getPotts();
  sim=simulator;
  cellFieldG = (WatchableField3D<CellG *> *)potts->getCellFieldG();
  
  simulator->registerSteerableObject(this);
  update(_xmlData);
  bool pluginAlreadyRegisteredFlag;
  
  pixelTrackerPlugin=(PixelTrackerPlugin*)Simulator::pluginManager.get("PixelTracker",&pluginAlreadyRegisteredFlag); //this will load VolumeTracker plugin if it is not already loaded
	if(!pluginAlreadyRegisteredFlag)
		pixelTrackerPlugin->init(simulator);

  pixelTrackerAccessorPtr=pixelTrackerPlugin->getPixelTrackerAccessorPtr();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SpheroidCentroid::extraInit(Simulator *simulator){

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SpheroidCentroid::start(){

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SpheroidCentroid::step(const unsigned int currentStep){
	
	CellInventory *cellInventoryPtr=& potts->getCellInventory();
	CellInventory::cellInventoryIterator cInvItr;
	CellG *cell;
	
	if(currentStep==0){
		for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
			cell=cellInventoryPtr->getCell(cInvItr);
			cell->SpheroidCentroidX=SpheroidCentroidVar[0];
			cell->SpheroidCentroidY=SpheroidCentroidVar[1];
			cell->SpheroidCentroidZ=SpheroidCentroidVar[2];
		}
	}
	
	vector <int> xcoord,ycoord,zcoord;
	for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
		cell=cellInventoryPtr->getCell(cInvItr);
		if((int)cell->type==SCCtype){
			 set<PixelTrackerData> cellPixels=pixelTrackerAccessorPtr->get(cell->extraAttribPtr)->pixelSet;
			for(set<PixelTrackerData>::iterator sitr=cellPixels.begin() ; sitr != cellPixels.end() ;++sitr){
				int xx=sitr->pixel.x;int yy=sitr->pixel.y;int zz=sitr->pixel.z;
				xcoord.push_back(xx);ycoord.push_back(yy);zcoord.push_back(zz);
			}
		}
    }
	float centroidx=float(std::accumulate(xcoord.begin(),xcoord.end(),0.0))/float(xcoord.size());
	float centroidy=float(std::accumulate(ycoord.begin(),ycoord.end(),0.0))/float(ycoord.size());
	float centroidz=float(std::accumulate(zcoord.begin(),zcoord.end(),0.0))/float(zcoord.size());
	
	for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
		cell=cellInventoryPtr->getCell(cInvItr);
		cell->SpheroidCentroidX=int(round(centroidx));
		cell->SpheroidCentroidY=int(round(centroidy));
		cell->SpheroidCentroidZ=int(round(centroidz));
	}
}






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SpheroidCentroid::update(CC3DXMLElement *_xmlData, bool _fullInitFlag){
     if(_xmlData->findElement("SCCtype")){SCCtype=_xmlData->getFirstElement("SCCtype")->getDouble();}
     else{SCCtype=1;cerr<<"SCCtype not specified, default value set to: "<<SCCtype<<endl;}
     
     if(_xmlData->findElement("VCAFtype")){VCAFtype=_xmlData->getFirstElement("VCAFtype")->getDouble();}
     else{VCAFtype=2;cerr<<"VCAFtype not specified, default value set to: "<<VCAFtype<<endl;}
     if(_xmlData->findElement("SpheroidCentroid_x"))
        SpheroidCentroidVar[0]=_xmlData->getFirstElement("SpheroidCentroid_x")->getDouble();
     if(_xmlData->findElement("SpheroidCentroid_y"))
        SpheroidCentroidVar[1]=_xmlData->getFirstElement("SpheroidCentroid_y")->getDouble();
     if(_xmlData->findElement("SpheroidCentroid_z"))
        SpheroidCentroidVar[2]=_xmlData->getFirstElement("SpheroidCentroid_z")->getDouble();   
}

std::string SpheroidCentroid::toString(){
   return "SpheroidCentroid";
}


std::string SpheroidCentroid::steerableName(){
   return toString();
}


