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

#include <CompuCell3D/Plugin.h>
#include <CompuCell3D/Simulator.h>
#include <CompuCell3D/Potts3D/Potts3D.h>
#include <CompuCell3D/Field3D/Field3D.h>
#include <CompuCell3D/Field3D/WatchableField3D.h>
using namespace CompuCell3D;
#include <iostream>
#include <numeric>
using namespace std;
#include <CompuCell3D/Potts3D/CellInventory.h>
#include "SpheroidCellBoundaryECMPenetration.h"
#include "CompuCell3D/plugins/BoundaryPixelTracker/BoundaryPixelTrackerPlugin.h"
#include "CompuCell3D/plugins/BoundaryPixelTracker/BoundaryPixelTracker.h"
#include "CompuCell3D/plugins/PixelTracker/PixelTrackerPlugin.h"
#include "CompuCell3D/plugins/PixelTracker/PixelTracker.h"

SpheroidCellBoundaryECMPenetration::SpheroidCellBoundaryECMPenetration() :  cellFieldG(0),sim(0),potts(0) {}

SpheroidCellBoundaryECMPenetration::~SpheroidCellBoundaryECMPenetration() {
}


void SpheroidCellBoundaryECMPenetration::init(Simulator *simulator, CC3DXMLElement *_xmlData) {
	potts = simulator->getPotts();
	sim=simulator;
	cellFieldG = (WatchableField3D<CellG *> *)potts->getCellFieldG();
	simulator->registerSteerableObject(this);
	update(_xmlData);
	bool pluginAlreadyRegisteredFlag;
	pixelTrackerPlugin=(PixelTrackerPlugin*)Simulator::pluginManager.get("PixelTracker",&pluginAlreadyRegisteredFlag);
	if(!pluginAlreadyRegisteredFlag)
		pixelTrackerPlugin->init(simulator);
	boundaryPixelTrackerPlugin=(BoundaryPixelTrackerPlugin*)Simulator::pluginManager.get("BoundaryPixelTracker",&pluginAlreadyRegisteredFlag);
	if(!pluginAlreadyRegisteredFlag) 
		boundaryPixelTrackerPlugin->init(simulator);
	pixelTrackerAccessorPtr=pixelTrackerPlugin->getPixelTrackerAccessorPtr();
	boundaryPixelTrackerAccessorPtr=boundaryPixelTrackerPlugin->getBoundaryPixelTrackerAccessorPtr();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SpheroidCellBoundaryECMPenetration::extraInit(Simulator *simulator){

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SpheroidCellBoundaryECMPenetration::start(){
	cout<<"within the start function of ECM"<<endl;
	const int x=dimx;
	const int y=dimy;
	const int z=dimz;
	cellIndexMatrix = new int**[x];
        timePointMatrix	 = new int**[x];
        cellBoundaryMatrix = new int**[x];
        timePointBoundaryMatrix = new int**[x];
	for(int i = 0; i < x; i++){
		cellIndexMatrix[i]= new int*[y];
		timePointMatrix[i]= new int*[y];
		cellBoundaryMatrix[i] = new int*[y];
        timePointBoundaryMatrix[i] = new int*[y];		
		for (int j = 0; j< y; j++){
			cellIndexMatrix[i][j] = new int[z];
            timePointMatrix[i][j] = new int[z];
            cellBoundaryMatrix[i][j] = new int[z];
        	timePointBoundaryMatrix[i][j] = new int[z];
		}
	}
	for (int i = 0; i< x; i++){
	for (int j = 0; j< y; j++){
	for (int k = 0; k< z; k++){
		cellIndexMatrix[i][j][k]=-1;
		timePointMatrix[i][j][k]=0;
		cellBoundaryMatrix[i][j][k] = -1;
    	timePointBoundaryMatrix[i][j][k] = 0;
	}}}
	
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SpheroidCellBoundaryECMPenetration::step(const unsigned int currentStep){
	CellInventory *cellInventoryPtr=& potts->getCellInventory();
	CellInventory::cellInventoryIterator cInvItr;
	CellG *cell;
	int neighbourOrder1=max(neighbourOrder-1,0);
	for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
		cell=cellInventoryPtr->getCell(cInvItr);
        cell->canremodel_x.clear();
        cell->canremodel_y.clear();
        cell->canremodel_z.clear();
		set<PixelTrackerData> cellPixels=pixelTrackerAccessorPtr->get(cell->extraAttribPtr)->pixelSet;
		for(set<PixelTrackerData>::iterator sitr=cellPixels.begin() ; sitr != cellPixels.end() ;++sitr){
			int xx=sitr->pixel.x;
			int yy=sitr->pixel.y;
			int zz=sitr->pixel.z;
			cellIndexMatrix[xx][yy][zz]=cell->id;
			timePointMatrix[xx][yy][zz]=currentStep;
		}
        }
	
	for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
        	cell=cellInventoryPtr->getCell(cInvItr);
        	set<BoundaryPixelTrackerData> cellBoundaryPixels=boundaryPixelTrackerAccessorPtr->get(cell->extraAttribPtr)->pixelSet;
		for(set<BoundaryPixelTrackerData>::iterator sitr=cellBoundaryPixels.begin() ; sitr != cellBoundaryPixels.end() ;++sitr){
			int xx=sitr->pixel.x;
			int yy=sitr->pixel.y;
			int zz=sitr->pixel.z;
			for(int i=-neighbourOrder1;i<neighbourOrder1+1;i++){
				for(int j=-neighbourOrder1;j<neighbourOrder1+1;j++){
					for(int k=-neighbourOrder1;k<neighbourOrder1+1;k++){
						int currz=zz+k;
						int curry=yy+j;
						int currx=xx+i;
						if(currz>0 && currz<dimz && curry>0 && curry<dimy && currx>0 && currx<dimx){
							if(timePointBoundaryMatrix[currx][curry][currz]!=currentStep || cellBoundaryMatrix[currx][curry][currz]!=cell->id ){
								if(timePointMatrix[currx][curry][currz]==currentStep && cellIndexMatrix[currx][curry][currz]==cell->id ){
									cell->canremodel_x.push_back(currx);
									cell->canremodel_y.push_back(curry);
									cell->canremodel_z.push_back(currz);
									cellBoundaryMatrix[currx][curry][currz]=cell->id;
									timePointBoundaryMatrix[currx][curry][currz]=currentStep;
								}
							}
							
						}
					}
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SpheroidCellBoundaryECMPenetration::finish(){	
}

void SpheroidCellBoundaryECMPenetration::update(CC3DXMLElement *_xmlData, bool _fullInitFlag){
     if(_xmlData->findElement("SCCtype")){SCCtype=_xmlData->getFirstElement("SCCtype")->getDouble();}
     else{SCCtype=1;cerr<<"SCCtype not specified, defaul t value set to: "<<SCCtype<<endl;}
     if(_xmlData->findElement("VCAFtype")){VCAFtype=_xmlData->getFirstElement("VCAFtype")->getDouble();}
     else{VCAFtype=2;cerr<<"VCAFtype not specified, defaul t value set to: "<<VCAFtype<<endl;}
     if(_xmlData->findElement("NeighbourOrder")){neighbourOrder=_xmlData->getFirstElement("NeighbourOrder")->getDouble();}
     else{neighbourOrder=1;cerr<<"NeighbourOrder not specified, default value set to: "<<neighbourOrder<<endl;}
     if(_xmlData->findElement("dimx")){dimx=_xmlData->getFirstElement("dimx")->getDouble();}
     if(_xmlData->findElement("dimy")){dimy=_xmlData->getFirstElement("dimy")->getDouble();}
     if(_xmlData->findElement("dimz")){dimz=_xmlData->getFirstElement("dimz")->getDouble();}
     cout<<"READ stuff"<<endl;      
}

std::string SpheroidCellBoundaryECMPenetration::toString(){
   return "SpheroidCellBoundaryECMPenetration";
}


std::string SpheroidCellBoundaryECMPenetration::steerableName(){
   return toString();
}


