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
#include "OrganotypicAdhesionECMConcentration.h"
#include "CompuCell3D/plugins/PixelTracker/PixelTrackerPlugin.h"
#include "CompuCell3D/plugins/PixelTracker/PixelTracker.h"
#include "CompuCell3D/plugins/BoundaryPixelTracker/BoundaryPixelTrackerPlugin.h"
#include "CompuCell3D/plugins/BoundaryPixelTracker/BoundaryPixelTracker.h"

OrganotypicAdhesionECMConcentration::OrganotypicAdhesionECMConcentration() :  cellFieldG(0),sim(0),potts(0) {}

OrganotypicAdhesionECMConcentration::~OrganotypicAdhesionECMConcentration() {
}


void OrganotypicAdhesionECMConcentration::init(Simulator *simulator, CC3DXMLElement *_xmlData) {
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
void OrganotypicAdhesionECMConcentration::extraInit(Simulator *simulator){

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void OrganotypicAdhesionECMConcentration::start(){
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

void OrganotypicAdhesionECMConcentration::step(const unsigned int currentStep){
	CellInventory *cellInventoryPtr=& potts->getCellInventory();
	CellInventory::cellInventoryIterator cInvItr;
	CellG *cell;
	if(currentStep==0){
		const char *inputname;
		inputname=CellAdhesionInputFileName.c_str();
		CellAdhesionInputFile.open(inputname,ifstream::in);
		if(!CellAdhesionInputFile.is_open()){cerr<<"could not open "<<CellAdhesionInputFile<<" !!"<<endl;}
		int counter=-1;
		int tempinput1;
		string tempinput2;
		float tempinput3;
		CellAdhesionInputFile>>tempinput2;
		CellAdhesionInputFile>>tempinput2;
		CellAdhesionInputFile>>tempinput2;
		CellAdhesionInputFile>>tempinput2;
		CellAdhesionInputFile>>tempinput2;
		CellAdhesionInputFile>>tempinput2;
		CellAdhesionInputFile>>tempinput2;
		CellAdhesionInputFile>>tempinput2;
		CellAdhesionInputFile>>tempinput2;
		CellAdhesionInputFile>>tempinput2;
		
		for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
			cell=cellInventoryPtr->getCell(cInvItr);
			if (cell->type == SCCtype){
				cell->MaximumVCAFECMAdhesion=0.0;
				cell->CurrentVCAFECMAdhesion=0.0;
				cell->ZeroVCAFECMAdhesion=0.0;
				CellAdhesionInputFile>>tempinput1;
				CellAdhesionInputFile>>tempinput2;
				CellAdhesionInputFile>>tempinput3;
				CellAdhesionInputFile>>tempinput3;
				CellAdhesionInputFile>>cell->MaximumSCCECMAdhesion;
				CellAdhesionInputFile>>tempinput3;
				CellAdhesionInputFile>>tempinput3;
				CellAdhesionInputFile>>tempinput3;
				CellAdhesionInputFile>>tempinput3;
				CellAdhesionInputFile>>tempinput3;
			}
			else if(cell->type == VCAFtype){
				cell->MaximumSCCECMAdhesion=0.0;
				cell->CurrentSCCECMAdhesion=0.0;
				cell->ZeroSCCECMAdhesion=0.0;
				CellAdhesionInputFile>>tempinput1;
				CellAdhesionInputFile>>tempinput2;
				CellAdhesionInputFile>>tempinput3;
				CellAdhesionInputFile>>tempinput3;
				CellAdhesionInputFile>>tempinput3;
				CellAdhesionInputFile>>tempinput3;
				CellAdhesionInputFile>>tempinput3;
				CellAdhesionInputFile>>cell->MaximumVCAFECMAdhesion;
				CellAdhesionInputFile>>tempinput3;
				CellAdhesionInputFile>>tempinput3;
			}
		}
		CellAdhesionInputFile.close();
	}
	for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
		cell=cellInventoryPtr->getCell(cInvItr);
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
			if(sitr==cellBoundaryPixels.begin()){
				int zz=sitr->pixel.z;
				cell->MinZExtent=zz;
			}	
			else{
				if(sitr->pixel.z<cell->MinZExtent){
					cell->MinZExtent=sitr->pixel.z;
				}
			}
				
		}
	}
	
	
	for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
        	cell=cellInventoryPtr->getCell(cInvItr);
        	cell->ECMNeighbourConcentrations.clear();
        	set<BoundaryPixelTrackerData> cellBoundaryPixels=boundaryPixelTrackerAccessorPtr->get(cell->extraAttribPtr)->pixelSet;
		for(set<BoundaryPixelTrackerData>::iterator sitr=cellBoundaryPixels.begin() ; sitr != cellBoundaryPixels.end() ;++sitr){
			int xx=sitr->pixel.x;
			int yy=sitr->pixel.y;
			int zz=sitr->pixel.z;
			for(int i=-1;i<2;i++){
				for(int j=-1;j<2;j++){
					for(int k=-1;k<2;k++){
						int currz=zz+k;
						if(currz>=0 && currz<dimz){
							int currx=xx+i;int curry=yy+j;
							if(currx>=dimx){
								currx=currx-dimx;
							}
							else if(currx<0){
								currx=currx+dimx;
							}
							if(curry>=dimy){
								curry=curry-dimy;
							}
							else if(curry<0){
								curry=curry+dimy;
							}
							if(timePointBoundaryMatrix[currx][curry][currz]!=currentStep || cellBoundaryMatrix[currx][curry][currz]!=cell->id ){
								if(timePointMatrix[currx][curry][currz]!=currentStep || cellIndexMatrix[currx][curry][currz]!=cell->id ){
									cell->ECMNeighbourConcentrations.push_back(cell->ECMfieldbyCell[currx][curry][currz]);
									cellBoundaryMatrix[currx][curry][currz]=cell->id;
									timePointBoundaryMatrix[currx][curry][currz]=currentStep;
								}
							}
							
						}
					}
				}
			}
		}
		cell->MeanSurroundingECMConc=std::accumulate(cell->ECMNeighbourConcentrations.begin(),cell->ECMNeighbourConcentrations.end(),0.0);
		cell->MeanSurroundingECMConc=cell->MeanSurroundingECMConc/cell->ECMNeighbourConcentrations.size();
	}
	
	for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
		cell=cellInventoryPtr->getCell(cInvItr);
		if(cell->MinZExtent>SpaceBoundary){
			if (cell->type == SCCtype){
				cell->CurrentSCCECMAdhesion=SCCSpaceAdhesion;
			}
			else if (cell->type == VCAFtype){
				cell->CurrentVCAFECMAdhesion=VCAFSpaceAdhesion;
			}
		}
		else{
			if (cell->type == SCCtype){
				MinAdh=SCCECMMinAdh;
				MaxAdh=cell->MaximumSCCECMAdhesion;
			}
			else if (cell->type == VCAFtype){
				MinAdh=VCAFECMMinAdh;
				MaxAdh=cell->MaximumVCAFECMAdhesion;
			}	
			quadraticB=MaxAdh/(MinAdh-MaxAdh);
			quadraticA=quadraticB*MinAdh;
			if (cell->type == SCCtype){
				cell->CurrentSCCECMAdhesion=quadraticA/(pow(cell->MeanSurroundingECMConc,RecipPower)+quadraticB);
				cell->CurrentVCAFECMAdhesion=0.0;
				if (cell->CurrentSCCECMAdhesion<cell->MaximumSCCECMAdhesion){
					cell->CurrentSCCECMAdhesion=cell->MaximumSCCECMAdhesion;
				}
			}
			else if (cell->type == VCAFtype){
				cell->CurrentSCCECMAdhesion=0.0;
				cell->CurrentVCAFECMAdhesion=quadraticA/(pow(cell->MeanSurroundingECMConc,RecipPower)+quadraticB);
				if (cell->CurrentVCAFECMAdhesion<cell->MaximumVCAFECMAdhesion){
					cell->CurrentVCAFECMAdhesion=cell->MaximumVCAFECMAdhesion;
				}
			}
		}
	}
	
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void OrganotypicAdhesionECMConcentration::finish(){	
}

void OrganotypicAdhesionECMConcentration::update(CC3DXMLElement *_xmlData, bool _fullInitFlag){
     if(_xmlData->findElement("SCCtype")){SCCtype=_xmlData->getFirstElement("SCCtype")->getDouble();}
     else{SCCtype=1;cerr<<"SCCtype not specified, default value set to: "<<SCCtype<<endl;}
     if(_xmlData->findElement("VCAFtype")){VCAFtype=_xmlData->getFirstElement("VCAFtype")->getDouble();}
     else{VCAFtype=2;cerr<<"VCAFtype not specified, default value set to: "<<VCAFtype<<endl;}
     if(_xmlData->findElement("dimx")){dimx=_xmlData->getFirstElement("dimx")->getDouble();}
     if(_xmlData->findElement("dimy")){dimy=_xmlData->getFirstElement("dimy")->getDouble();}
     if(_xmlData->findElement("dimz")){dimz=_xmlData->getFirstElement("dimz")->getDouble();}
     if(_xmlData->findElement("SCCECMZeroConcAdhesion")){SCCECMMinAdh=_xmlData->getFirstElement("SCCECMZeroConcAdhesion")->getDouble();}
     else{SCCECMMinAdh=50;cerr<<"SCCECMZeroConcAdhesion not specified, default value set to: "<<SCCECMMinAdh<<endl;}
     if(_xmlData->findElement("VCAFECMZeroConcAdhesion")){VCAFECMMinAdh=_xmlData->getFirstElement("VCAFECMZeroConcAdhesion")->getDouble();}
     else{VCAFECMMinAdh=50;cerr<<"VCAFECMZeroConcAdhesion not specified, default value set to: "<<VCAFECMMinAdh<<endl;}
     if(_xmlData->findElement("ReciprocalPower")){RecipPower=_xmlData->getFirstElement("ReciprocalPower")->getDouble();}
     else{RecipPower=0.25;cerr<<"ReciprocalPower not specified, default value set to: "<<RecipPower<<endl;}
     if(_xmlData->findElement("SCCSpaceAdhesion")){SCCSpaceAdhesion=_xmlData->getFirstElement("SCCSpaceAdhesion")->getDouble();}
     else{SCCSpaceAdhesion=50;cerr<<"SCCSpaceAdhesion not specified, default value set to: "<<SCCSpaceAdhesion<<endl;}
     if(_xmlData->findElement("VCAFSpaceAdhesion")){VCAFSpaceAdhesion=_xmlData->getFirstElement("VCAFSpaceAdhesion")->getDouble();}
     else{VCAFSpaceAdhesion=50;cerr<<"VCAFSpaceAdhesion not specified, default value set to: "<<VCAFSpaceAdhesion<<endl;}
     if(_xmlData->findElement("ECMSpaceBoundary")){SpaceBoundary=_xmlData->getFirstElement("ECMSpaceBoundary")->getDouble();}
     else{SpaceBoundary=10000;cerr<<"ECMSpaceBoundary not specified, default value set to: "<<SpaceBoundary<<endl;}
     if(_xmlData->findElement("VariableCellAdhesionInputFileName")){CellAdhesionInputFileName=_xmlData->getFirstElement("VariableCellAdhesionInputFileName")->getText();}
     cout<<"READ stuff"<<endl;      
}

std::string OrganotypicAdhesionECMConcentration::toString(){
   return "OrganotypicAdhesionECMConcentration";
}


std::string OrganotypicAdhesionECMConcentration::steerableName(){
   return toString();
}


