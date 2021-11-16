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


#include <CompuCell3D/Potts3D/CellInventory.h>

#include <iostream>
using namespace std;


#include "VCAFPositionListing.h"

VCAFPositionListing::VCAFPositionListing() :  cellFieldG(0),sim(0),potts(0),VCAFtype(2.0),VCAFListroot(NULL) {}

VCAFPositionListing::~VCAFPositionListing() {
}


void VCAFPositionListing::init(Simulator *simulator, CC3DXMLElement *_xmlData) {
	potts = simulator->getPotts();
    sim=simulator;
    cellFieldG = (WatchableField3D<CellG *> *)potts->getCellFieldG();
	simulator->registerSteerableObject(this);
  	update(_xmlData);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VCAFPositionListing::extraInit(Simulator *simulator){

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void VCAFPositionListing::start(){
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VCAFPositionListing::step(const unsigned int currentStep){
	CellInventory *cellInventoryPtr=& potts->getCellInventory();
	CellInventory::cellInventoryIterator cInvItr;
	CellG *cell;
	if(currentStep==0){
		VCAFListroot = new VCAFpositionnode;
		VCAFpositionnode *conductor; 
		conductor = VCAFListroot;
		bool initiatedfirstcell=false;
   		for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
	      		cell=cellInventoryPtr->getCell(cInvItr);
      			if((int)cell->type==VCAFtype){
      				if(!initiatedfirstcell){
      					initiatedfirstcell=true;
	      			}
      				else{
	      				conductor->next = new VCAFpositionnode;
	      				conductor = conductor->next;
		      		}
		      		cell->VCAFListroot=VCAFListroot;
      				conductor->id=cell->id;
      				conductor->pos[0]=cell->xCM/cell->volume;
      				conductor->pos[1]=cell->yCM/cell->volume;
	      			conductor->pos[2]=cell->zCM/cell->volume;
      			}
      		}
	}
	else{
		for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
			cell=cellInventoryPtr->getCell(cInvItr);
			if((int)cell->type==VCAFtype){
				VCAFpositionnode *conductor = VCAFListroot;
				while ( conductor!= NULL ){
   					if(conductor->id==cell->id){
   						conductor->pos[0]=cell->xCM/cell->volume;
						conductor->pos[1]=cell->yCM/cell->volume;
						conductor->pos[2]=cell->zCM/cell->volume;
   					}
   					conductor = conductor->next;
   				}
			}		
	   	}
	}
}






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VCAFPositionListing::finish(){
	VCAFpositionnode *conductor;
	conductor = VCAFListroot;
	
	while (conductor!=NULL){
		VCAFpositionnode *current=conductor;
		conductor=current->next;
		delete current;
	}	
}


void VCAFPositionListing::update(CC3DXMLElement *_xmlData, bool _fullInitFlag){
	 if(_xmlData->findElement("VCAFtype")) VCAFtype=_xmlData->getFirstElement("VCAFtype")->getDouble();
}

std::string VCAFPositionListing::toString(){
   return "VCAFPositionListing";
}


std::string VCAFPositionListing::steerableName(){
   return toString();
}


