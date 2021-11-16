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



#include <CompuCell3D/Simulator.h>
#include <CompuCell3D/Potts3D/Potts3D.h>
#include <CompuCell3D/Automaton/Automaton.h>
#include <CompuCell3D/Potts3D/Cell.h>
using namespace CompuCell3D;

#include <iostream>
#include <string>
#include <algorithm>
using namespace std;


#include "ECMPenetrationPlugin.h"

void ECMPenetrationPlugin::init(Simulator *simulator, CC3DXMLElement *_xmlData)
{
	potts = simulator->getPotts();
	bool pluginAlreadyRegisteredFlag;
	Plugin *plugin=Simulator::pluginManager.get("VolumeTracker",&pluginAlreadyRegisteredFlag);
	if(!pluginAlreadyRegisteredFlag)
		plugin->init(simulator);
	potts->registerEnergyFunctionWithName(this,toString());
	xmlData=_xmlData;
	simulator->registerSteerableObject(this);
	update(_xmlData);
}

void ECMPenetrationPlugin::update(CC3DXMLElement *_xmlData, bool _fullInitFlag)
     if(_xmlData->findElement("ECMDensityThreshold")){ecmDensityThreshold=_xmlData->getFirstElement("ECMDensityThreshold")->getDouble();}
     else{ecmDensityThreshold=0.9;cerr<<"ECMDensityThreshold not specified, default value set to: "<<ecmDensityThreshold<<endl;}
     if(_xmlData->findElement("MinimumDegThreshold")){minDeg=_xmlData->getFirstElement("MinimumDegThreshold")->getDouble();}
     else{minDeg=0.0001;cerr<<"MinimumDegThreshold not specified, default value set to: "<<minDeg<<endl;}
     if(_xmlData->findElement("MinimumPushThreshold")){minPush=_xmlData->getFirstElement("MinimumPushThreshold")->getDouble();}
     else{minPush=0.001;cerr<<"MinimumPushThreshold not specified, default value set to: "<<minPush<<endl;}
}


double ECMPenetrationPlugin::changeEnergy(const Point3D &pt,const CellG *newCell,const CellG *oldCell) 
{
	double energy = 0.0;
	if (oldCell == newCell) return 0;
	CellG *cellforECMfield;
	CellInventory *cellInventoryPtr=& potts->getCellInventory();
	bool foundPointerToECMField = false;
	if (newCell && newCell->steppablesinitialised==1){
		if (newCell->just_divided){
			if(oldCell && oldCell->steppablesinitialised==1 && oldCell->just_divided==false){
				CellG *cell2;
				for(CellInventory::cellInventoryIterator cInvItr2=cellInventoryPtr->cellInventoryBegin() ; cInvItr2 !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr2 ){	
					cell2=cellInventoryPtr->getCell(cInvItr2);
					if (cell2->id == newCell->id){
						cell2->ECMfieldbyCell=oldCell->ECMfieldbyCell;
						cell2->just_divided=false;
						foundPointerToECMField = true;
						break;
					}
				}
			}
			else{
				for(CellInventory::cellInventoryIterator cInvItr2=cellInventoryPtr->cellInventoryBegin() ; cInvItr2 !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr2 ){
					cellforECMfield=cellInventoryPtr->getCell(cInvItr2);
					if (cellforECMfield->just_divided==false){
						CellG *cell2;
						for(CellInventory::cellInventoryIterator cInvItr2=cellInventoryPtr->cellInventoryBegin() ; cInvItr2 !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr2 ){	
							cell2=cellInventoryPtr->getCell(cInvItr2);
							if (cell2->id == newCell->id){
								cell2->ECMfieldbyCell=cellforECMfield->ECMfieldbyCell;
								cell2->just_divided=false;
								foundPointerToECMField = true;
								break;
							}
						}
						break;
					}
				}	
			}
			if (!foundPointerToECMField){cerr<<"COULD NOT FIND A POINTER TO ECM FIELD!!! WILL CRASH!!!"<<endl;}
		}
		bool inNewCellBoundaryList=false;
		int flipx=pt.x;int flipy=pt.y;int flipz=pt.z;
		int nx = newCell->canremodel_x.size();
		for (int ii=0;ii<nx;ii++){
			if(newCell->canremodel_x[ii]==flipx && newCell->canremodel_y[ii]==flipy && newCell->canremodel_z[ii]==flipz){
				inNewCellBoundaryList=true;
				break;
			}
		}
		
		if (((newCell->kdeg>minDeg)||(newCell->wpush>minPush))&&(inNewCellBoundaryList==true)){
			energy=energy+newCell->ECMfieldbyCell[pt.x][pt.y][pt.z]*newCell->lambdaECMPenetration;
		}
		else{
			if (newCell->ECMfieldbyCell[pt.x][pt.y][pt.z]<ecmDensityThreshold){
				energy=energy+newCell->ECMfieldbyCell[pt.x][pt.y][pt.z]*newCell->lambdaECMPenetration;
			}
			else{
				energy=energy+1e10;
			}
		}
	}
	if (oldCell && oldCell->steppablesinitialised==1){
		if (oldCell->just_divided){
			if (newCell && newCell->steppablesinitialised==1 && newCell->just_divided==false){
				CellG *cell2;
				for(CellInventory::cellInventoryIterator cInvItr2=cellInventoryPtr->cellInventoryBegin() ; cInvItr2 !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr2 ){	
					cell2=cellInventoryPtr->getCell(cInvItr2);
					if (cell2->id == oldCell->id){
						cell2->ECMfieldbyCell=newCell->ECMfieldbyCell;
						cell2->just_divided=false;
						foundPointerToECMField = true;
						break;
					}
				}
			}
			else{
				for(CellInventory::cellInventoryIterator cInvItr2=cellInventoryPtr->cellInventoryBegin() ; cInvItr2 !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr2 ){
					cellforECMfield=cellInventoryPtr->getCell(cInvItr2);
					if (cellforECMfield->just_divided==false){
						CellG *cell2;
						for(CellInventory::cellInventoryIterator cInvItr2=cellInventoryPtr->cellInventoryBegin() ; cInvItr2 !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr2 ){	
							cell2=cellInventoryPtr->getCell(cInvItr2);
							if (cell2->id == oldCell->id){
								cell2->ECMfieldbyCell=cellforECMfield->ECMfieldbyCell;
								cell2->just_divided=false;
								foundPointerToECMField = true;
								break;
							}
						}
						break;
					}
				}
			}
			if (!foundPointerToECMField){cerr<<"COULD NOT FIND A POINTER TO ECM FIELD!!! WILL CRASH!!!"<<endl;}
		}
		bool inOldCellBoundaryList=false;
		int flipx=pt.x;int flipy=pt.y;int flipz=pt.z;
		int nx = oldCell->canremodel_x.size();
		for (int ii=0;ii<nx;ii++){
			if(oldCell->canremodel_x[ii]==flipx && oldCell->canremodel_y[ii]==flipy && oldCell->canremodel_z[ii]==flipz){
				inOldCellBoundaryList=true;
				break;
			}
		}
		
		if (((oldCell->kdeg>minDeg)||(oldCell->wpush>minPush))&&(inOldCellBoundaryList==true)){
			energy=energy-oldCell->ECMfieldbyCell[pt.x][pt.y][pt.z]*oldCell->lambdaECMPenetration;
		}
		else{
			if (oldCell->ECMfieldbyCell[pt.x][pt.y][pt.z]<ecmDensityThreshold){
				energy=energy-oldCell->ECMfieldbyCell[pt.x][pt.y][pt.z]*oldCell->lambdaECMPenetration;
			}
			else{
				energy=energy-1e10;
			}
		}
	}
	return energy;
}




std::string ECMPenetrationPlugin::steerableName()
{
	return toString();
}

std::string ECMPenetrationPlugin::toString()
{
	return "ECMPenetration";
}

