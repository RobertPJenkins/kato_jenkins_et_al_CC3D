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
using namespace std;


#include "OrganotypicKinesisDirectionality.h"

OrganotypicKinesisDirectionality::OrganotypicKinesisDirectionality() :  cellFieldG(0),sim(0),potts(0) {}

OrganotypicKinesisDirectionality::~OrganotypicKinesisDirectionality() {
}


void OrganotypicKinesisDirectionality::init(Simulator *simulator, CC3DXMLElement *_xmlData) {
	potts = simulator->getPotts();
    sim=simulator;
    cellFieldG = (WatchableField3D<CellG *> *)potts->getCellFieldG();
	simulator->registerSteerableObject(this);
    update(_xmlData);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void OrganotypicKinesisDirectionality::extraInit(Simulator *simulator){

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void OrganotypicKinesisDirectionality::start(){

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void OrganotypicKinesisDirectionality::step(const unsigned int currentStep){
	CellInventory *cellInventoryPtr=& potts->getCellInventory();
	CellInventory::cellInventoryIterator cInvItr;
	CellG *cell;
	double mean=0.0;
	unsigned int cellCounter=0;
	for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
		cell=cellInventoryPtr->getCell(cInvItr);
		if(currentStep==0){
			if((int)cell->type==SCCtype){
				cell->lambdataxis=SCClambdataxis;
				double mag=SCCchemotaxisdir[0]*SCCchemotaxisdir[0]+SCCchemotaxisdir[1]*SCCchemotaxisdir[1]+SCCchemotaxisdir[2]*SCCchemotaxisdir[2];
				mag=pow(mag,0.5);
				float SCCchemotaxisnorm[3];
				SCCchemotaxisnorm[0]=SCCchemotaxisdir[0]/mag;
				SCCchemotaxisnorm[1]=SCCchemotaxisdir[1]/mag;
				SCCchemotaxisnorm[2]=SCCchemotaxisdir[2]/mag;
				float r[3];
				r[0]=rand()%2000;r[0]=r[0]/1000.0;r[0]=r[0]-1;
				r[1]=rand()%2000;r[1]=r[1]/1000.0;r[1]=r[1]-1;
				r[2]=rand()%2000;r[2]=r[2]/1000.0;r[2]=r[2]-1;
				float rmag=pow((r[0]*r[0] + r[1]*r[1] + r[2]*r[2]),(float) 0.5);
				r[0]=r[0]/rmag;r[1]=r[1]/rmag;r[2]=r[2]/rmag;
				float wsccpre=1.0-wsccrand;
				cell->taxisidir[0]=wsccpre*SCCchemotaxisnorm[0]+wsccrand*r[0];
				cell->taxisidir[1]=wsccpre*SCCchemotaxisnorm[1]+wsccrand*r[1];
				cell->taxisidir[2]=wsccpre*SCCchemotaxisnorm[2]+wsccrand*r[2];
				double taxisdirmag = cell->taxisidir[0]*cell->taxisidir[0]+cell->taxisidir[1]*cell->taxisidir[1]+cell->taxisidir[2]*cell->taxisidir[2];
				taxisdirmag=pow(taxisdirmag,(double) 0.5);
				cell->taxisidir[0]=cell->taxisidir[0]/taxisdirmag;
				cell->taxisidir[1]=cell->taxisidir[1]/taxisdirmag;
				cell->taxisidir[2]=cell->taxisidir[2]/taxisdirmag;
				cell->previousCM[0]=cell->xCOM;
				cell->previousCM[1]=cell->yCOM;
				cell->previousCM[2]=cell->zCOM;
			}
			else if((int)cell->type==VCAFtype){
				cell->lambdataxis=VCAFlambdataxis_min;
				double mag=VCAFchemotaxisdir[0]*VCAFchemotaxisdir[0]+VCAFchemotaxisdir[1]*VCAFchemotaxisdir[1]+VCAFchemotaxisdir[2]*VCAFchemotaxisdir[2];
				mag=pow(mag,0.5);
				float VCAFchemotaxisnorm[3];
				VCAFchemotaxisnorm[0]=VCAFchemotaxisdir[0]/mag;
				VCAFchemotaxisnorm[1]=VCAFchemotaxisdir[1]/mag;
				VCAFchemotaxisnorm[2]=VCAFchemotaxisdir[2]/mag;
				float r[3];
				r[0]=rand()%2000;r[0]=r[0]/1000.0;r[0]=r[0]-1;
				r[1]=rand()%2000;r[1]=r[1]/1000.0;r[1]=r[1]-1;
				r[2]=rand()%2000;r[2]=r[2]/1000.0;r[2]=r[2]-1;
				float rmag=pow((r[0]*r[0] + r[1]*r[1] + r[2]*r[2]),(float) 0.5);
				r[0]=r[0]/rmag;r[1]=r[1]/rmag;r[2]=r[2]/rmag;				
				float wrandinitial=1.0-wchem;
				cell->taxisidir[0]=wchem*VCAFchemotaxisnorm[0]+wrandinitial*r[0];
				cell->taxisidir[1]=wchem*VCAFchemotaxisnorm[1]+wrandinitial*r[1];
				cell->taxisidir[2]=wchem*VCAFchemotaxisnorm[2]+wrandinitial*r[2];
				double taxisdirmag = cell->taxisidir[0]*cell->taxisidir[0]+cell->taxisidir[1]*cell->taxisidir[1]+cell->taxisidir[2]*cell->taxisidir[2];
				taxisdirmag=pow(taxisdirmag,(double) 0.5);
				cell->taxisidir[0]=cell->taxisidir[0]/taxisdirmag;
				cell->taxisidir[1]=cell->taxisidir[1]/taxisdirmag;
				cell->taxisidir[2]=cell->taxisidir[2]/taxisdirmag;
				cell->previousCM[0]=cell->xCOM;
				cell->previousCM[1]=cell->yCOM;
				cell->previousCM[2]=cell->zCOM;
			}
		}
		else if((int)cell->type==SCCtype){
			cell->lambdataxis=SCClambdataxis;
			double mag=SCCchemotaxisdir[0]*SCCchemotaxisdir[0]+SCCchemotaxisdir[1]*SCCchemotaxisdir[1]+SCCchemotaxisdir[2]*SCCchemotaxisdir[2];
			mag=pow(mag,0.5);
			float SCCchemotaxisnorm[3];
			SCCchemotaxisnorm[0]=SCCchemotaxisdir[0]/mag;
			SCCchemotaxisnorm[1]=SCCchemotaxisdir[1]/mag;
			SCCchemotaxisnorm[2]=SCCchemotaxisdir[2]/mag;
			float r[3];
			r[0]=rand()%2000;r[0]=r[0]/1000.0;r[0]=r[0]-1;
			r[1]=rand()%2000;r[1]=r[1]/1000.0;r[1]=r[1]-1;
			r[2]=rand()%2000;r[2]=r[2]/1000.0;r[2]=r[2]-1;
			float rmag=pow((r[0]*r[0] + r[1]*r[1] + r[2]*r[2]),(float) 0.5);
			r[0]=r[0]/rmag;r[1]=r[1]/rmag;r[2]=r[2]/rmag;
			float wsccpre=1.0-wsccrand;
			cell->taxisidir[0]=wsccpre*SCCchemotaxisnorm[0]+wsccrand*r[0];
			cell->taxisidir[1]=wsccpre*SCCchemotaxisnorm[1]+wsccrand*r[1];
			cell->taxisidir[2]=wsccpre*SCCchemotaxisnorm[2]+wsccrand*r[2];
			double taxisdirmag = cell->taxisidir[0]*cell->taxisidir[0]+cell->taxisidir[1]*cell->taxisidir[1]+cell->taxisidir[2]*cell->taxisidir[2];
			taxisdirmag=pow(taxisdirmag,(double) 0.5);
			cell->taxisidir[0]=cell->taxisidir[0]/taxisdirmag;
			cell->taxisidir[1]=cell->taxisidir[1]/taxisdirmag;
			cell->taxisidir[2]=cell->taxisidir[2]/taxisdirmag;
			if (cell->just_divided_taxis_label == true ){
				cell->previousCM[0]=cell->xCOM;
				cell->previousCM[1]=cell->yCOM;
				cell->previousCM[2]=cell->zCOM;
				cell->just_divided_taxis_label=false;
			}
		}
		else if((int)cell->type==VCAFtype){
			cell->lambdataxis=cell->lambdataxis-rateoflambdadecay;
			double repulsion=0.0;
			CellInventory::cellInventoryIterator cInvItr2;
			CellG *cell2;
			for(cInvItr2=cellInventoryPtr->cellInventoryBegin() ; cInvItr2 !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr2 ){
				cell2=cellInventoryPtr->getCell(cInvItr2);
				if (((int)cell2->type==SCCtype || (int)cell2->type==VCAFtype) && cell2->id!=cell->id){
					float x=cell2->xCOM,y=cell2->yCOM,z=cell2->zCOM;
					float dx=(cell->xCOM-x),dy=(cell->yCOM-y),dz=(cell->zCOM-z);
					float d2=dx*dx+dy*dy+dz*dz;
					if(d2<thres2){
						float d=pow(d2,(float)0.5);
						repulsion=repulsion+(-(1.0/thres_diameter)*d+1);
					}
				}
			}
			cell->lambdataxis=cell->lambdataxis+repulsionscalefactor*repulsion;
			if(cell->lambdataxis<VCAFlambdataxis_min){
				cell->lambdataxis=VCAFlambdataxis_min;
			}
			else if(cell->lambdataxis>VCAFlambdataxis_max){
				cell->lambdataxis=VCAFlambdataxis_max;
			}
			double mag=VCAFchemotaxisdir[0]*VCAFchemotaxisdir[0]+VCAFchemotaxisdir[1]*VCAFchemotaxisdir[1]+VCAFchemotaxisdir[2]*VCAFchemotaxisdir[2];
			mag=pow(mag,0.5);
			float VCAFchemotaxisnorm[3];
			VCAFchemotaxisnorm[0]=VCAFchemotaxisdir[0]/mag;
			VCAFchemotaxisnorm[1]=VCAFchemotaxisdir[1]/mag;
			VCAFchemotaxisnorm[2]=VCAFchemotaxisdir[2]/mag;
			float r[3];
			r[0]=rand()%2000;r[0]=r[0]/1000.0;r[0]=r[0]-1;
			r[1]=rand()%2000;r[1]=r[1]/1000.0;r[1]=r[1]-1;
			r[2]=rand()%2000;r[2]=r[2]/1000.0;r[2]=r[2]-1;
			float rmag=pow((r[0]*r[0] + r[1]*r[1] + r[2]*r[2]),(float) 0.5);
			r[0]=r[0]/rmag;r[1]=r[1]/rmag;r[2]=r[2]/rmag;
			float vprev[3]={cell->xCOM-cell->previousCM[0],cell->yCOM-cell->previousCM[1],cell->zCOM-cell->previousCM[2]};
			float vprevmag=pow((vprev[0]*vprev[0] + vprev[1]*vprev[1] + vprev[2]*vprev[2]),(float)0.5);
			if(vprevmag==0){
				vprev[0]=0.0;vprev[1]=0.0;vprev[2]=0.0;
			}
			else{
				vprev[0]=vprev[0]/vprevmag;vprev[1]=vprev[1]/vprevmag;vprev[2]=vprev[2]/vprevmag;
			}
			cell->previousCM[0]=cell->xCOM;cell->previousCM[1]=cell->yCOM;cell->previousCM[2]=cell->zCOM;
			float wpre=1.0-wrand-wprev-wchem;
			cell->taxisidir[0]=wpre*cell->taxisidir[0]+wrand*r[0]+wprev*vprev[0]+wchem*VCAFchemotaxisnorm[0];
			cell->taxisidir[1]=wpre*cell->taxisidir[1]+wrand*r[1]+wprev*vprev[1]+wchem*VCAFchemotaxisnorm[1];
			cell->taxisidir[2]=wpre*cell->taxisidir[2]+wrand*r[2]+wprev*vprev[2]+wchem*VCAFchemotaxisnorm[2];
			double taxisdirmag = cell->taxisidir[0]*cell->taxisidir[0]+cell->taxisidir[1]*cell->taxisidir[1]+cell->taxisidir[2]*cell->taxisidir[2];
			taxisdirmag=pow(taxisdirmag,(double) 0.5);
			cell->taxisidir[0]=cell->taxisidir[0]/taxisdirmag;
			cell->taxisidir[1]=cell->taxisidir[1]/taxisdirmag;
			cell->taxisidir[2]=cell->taxisidir[2]/taxisdirmag;
		}
	}
}






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void OrganotypicKinesisDirectionality::update(CC3DXMLElement *_xmlData, bool _fullInitFlag){
     if(_xmlData->findElement("SCCtype")){SCCtype=_xmlData->getFirstElement("SCCtype")->getDouble();}
     else{SCCtype=1;cerr<<"SCCtype not specified, default value set to: "<<SCCtype<<endl;}
     if(_xmlData->findElement("VCAFtype")){VCAFtype=_xmlData->getFirstElement("VCAFtype")->getDouble();}
     else{VCAFtype=2;cerr<<"VCAFtype not specified, default value set to: "<<VCAFtype<<endl;}
     if(_xmlData->findElement("SCCLambdaTaxis"))
        SCClambdataxis=_xmlData->getFirstElement("SCCLambdaTaxis")->getDouble();
     if(_xmlData->findElement("VCAFLambdaTaxis_min"))
        VCAFlambdataxis_min=_xmlData->getFirstElement("VCAFLambdaTaxis_min")->getDouble();
     if(_xmlData->findElement("VCAFLambdaTaxis_max"))
        VCAFlambdataxis_max=_xmlData->getFirstElement("VCAFLambdaTaxis_max")->getDouble();
     if(_xmlData->findElement("SCCchemotaxisdir_x"))
        SCCchemotaxisdir[0]=_xmlData->getFirstElement("SCCchemotaxisdir_x")->getDouble();
     if(_xmlData->findElement("SCCchemotaxisdir_y"))
        SCCchemotaxisdir[1]=_xmlData->getFirstElement("SCCchemotaxisdir_y")->getDouble();
     if(_xmlData->findElement("SCCchemotaxisdir_z"))
        SCCchemotaxisdir[2]=_xmlData->getFirstElement("SCCchemotaxisdir_z")->getDouble();
     if(_xmlData->findElement("VCAFchemotaxisdir_x")){
        VCAFchemotaxisdir[0]=_xmlData->getFirstElement("VCAFchemotaxisdir_x")->getDouble();
     }
     else{
      	VCAFchemotaxisdir[0]=0.0;cerr<<"VCAFchemotaxisdir_x not specified, default value set to: "<<VCAFchemotaxisdir[0]<<"."<<endl;
     }
     if(_xmlData->findElement("VCAFchemotaxisdir_y")){
        VCAFchemotaxisdir[1]=_xmlData->getFirstElement("VCAFchemotaxisdir_y")->getDouble();
     }
     else{
      	VCAFchemotaxisdir[1]=0.0;cerr<<"VCAFchemotaxisdir_y not specified, default value set to: "<<VCAFchemotaxisdir[1]<<"."<<endl;
     }
     if(_xmlData->findElement("VCAFchemotaxisdir_z")){
        VCAFchemotaxisdir[2]=_xmlData->getFirstElement("VCAFchemotaxisdir_z")->getDouble();
     }
     else{
      	VCAFchemotaxisdir[2]=-1.0;cerr<<"VCAFchemotaxisdir_z not specified, default value set to: "<<VCAFchemotaxisdir[2]<<"."<<endl;
     }      
     if(_xmlData->findElement("LambdaDecayRate"))
        rateoflambdadecay=_xmlData->getFirstElement("LambdaDecayRate")->getDouble();    
     if(_xmlData->findElement("RepulsionScalingFactor"))
        repulsionscalefactor=_xmlData->getFirstElement("RepulsionScalingFactor")->getDouble();
     if(_xmlData->findElement("Weight_RandomNoise"))
        wrand=_xmlData->getFirstElement("Weight_RandomNoise")->getDouble();
     if(_xmlData->findElement("Weight_PreviousVelocity"))
        wprev=_xmlData->getFirstElement("Weight_PreviousVelocity")->getDouble();
     if(_xmlData->findElement("Weight_Chemotactic")){ 
     	wchem=_xmlData->getFirstElement("Weight_Chemotactic")->getDouble();
     }
     else{
     	wchem=0.0;cerr<<"Weight_Chemotactic not specified, default value set to: "<<wchem<<"."<<endl;
     }   
     if(_xmlData->findElement("Weight_RandomSCCNoise"))
        wsccrand=_xmlData->getFirstElement("Weight_RandomSCCNoise")->getDouble();
     if(_xmlData->findElement("ThresholdDiameter"))
        thres_diameter=_xmlData->getFirstElement("ThresholdDiameter")->getDouble();
      thres2=thres_diameter*thres_diameter;     
}

std::string OrganotypicKinesisDirectionality::toString(){
   return "OrganotypicKinesisDirectionality";
}


std::string OrganotypicKinesisDirectionality::steerableName(){
   return toString();
}


