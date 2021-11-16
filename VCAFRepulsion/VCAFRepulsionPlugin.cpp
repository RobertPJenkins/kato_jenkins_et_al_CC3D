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


#include "VCAFRepulsionPlugin.h"

void VCAFRepulsionPlugin::init(Simulator *simulator, CC3DXMLElement *_xmlData)
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

void VCAFRepulsionPlugin::update(CC3DXMLElement *_xmlData, bool _fullInitFlag)
{
      if(_xmlData->findElement("VCAFtype"))
        VCAFtype=_xmlData->getFirstElement("VCAFtype")->getDouble();
      if(_xmlData->findElement("LambdaVCAFrepulsion"))
        lambdaVCAFrepulsion=_xmlData->getFirstElement("LambdaVCAFrepulsion")->getDouble();
      if(_xmlData->findElement("ThresholdDiameter"))
        thres_diameter=_xmlData->getFirstElement("ThresholdDiameter")->getDouble();
      thres2=thres_diameter*thres_diameter;       
}

double VCAFRepulsionPlugin::changeEnergy(const Point3D &pt,const CellG *newCell,const CellG *oldCell) 
{
	double energy = 0.0;
	if (oldCell == newCell) return 0;
	if (newCell){
		if((int)newCell->type==VCAFtype){
			float repulsion=0.0;
			float CMpre[3];
			float CMpos[3];
			CMpre[0]=newCell->xCM/newCell->volume;
			CMpre[1]=newCell->yCM/newCell->volume;
			CMpre[2]=newCell->zCM/newCell->volume;	
			CMpos[0]=(newCell->xCM+pt.x)/(newCell->volume+1);
			CMpos[1]=(newCell->yCM+pt.y)/(newCell->volume+1);
			CMpos[2]=(newCell->zCM+pt.z)/(newCell->volume+1);
			VCAFpositionnode *conductor = newCell->VCAFListroot;
			while ( conductor!= NULL){
   				if(conductor->id!=newCell->id){
   					float dxpre=conductor->pos[0]-CMpre[0];
					float dypre=conductor->pos[1]-CMpre[1];
					float dzpre=conductor->pos[2]-CMpre[2];
					if(dxpre<thres_diameter && dypre<thres_diameter && dzpre<thres_diameter){
						float d2pre=dxpre*dxpre+dypre*dypre+dzpre*dzpre;
						if(d2pre<thres2){
							float dxpos=conductor->pos[0]-CMpos[0];
							float dypos=conductor->pos[1]-CMpos[1];
							float dzpos=conductor->pos[2]-CMpos[2];
							float d2pos=dxpos*dxpos+dypos*dypos+dzpos*dzpos;
							float dpre=pow(d2pre,(float)0.5);
							float dpos=pow(d2pos,(float)0.5);
							float repulsionpre=(-(1.0/thres_diameter)*dpre+1);
							float repulsionpos=(-(1.0/thres_diameter)*dpos+1);
							repulsion=repulsion+(repulsionpos-repulsionpre);
						}
					}
   				}
   				conductor = conductor->next;
			}
			energy=energy+repulsion*lambdaVCAFrepulsion;
		}
	}
	if (oldCell){
		if((int)oldCell->type==VCAFtype){
			float repulsion=0.0;
			float CMpre[3];
			float CMpos[3];
			CMpre[0]=oldCell->xCM/oldCell->volume;
			CMpre[1]=oldCell->yCM/oldCell->volume;
			CMpre[2]=oldCell->zCM/oldCell->volume;	
			CMpos[0]=(oldCell->xCM-pt.x)/(oldCell->volume-1);
			CMpos[1]=(oldCell->yCM-pt.y)/(oldCell->volume-1);
			CMpos[2]=(oldCell->zCM-pt.z)/(oldCell->volume-1);
			VCAFpositionnode *conductor = oldCell->VCAFListroot;
			while ( conductor!= NULL){
   				if(conductor->id!=oldCell->id){
   					float dxpre=conductor->pos[0]-CMpre[0];
					float dypre=conductor->pos[1]-CMpre[1];
					float dzpre=conductor->pos[2]-CMpre[2];
					if(dxpre<thres_diameter && dypre<thres_diameter && dzpre<thres_diameter){
						float d2pre=dxpre*dxpre+dypre*dypre+dzpre*dzpre;
						if(d2pre<thres2){
							float dxpos=conductor->pos[0]-CMpos[0];
							float dypos=conductor->pos[1]-CMpos[1];
							float dzpos=conductor->pos[2]-CMpos[2];
							float d2pos=dxpos*dxpos+dypos*dypos+dzpos*dzpos;
							float dpre=pow(d2pre,(float)0.5);
							float dpos=pow(d2pos,(float)0.5);
							float repulsionpre=(-(1.0/thres_diameter)*dpre+1);
							float repulsionpos=(-(1.0/thres_diameter)*dpos+1);
							repulsion=repulsion+(repulsionpos-repulsionpre);
						}
					}
   				}
   				conductor = conductor->next;
			}
			energy=energy+repulsion*lambdaVCAFrepulsion;
		}
	}
	return energy;
}

std::string VCAFRepulsionPlugin::steerableName()
{
	return toString();
}

std::string VCAFRepulsionPlugin::toString()
{
	return "VCAFRepulsion";
}

