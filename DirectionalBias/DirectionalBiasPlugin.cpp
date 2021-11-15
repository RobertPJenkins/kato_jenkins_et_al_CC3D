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


#include "DirectionalBiasPlugin.h"

void DirectionalBiasPlugin::init(Simulator *simulator, CC3DXMLElement *_xmlData)
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

void DirectionalBiasPlugin::update(CC3DXMLElement *_xmlData, bool _fullInitFlag)
{
     if(_xmlData->findElement("SCCtype")){SCCtype=_xmlData->getFirstElement("SCCtype")->getDouble();}
     else{SCCtype=1;cerr<<"SCCtype not specified, defaul t value set to: "<<SCCtype<<endl;}
     if(_xmlData->findElement("VCAFtype")){VCAFtype=_xmlData->getFirstElement("VCAFtype")->getDouble();}
     else{VCAFtype=2;cerr<<"VCAFtype not specified, defaul t value set to: "<<VCAFtype<<endl;}
}

double DirectionalBiasPlugin::changeEnergy(const Point3D &pt,const CellG *newCell,const CellG *oldCell) 
{
	double energy = 0.0;
	if (oldCell == newCell) return 0;

	if (newCell){
		if((int)newCell->type==SCCtype || (int)newCell->type==VCAFtype){
			float Cxpre=newCell->xCM/newCell->volume; 
			float Cypre=newCell->yCM/newCell->volume;
			float Czpre=newCell->zCM/newCell->volume;
			float Cxpost=(newCell->xCM+pt.x)/(newCell->volume+1);
			float Cypost=(newCell->yCM+pt.y)/(newCell->volume+1);
			float Czpost=(newCell->zCM+pt.z)/(newCell->volume+1);
			float dCM[3]={Cxpost-Cxpre,Cypost-Cypre,Czpost-Czpre};
			float dCMmag=pow((dCM[0]*dCM[0]+dCM[1]*dCM[1]+dCM[2]*dCM[2]), (float) 0.5);
			dCM[0]=dCM[0]/dCMmag;dCM[1]=dCM[1]/dCMmag;dCM[2]=dCM[2]/dCMmag;
			float dotp=dCM[0]*newCell->taxisidir[0]+dCM[1]*newCell->taxisidir[1]+dCM[2]*newCell->taxisidir[2];
			energy=energy-newCell->lambdataxis*dotp;
		}
	}
	if (oldCell){
		if((int)oldCell->type==SCCtype || (int)oldCell->type==VCAFtype){
			float Cxpre=oldCell->xCM/oldCell->volume; 
			float Cypre=oldCell->yCM/oldCell->volume;
			float Czpre=oldCell->zCM/oldCell->volume;
			float Cxpost=(oldCell->xCM-pt.x)/(oldCell->volume-1);
			float Cypost=(oldCell->yCM-pt.y)/(oldCell->volume-1);
			float Czpost=(oldCell->zCM-pt.z)/(oldCell->volume-1);
			float dCM[3]={Cxpost-Cxpre,Cypost-Cypre,Czpost-Czpre};
			float dCMmag=pow((dCM[0]*dCM[0]+dCM[1]*dCM[1]+dCM[2]*dCM[2]), (float) 0.5);
			dCM[0]=dCM[0]/dCMmag;dCM[1]=dCM[1]/dCMmag;dCM[2]=dCM[2]/dCMmag;
			float dotp=dCM[0]*oldCell->taxisidir[0]+dCM[1]*oldCell->taxisidir[1]+dCM[2]*oldCell->taxisidir[2];
			energy=energy-oldCell->lambdataxis*dotp;
		}
	}
	return energy;
}

std::string DirectionalBiasPlugin::steerableName()
{
	return toString();
}

std::string DirectionalBiasPlugin::toString()
{
	return "DirectionalBias";
}

