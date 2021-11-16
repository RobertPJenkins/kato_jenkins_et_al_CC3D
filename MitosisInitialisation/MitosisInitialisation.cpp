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
#include <math.h>
#include <algorithm> 
#include <numeric>
#include <vector>


using namespace CompuCell3D;

#include <iostream>
using namespace std;
#include <fstream>
#include <iomanip>
#include <sstream>
//This is a random number generator written by Xioafen Li.
#include "RandomGenerator.h"
#include <CompuCell3D/Potts3D/CellInventory.h>



#include "MitosisInitialisation.h"
#include "CompuCell3D/plugins/PixelTracker/PixelTrackerPlugin.h"
#include "CompuCell3D/plugins/PixelTracker/PixelTracker.h"



MitosisInitialisation::MitosisInitialisation() :  cellFieldG(0),sim(0),potts(0) {}
MitosisInitialisation::~MitosisInitialisation() {
}


void MitosisInitialisation::init(Simulator *simulator, CC3DXMLElement *_xmlData) {
	potts = simulator->getPotts();
	sim=simulator;
	cellFieldG = (WatchableField3D<CellG *> *)potts->getCellFieldG();
	simulator->registerSteerableObject(this);
	update(_xmlData);
	bool pluginAlreadyRegisteredFlag;
	pixelTrackerPlugin=(PixelTrackerPlugin*)Simulator::pluginManager.get("PixelTracker",&pluginAlreadyRegisteredFlag);
	if(!pluginAlreadyRegisteredFlag)
		pixelTrackerPlugin->init(simulator);
	pixelTrackerAccessorPtr=pixelTrackerPlugin->getPixelTrackerAccessorPtr();
}

void MitosisInitialisation::extraInit(Simulator *simulator){

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MitosisInitialisation::start(){
	cout<<"within the start function of MitosisInitialisation"<<endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MitosisInitialisation::step(const unsigned int currentStep){
	CellInventory *cellInventoryPtr=& potts->getCellInventory();
	CellInventory::cellInventoryIterator cInvItr;
	CellG *cell;
	const char *distributionname;
	distributionname=DivisionTimeDistribution.c_str();
	const char *Dis1 = "Exponential";
	const char *Dis2 = "Erlang";
	const char *SurfVolRatio;
	SurfVolRatio=SurfacevolumeRatio.c_str();
	const char *Ratio1 = "Sphere";
	const char *Ratio2 = "User";
	
	if(currentStep==0){
		for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
			cell=cellInventoryPtr->getCell(cInvItr);
			if (cell->type == SCCtype){
				cell->InitialConditionCell=true;
				cell->targetVolume=SCCVolume;
				cell->lambdaVolume=SCCLambdaVolume;
				if(strcmp(SurfVolRatio, Ratio1)==0){
					double radius=pow(3.0*(double)cell->targetVolume/4.0/3.14,1/3.0);	
					double SurfVolRatioMultiplier=1.0;
					cell->targetSurface=SurfVolRatioMultiplier*4*3.14*pow(radius,2);
				}
				else if(strcmp(SurfVolRatio, Ratio2)==0){
					cell->targetSurface=SCCSurface;
				}
				cell->lambdaSurface=SCCLambdaSurface;
				double t_rand=RandomGenerator::Obj().getUniformRV( 0.0, 1.0 );
				cell->InitialTimeToMitosis = (int) (t_rand*InitialTimeToMitosisSCC);
				cell->OriginalTargetVolume = cell->targetVolume;
				cell->OriginalTargetSurface = cell->targetSurface;
			}
			else if(cell->type== VCAFtype){
				cell->targetVolume=VCAFVolume;
				cell->lambdaVolume=VCAFLambdaVolume;
				cell->targetSurface=VCAFSurface;
				cell->lambdaSurface=VCAFLambdaSurface;
			}
		}
	}
	else{
		for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
			cell=cellInventoryPtr->getCell(cInvItr);
			if (cell->InitialConditionCell==true){
				cell->InitialTimeToMitosis--;
				if (cell->InitialTimeToMitosis==0){
					cell->InitialConditionCell=false;
					if (cell->type == SCCtype){
						if(strcmp(distributionname, Dis1)==0){
							double t_rand=RandomGenerator::Obj().getUniformRV( 0.0, 1.0 );
							double t_rand_exp=-log(t_rand); 
							cell->TimeToMitosisCheck = (int) (t_rand_exp*MeanTimeToMitosisCheckSCC);
						}
						else if(strcmp(distributionname, Dis2)==0){
							double t_rand=1.0;
							for( int a = 1; a < (int)ErlangShapeParameter+1; a = a + 1 ){
								double t_randsingle=RandomGenerator::Obj().getUniformRV( 0.0, 1.0 );
								t_rand=t_rand*t_randsingle;
							}
							double t_rand_exp=-log(t_rand)/ErlangShapeParameter;
							cell->TimeToMitosisCheck = (int) (t_rand_exp*MeanTimeToMitosisCheckSCC);
						}
					}
				}
			}
			else if(cell->InitialConditionCell==false){
				double cellvolumeplus=1.1*cell->volume;
				cell->mitosistimer++;
				cell->MaxCellVolume=cell->OriginalTargetVolume*VolumeMultiplierForGrowth;
				double volumegradient=cell->OriginalTargetVolume+(VolumeMultiplierForGrowth-1)*cell->OriginalTargetVolume/cell->TimeToMitosisCheck*cell->mitosistimer;
				cell->targetVolume=(int)(min(min(volumegradient,max(cellvolumeplus,double(cell->OriginalTargetVolume))),double(cell->MaxCellVolume)));
				if(strcmp(SurfVolRatio, Ratio1)==0){
					double radius=pow(3.0*(double)cell->targetVolume/4.0/3.14,1/3.0);	
					double SurfVolRatioMultiplier=1.0;
					cell->targetSurface=SurfVolRatioMultiplier*4*3.14*pow(radius,2);
				}
				else if(strcmp(SurfVolRatio, Ratio2)==0){
					cell->targetSurface=(int)(cell->targetVolume/cell->OriginalTargetVolume*cell->OriginalTargetSurface);
				}
				if(cell->volume>cell->MaxCellVolume*DivisionThresholdMultiplier){
					cell->DivideThisCellInPython=true;
					cell->targetVolume=cell->OriginalTargetVolume;
					cell->targetSurface=cell->OriginalTargetSurface;
					cell->mitosistimer=0;
					if (cell->type == SCCtype){
					
						if(strcmp(distributionname, Dis1)==0){
							double t_rand=RandomGenerator::Obj().getUniformRV( 0.0, 1.0 );
							double t_rand_exp=-log(t_rand); 
							cell->TimeToMitosisCheck = (int) (t_rand_exp*MeanTimeToMitosisCheckSCC);
							double t_rand2=RandomGenerator::Obj().getUniformRV( 0.0, 1.0 );
							double t_rand_exp2=-log(t_rand2); 
							cell->DaughterCellTimer = (int) (t_rand_exp2*MeanTimeToMitosisCheckSCC);
						}
						else if(strcmp(distributionname, Dis2)==0){
							double t_rand=1.0;
							for( int a = 1; a < (int)ErlangShapeParameter+1; a = a + 1 ){
								double t_randsingle=RandomGenerator::Obj().getUniformRV( 0.0, 1.0 );
								t_rand=t_rand*t_randsingle;
							}
							double t_rand_exp=-log(t_rand)/ErlangShapeParameter;
							cell->TimeToMitosisCheck = (int) (t_rand_exp*MeanTimeToMitosisCheckSCC);
						
							double t_rand2=1.0;
							for( int a = 1; a < (int)ErlangShapeParameter+1; a = a + 1 ){
								double t_randsingle2=RandomGenerator::Obj().getUniformRV( 0.0, 1.0 );
								t_rand2=t_rand2*t_randsingle2;
							}
							double t_rand_exp2=-log(t_rand2)/ErlangShapeParameter;
							cell->DaughterCellTimer = (int)(t_rand_exp2*MeanTimeToMitosisCheckSCC);	
						}				
					}
				}
			}
		}
	}
	if(currentStep==TimeToWriteCells){
		const char *outputname_pif;
		outputname_pif=InitialCellOutputFileName.c_str();
		InitialCellOutputFile.open(outputname_pif,ofstream::out);
		if (!InitialCellOutputFile.is_open()){cerr<<"could not open "<<InitialCellOutputFileName<<"!!"<<endl;}
		
		const char *outputname_time;
		outputname_time=TimeToCellDivisionFileName.c_str();
		TimeToCellDivisionFile.open(outputname_time,ofstream::out);
		if (!TimeToCellDivisionFile.is_open()){cerr<<"could not open "<<TimeToCellDivisionFile<<"!!"<<endl;}
		for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
			cell=cellInventoryPtr->getCell(cInvItr);
			if(cell->type == SCCtype){
				TimeToCellDivisionFile<<setw(5)<<cell->id<<setw(15)<<cell->mitosistimer<<setw(15)<<cell->TimeToMitosisCheck<<setw(10)<<cell->MaxCellVolume<<setw(10)<<cell->targetVolume<<setw(10)<<cell->targetSurface<<endl;
            }
            set<PixelTrackerData> cellPixels=pixelTrackerAccessorPtr->get(cell->extraAttribPtr)->pixelSet;
			for(set<PixelTrackerData>::iterator sitr=cellPixels.begin() ; sitr != cellPixels.end() ;++sitr){
				int xx=sitr->pixel.x;
				int yy=sitr->pixel.y;
				int zz=sitr->pixel.z;
				int ct = cell->type;
				if(ct==1){
					 InitialCellOutputFile <<setw(5)<<cell->id<<setw(7)<<"SCC"<<setw(5)<< xx<<setw(5)<< xx<<setw(5)<< yy <<setw(5)<< yy<<setw(5)<<zz<<setw(5)<<zz<<endl;
				}
				else if(ct==2){
					InitialCellOutputFile <<setw(5)<<cell->id<<setw(7)<<"VCAF"<<setw(5)<< xx<<setw(5)<< xx<<setw(5)<< yy <<setw(5)<< yy<<setw(5)<<zz<<setw(5)<<zz<<endl;
				}
			}
        }
        InitialCellOutputFile.close();
		TimeToCellDivisionFile.close();
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MitosisInitialisation::finish(){	
}


void MitosisInitialisation::update(CC3DXMLElement *_xmlData, bool _fullInitFlag){
     if(_xmlData->findElement("SCCtype")){SCCtype=_xmlData->getFirstElement("SCCtype")->getDouble();}
     else{SCCtype=1;cerr<<"SCCtype not specified, default value set to: "<<SCCtype<<endl;}
     if(_xmlData->findElement("VCAFtype")){VCAFtype=_xmlData->getFirstElement("VCAFtype")->getDouble();}
     else{VCAFtype=2;cerr<<"VCAFtype not specified, default value set to: "<<VCAFtype<<endl;}
     if(_xmlData->findElement("SCCVolume")){SCCVolume=_xmlData->getFirstElement("SCCVolume")->getDouble();}
     else{SCCVolume=550;cerr<<"SCCVolume not specified, default value set to: "<<SCCVolume<<endl;}
     if(_xmlData->findElement("VCAFVolume")){VCAFVolume=_xmlData->getFirstElement("VCAFVolume")->getDouble();}
     else{VCAFVolume=800;cerr<<"VCAFVolume not specified, default value set to: "<<VCAFVolume<<endl;}  
     if(_xmlData->findElement("SCCLambdaVolume")){SCCLambdaVolume=_xmlData->getFirstElement("SCCLambdaVolume")->getDouble();}
     else{SCCLambdaVolume=1;cerr<<"SCCLambdaVolume not specified, default value set to: "<<SCCLambdaVolume<<endl;}
     if(_xmlData->findElement("VCAFLambdaVolume")){VCAFLambdaVolume=_xmlData->getFirstElement("VCAFLambdaVolume")->getDouble();}
     else{VCAFLambdaVolume=1;cerr<<"VCAFLambdaVolume not specified, default value set to: "<<VCAFLambdaVolume<<endl;}
     if(_xmlData->findElement("SCCSurface")){SCCSurface=_xmlData->getFirstElement("SCCSurface")->getDouble();}
     else{SCCSurface=550;cerr<<"SCCSurface not specified, default value set to: "<<SCCSurface<<endl;}
     if(_xmlData->findElement("VCAFSurface")){VCAFSurface=_xmlData->getFirstElement("VCAFSurface")->getDouble();}
     else{VCAFSurface=700;cerr<<"VCAFSurface not specified, default value set to: "<<VCAFSurface<<endl;}  
     if(_xmlData->findElement("SCCLambdaSurface")){SCCLambdaSurface=_xmlData->getFirstElement("SCCLambdaSurface")->getDouble();}
     else{SCCLambdaSurface=1;cerr<<"SCCLambdaSurface not specified, default value set to: "<<SCCLambdaSurface<<endl;}
     if(_xmlData->findElement("VCAFLambdaSurface")){VCAFLambdaSurface=_xmlData->getFirstElement("VCAFLambdaSurface")->getDouble();}
     else{VCAFLambdaSurface=1;cerr<<"VCAFLambdaSurface not specified, default value set to: "<<VCAFLambdaSurface<<endl;}
     if(_xmlData->findElement("InitialTimeToMitosisSCC")){InitialTimeToMitosisSCC=_xmlData->getFirstElement("InitialTimeToMitosisSCC")->getDouble();}
     else{InitialTimeToMitosisSCC=5760;cerr<<"InitialTimeToMitosisSCC not specified, default value set to: "<<InitialTimeToMitosisSCC<<"steps"<<endl;}
     if(_xmlData->findElement("MeanTimeToMitosisCheckSCC")){MeanTimeToMitosisCheckSCC=_xmlData->getFirstElement("MeanTimeToMitosisCheckSCC")->getDouble();}
     else{MeanTimeToMitosisCheckSCC=5760;cerr<<"MeanTimeToMitosisCheckSCC not specified, default value set to: "<<MeanTimeToMitosisCheckSCC<<"steps"<<endl;}
     if(_xmlData->findElement("VolumeMultiplierForGrowth")){ VolumeMultiplierForGrowth=_xmlData->getFirstElement("VolumeMultiplierForGrowth")->getDouble();}
     else{ VolumeMultiplierForGrowth=2;cerr<<"VolumeMultiplierForGrowth not specified, default value set to: "<<VolumeMultiplierForGrowth<<"steps"<<endl;}
     if(_xmlData->findElement("DivisionThresholdMultiplier")){ DivisionThresholdMultiplier=_xmlData->getFirstElement("DivisionThresholdMultiplier")->getDouble();}
     else{ DivisionThresholdMultiplier=0.9;cerr<<"DivisionThresholdMultiplier not specified, default value set to: "<<DivisionThresholdMultiplier<<"steps"<<endl;}
     if(_xmlData->findElement("InitialCellOutputFileName")){InitialCellOutputFileName=_xmlData->getFirstElement("InitialCellOutputFileName")->getText();}
     if(_xmlData->findElement("TimeToCellDivisionFileName")){TimeToCellDivisionFileName=_xmlData->getFirstElement("TimeToCellDivisionFileName")->getText();}
     if(_xmlData->findElement("TimeToWriteCells")){TimeToWriteCells=_xmlData->getFirstElement("TimeToWriteCells")->getDouble();}
     else{TimeToWriteCells=5800;cerr<<"TimeToWriteCells not specified, default value set to: "<<TimeToWriteCells<<"steps"<<endl;}
     if(_xmlData->findElement("DivisionTimeDistribution")){DivisionTimeDistribution=_xmlData->getFirstElement("DivisionTimeDistribution")->getText();}
     else{std::string str ("Exponential");DivisionTimeDistribution=str;cerr<<"DivisionTimeDistribution not specified, default distribution set to: "<<DivisionTimeDistribution<<endl;}
     if(_xmlData->findElement("ErlangShapeParameter")){ErlangShapeParameter=_xmlData->getFirstElement("ErlangShapeParameter")->getDouble();}
     else{ErlangShapeParameter=1;cerr<<"ErlangShapeParameter not specified, default value set to: "<<ErlangShapeParameter<<endl;}
     if(_xmlData->findElement("SurfacevolumeRatio")){SurfacevolumeRatio=_xmlData->getFirstElement("SurfacevolumeRatio")->getText();}
     else{std::string str ("Sphere");DivisionTimeDistribution=str;cerr<<"DivisionTimeDistribution not specified, default distribution set to: "<<SurfacevolumeRatio<<endl;}  
}

std::string MitosisInitialisation::toString(){
   return "MitosisInitialisation";
}


std::string MitosisInitialisation::steerableName(){
   return toString();
}


