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
#include <vector> 
#include <algorithm> 
#include <numeric>
using namespace CompuCell3D;
#include <iostream>
using namespace std;
#include <fstream>
#include <iomanip>
#include <sstream>
//This is a random number generator written by Xioafen Li.
#include "RandomGenerator.h"
#include <CompuCell3D/Potts3D/CellInventory.h>


#include "MitosisPostInitialise.h"
MitosisPostInitialise::MitosisPostInitialise() :  cellFieldG(0),sim(0),potts(0) {}
MitosisPostInitialise::~MitosisPostInitialise() {
}


void MitosisPostInitialise::init(Simulator *simulator, CC3DXMLElement *_xmlData) {
	potts = simulator->getPotts();
	sim=simulator;
	cellFieldG = (WatchableField3D<CellG *> *)potts->getCellFieldG();
	simulator->registerSteerableObject(this);
	update(_xmlData);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MitosisPostInitialise::extraInit(Simulator *simulator){

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MitosisPostInitialise::start(){
	cout<<"within the start function of MitosisPostInitialise"<<endl;
	const char *outputname;
	outputname=TimeOfCellDivisionFileName.c_str();
	TimeOfCellDivisionFile.open(outputname,ofstream::out);
	if (!TimeOfCellDivisionFile.is_open()){cerr<<"could not open "<<TimeOfCellDivisionFile<<"!!"<<endl;}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MitosisPostInitialise::step(const unsigned int currentStep){
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
		const char *inputname;
		inputname=MitosisStageInputFileName.c_str();
		MitosisStageInputFile.open(inputname,ifstream::in);
		if(!MitosisStageInputFile.is_open()){cerr<<"could not open "<<MitosisStageInputFile<<" !!"<<endl;}
		int counter=-1;
		for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
			cell=cellInventoryPtr->getCell(cInvItr);
			if (cell->type == SCCtype){
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
				cell->OriginalTargetVolume = cell->targetVolume;
				cell->OriginalTargetSurface = cell->targetSurface;
				MitosisStageInputFile>>cell->id;
				MitosisStageInputFile>>cell->mitosistimer;
				MitosisStageInputFile>>cell->TimeToMitosisCheck;
				MitosisStageInputFile>>cell->MaxCellVolume;
				MitosisStageInputFile>>cell->targetVolume;
				MitosisStageInputFile>>cell->targetSurface;
			}
			else if(cell->type == VCAFtype){
				cell->targetVolume=VCAFVolume;
				cell->lambdaVolume=VCAFLambdaVolume;
				cell->targetSurface=VCAFSurface;
				cell->lambdaSurface=VCAFLambdaSurface;
			}			
		}
		MitosisStageInputFile.close();
	}
	else{
		for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
			cell=cellInventoryPtr->getCell(cInvItr);
			if (cell->type == SCCtype){
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
					TimeOfCellDivisionFile<<setw(5)<<cell->id<<setw(15)<<cell->mitosistimer<<setw(15)<<cell->TimeToMitosisCheck<<setw(10)<<cell->volume<<setw(10)<<cell->targetVolume<<setw(10)<<endl;
					cell->DivideThisCellInPython=true;
					cell->targetVolume=cell->OriginalTargetVolume;
					cell->targetSurface=cell->OriginalTargetSurface;
					cell->mitosistimer=0;
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


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MitosisPostInitialise::finish(){	
}

void MitosisPostInitialise::update(CC3DXMLElement *_xmlData, bool _fullInitFlag){
     if(_xmlData->findElement("SCCtype")){SCCtype=_xmlData->getFirstElement("SCCtype")->getDouble();}
     else{SCCtype=1;cerr<<"SCCtype not specified, defaul t value set to: "<<SCCtype<<endl;}
     if(_xmlData->findElement("VCAFtype")){VCAFtype=_xmlData->getFirstElement("VCAFtype")->getDouble();}
     else{VCAFtype=2;cerr<<"VCAFtype not specified, defaul t value set to: "<<VCAFtype<<endl;}
     if(_xmlData->findElement("SCCVolume")){SCCVolume=_xmlData->getFirstElement("SCCVolume")->getDouble();}
     else{SCCVolume=550;cerr<<"SCCVolume not specified, defaul t value set to: "<<SCCVolume<<endl;}
     if(_xmlData->findElement("VCAFVolume")){VCAFVolume=_xmlData->getFirstElement("VCAFVolume")->getDouble();}
     else{VCAFVolume=800;cerr<<"VCAFVolume not specified, defaul t value set to: "<<VCAFVolume<<endl;}  
     if(_xmlData->findElement("SCCLambdaVolume")){SCCLambdaVolume=_xmlData->getFirstElement("SCCLambdaVolume")->getDouble();}
     else{SCCLambdaVolume=1;cerr<<"SCCLambdaVolume not specified, defaul t value set to: "<<SCCLambdaVolume<<endl;}
     if(_xmlData->findElement("VCAFLambdaVolume")){VCAFLambdaVolume=_xmlData->getFirstElement("VCAFLambdaVolume")->getDouble();}
     else{VCAFLambdaVolume=1;cerr<<"VCAFLambdaVolume not specified, defaul t value set to: "<<VCAFLambdaVolume<<endl;}
     if(_xmlData->findElement("SCCSurface")){SCCSurface=_xmlData->getFirstElement("SCCSurface")->getDouble();}
     else{SCCSurface=550;cerr<<"SCCSurface not specified, default value set to: "<<SCCSurface<<endl;}
     if(_xmlData->findElement("VCAFSurface")){VCAFSurface=_xmlData->getFirstElement("VCAFSurface")->getDouble();}
     else{VCAFSurface=700;cerr<<"VCAFSurface not specified, default value set to: "<<VCAFSurface<<endl;}  
     if(_xmlData->findElement("SCCLambdaSurface")){SCCLambdaSurface=_xmlData->getFirstElement("SCCLambdaSurface")->getDouble();}
     else{SCCLambdaSurface=1;cerr<<"SCCLambdaSurface not specified, default value set to: "<<SCCLambdaSurface<<endl;}
     if(_xmlData->findElement("VCAFLambdaSurface")){VCAFLambdaSurface=_xmlData->getFirstElement("VCAFLambdaSurface")->getDouble();}
     else{VCAFLambdaSurface=1;cerr<<"VCAFLambdaSurface not specified, default value set to: "<<VCAFLambdaSurface<<endl;}
     if(_xmlData->findElement("MeanTimeToMitosisCheckSCC")){MeanTimeToMitosisCheckSCC=_xmlData->getFirstElement("MeanTimeToMitosisCheckSCC")->getDouble();}
     else{MeanTimeToMitosisCheckSCC=2880;cerr<<"MeanTimeToMitosisCheckSCC not specified, default value set to: "<<MeanTimeToMitosisCheckSCC<<"steps"<<endl;}
     if(_xmlData->findElement("VolumeMultiplierForGrowth")){ VolumeMultiplierForGrowth=_xmlData->getFirstElement("VolumeMultiplierForGrowth")->getDouble();}
     else{ VolumeMultiplierForGrowth=2;cerr<<"VolumeMultiplierForGrowth not specified, default value set to: "<<VolumeMultiplierForGrowth<<"steps"<<endl;}
     if(_xmlData->findElement("DivisionThresholdMultiplier")){ DivisionThresholdMultiplier=_xmlData->getFirstElement("DivisionThresholdMultiplier")->getDouble();}
     else{ DivisionThresholdMultiplier=0.8;cerr<<"DivisionThresholdMultiplier not specified, default value set to: "<<DivisionThresholdMultiplier<<"steps"<<endl;}
     if(_xmlData->findElement("MitosisStageInputFileName")){MitosisStageInputFileName=_xmlData->getFirstElement("MitosisStageInputFileName")->getText();}
     if(_xmlData->findElement("TimeOfCellDivisionFileName")){TimeOfCellDivisionFileName=_xmlData->getFirstElement("TimeOfCellDivisionFileName")->getText();}
     if(_xmlData->findElement("DivisionTimeDistribution")){DivisionTimeDistribution=_xmlData->getFirstElement("DivisionTimeDistribution")->getText();}
     else{std::string str ("Exponential");DivisionTimeDistribution=str;cerr<<"DivisionTimeDistribution not specified, default distribution set to: "<<DivisionTimeDistribution<<endl;}
     if(_xmlData->findElement("ErlangShapeParameter")){ErlangShapeParameter=_xmlData->getFirstElement("ErlangShapeParameter")->getDouble();}
     else{ErlangShapeParameter=1;cerr<<"ErlangShapeParameter not specified, default value set to: "<<ErlangShapeParameter<<endl;}
     if(_xmlData->findElement("SurfacevolumeRatio")){SurfacevolumeRatio=_xmlData->getFirstElement("SurfacevolumeRatio")->getText();}
     else{std::string str ("Sphere");DivisionTimeDistribution=str;cerr<<"DivisionTimeDistribution not specified, default distribution set to: "<<SurfacevolumeRatio<<endl;}       
}

std::string MitosisPostInitialise::toString(){
   return "MitosisPostInitialise";
}


std::string MitosisPostInitialise::steerableName(){
   return toString();
}


