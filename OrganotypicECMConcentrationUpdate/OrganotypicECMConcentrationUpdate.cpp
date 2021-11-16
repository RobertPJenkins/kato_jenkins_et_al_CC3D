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
#include <fstream>
#include <iomanip>
#include <sstream>


#include "OrganotypicECMConcentrationUpdate.h"
#include "CompuCell3D/plugins/PixelTracker/PixelTrackerPlugin.h"
#include "CompuCell3D/plugins/PixelTracker/PixelTracker.h"

#include "CompuCell3D/plugins/BoundaryPixelTracker/BoundaryPixelTrackerPlugin.h"
#include "CompuCell3D/plugins/BoundaryPixelTracker/BoundaryPixelTracker.h"

OrganotypicECMConcentrationUpdate::OrganotypicECMConcentrationUpdate() :  cellFieldG(0),sim(0),potts(0) {}

OrganotypicECMConcentrationUpdate::~OrganotypicECMConcentrationUpdate() {
}


void OrganotypicECMConcentrationUpdate::init(Simulator *simulator, CC3DXMLElement *_xmlData) {

	potts = simulator->getPotts();
	sim=simulator;
	cellFieldG = (WatchableField3D<CellG *> *)potts->getCellFieldG();
	simulator->registerSteerableObject(this);
	update(_xmlData);
	bool pluginAlreadyRegisteredFlag;
	Plugin *plugin=Simulator::pluginManager.get("VolumeTracker",&pluginAlreadyRegisteredFlag); 
	cerr<<"GOT HERE BEFORE CALLING INIT"<<endl;
	if(!pluginAlreadyRegisteredFlag)
		plugin->init(simulator);
	Plugin *pluginCOM=Simulator::pluginManager.get("CenterOfMass",&pluginAlreadyRegisteredFlag); 
	cerr<<"GOT HERE BEFORE CALLING INIT"<<endl;
	if(!pluginAlreadyRegisteredFlag)
		pluginCOM->init(simulator);
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

void OrganotypicECMConcentrationUpdate::extraInit(Simulator *simulator){

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void OrganotypicECMConcentrationUpdate::start(){
	cout<<"within the start function of ECM"<<endl;
	const int x=dimx;
	const int y=dimy;
	const int z=dimz;
	ECMconcfield = new float**[x];
	for(int i = 0; i < x; i++){
		ECMconcfield[i]= new float*[y];
		for (int j = 0; j< y; j++){
			ECMconcfield[i][j] = new float[z];
		}
	}
	pixelupdatestep = new int**[x];
	for(int i = 0; i < x; i++){
		pixelupdatestep[i]= new int*[y];
		for (int j = 0; j< y; j++){
			pixelupdatestep[i][j] = new int[z];
		}
	}
	for (int i = 0; i< x; i++){
	for (int j = 0; j< y; j++){
	for (int k = 0; k< z; k++){
		ECMconcfield[i][j][k]=1.0;
		pixelupdatestep[i][j][k]=0;
	}}}
	
	const char *inputname;
	inputname=ECMInputFileName.c_str();
	ECMInputFile.open(inputname,ifstream::in);
	if(!ECMInputFile.is_open()){cerr<<"could not open "<<ECMInputFileName<<" !!"<<endl;}
	int xcoord, ycoord, zcoord;
	float concentration;
	while( !ECMInputFile.eof() ) {
		ECMInputFile>>xcoord;
		ECMInputFile>>ycoord;
		ECMInputFile>>zcoord;
		ECMInputFile>>concentration;
		ECMconcfield[xcoord][ycoord][zcoord]=concentration;
		pixelupdatestep[xcoord][ycoord][zcoord]=0;
	}
	ECMInputFile.close();
	
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void OrganotypicECMConcentrationUpdate::step(const unsigned int currentStep){
	const int ECMSavePeriod=ECMOutputSavePeriod;
	CellInventory *cellInventoryPtr=& potts->getCellInventory();
	CellInventory::cellInventoryIterator cInvItr;
	CellG *cell;
	if(currentStep==0){
		for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
			cell=cellInventoryPtr->getCell(cInvItr);
			cell->ECMfieldbyCell=ECMconcfield;
			if((int)cell->type==VCAFtype){
				cell->lambdaECMPenetration=VCAFLambdaECMpenetration;
				cell->kdeg=ecmDegradationRate_max_vcaf;
				cell->wpush=ecmPushingRate_max_vcaf;
			}
			else if((int)cell->type==SCCtype){
				cell->lambdaECMPenetration=SCCLambdaECMpenetration;
				cell->kdeg=ecmDegradationRate_max_scc;
				cell->wpush=ecmPushingRate_max_scc;
			}
			set<PixelTrackerData> cellPixels=pixelTrackerAccessorPtr->get(cell->extraAttribPtr)->pixelSet;
			for(set<PixelTrackerData>::iterator sitr=cellPixels.begin() ; sitr != cellPixels.end() ;++sitr){
				int xx=sitr->pixel.x;int yy=sitr->pixel.y;int zz=sitr->pixel.z;
				if(zz<=SpaceBoundary){
					ECMconcfield[xx][yy][zz]=0.0;
				}
				else{
					ECMconcfield[xx][yy][zz]=CellECMConc;
				}
			}
			cell->steppablesinitialised=1;
		}
	}
	else{
		for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
			cell=cellInventoryPtr->getCell(cInvItr);
			if (cell->just_divided == true ){
				cell->ECMfieldbyCell=ECMconcfield;
				if((int)cell->type==VCAFtype){cell->lambdaECMPenetration=VCAFLambdaECMpenetration;}
				else if((int)cell->type==SCCtype){cell->lambdaECMPenetration=SCCLambdaECMpenetration;}
				cell->steppablesinitialised=1;
				cell->just_divided=false;
			}
		}
		vector <int> xcoord,ycoord,zcoord;
		vector <float> kdeglist,wpushlist;
		for(cInvItr=cellInventoryPtr->cellInventoryBegin() ; cInvItr !=cellInventoryPtr->cellInventoryEnd() ;++cInvItr ){
			cell=cellInventoryPtr->getCell(cInvItr);
			float weightlambda;
			float currkdeg;
			float currwpush;
			if((int)cell->type==VCAFtype){
				weightlambda=(cell->lambdataxis-lambdamin)/(lambdamax-lambdamin);	
			}
			else if((int)cell->type==SCCtype){
				weightlambda=(lambdaSCC_VCAFequiv-lambdamin)/(lambdamax-lambdamin);	
			}
			currkdeg=cell->kdeg*weightlambda;
			currwpush=cell->wpush*weightlambda;
			set<BoundaryPixelTrackerData> cellBoundaryPixels=boundaryPixelTrackerAccessorPtr->get(cell->extraAttribPtr)->pixelSet;
			for(set<BoundaryPixelTrackerData>::iterator sitr=cellBoundaryPixels.begin() ; sitr != cellBoundaryPixels.end() ;++sitr){
				int xx=sitr->pixel.x;int yy=sitr->pixel.y;int zz=sitr->pixel.z;
				for(int i=-1;i<2;i++){
					for(int j=-1;j<2;j++){
						for(int k=-1;k<2;k++){
							int currz=zz+k;
							if(currz>0 && currz<dimz){
								int currx=xx+i;int curry=yy+j;
								if(currx>=dimx){currx=currx-dimx;}else if(currx<0){currx=currx+dimx;}
								if(curry>=dimy){curry=curry-dimy;}else if(curry<0){curry=curry+dimy;}
								if(pixelupdatestep[currx][curry][currz]!=currentStep){
									xcoord.push_back(currx);ycoord.push_back(curry);zcoord.push_back(currz);
									kdeglist.push_back(currkdeg);wpushlist.push_back(currwpush);
									pixelupdatestep[currx][curry][currz]=currentStep;
								}
							}
						}
					}	
				}	
			}
			set<PixelTrackerData> cellPixels=pixelTrackerAccessorPtr->get(cell->extraAttribPtr)->pixelSet;
			for(set<PixelTrackerData>::iterator sitr=cellPixels.begin() ; sitr != cellPixels.end() ;++sitr){
				int xx=sitr->pixel.x;int yy=sitr->pixel.y;int zz=sitr->pixel.z;
				if(pixelupdatestep[xx][yy][zz]!=currentStep){
					xcoord.push_back(xx);ycoord.push_back(yy);zcoord.push_back(zz);
					kdeglist.push_back(currkdeg);wpushlist.push_back(currwpush);
					pixelupdatestep[xx][yy][zz]=currentStep;
				}
			}
		}
		const int N=xcoord.size();
		int randomarray[N];
		for (int i=0;i<N;i++){randomarray[i]=i;}
		for (int i=0; i<(N-1); i++){
 			int r=i + (rand() % (N-i));
			int temp = randomarray[i]; 
			randomarray[i] = randomarray[r]; 
			randomarray[r] = temp;
		}
		vector<int>::iterator intItr;
		vector<float>::iterator floatItr;
		for (int m=0;m<N;m++){
			intItr=xcoord.begin();intItr=intItr+randomarray[m];int xx=(*intItr);
			intItr=ycoord.begin();intItr=intItr+randomarray[m];int yy=(*intItr);
			intItr=zcoord.begin();intItr=intItr+randomarray[m];int zz=(*intItr);
			floatItr=kdeglist.begin();floatItr=floatItr+randomarray[m];float currkdeg=(*floatItr);
			floatItr=wpushlist.begin();floatItr=floatItr+randomarray[m];float currwpush=(*floatItr);
			if(ECMconcfield[xx][yy][zz]<0){
				cout<<"current concentration read: "<<ECMconcfield[xx][yy][zz]<<endl;
			}
			ECMconcfield[xx][yy][zz]=ECMconcfield[xx][yy][zz]-currkdeg;if(ECMconcfield[xx][yy][zz]<0.0){ECMconcfield[xx][yy][zz]=0.0;}
			float pushedvalue=ECMconcfield[xx][yy][zz]*currwpush/26.0;
			for(int i=-1;i<2;i++){
			for(int j=-1;j<2;j++){
			for(int k=-1;k<2;k++){
				if (i!=0 || j!=0 || k!=0){
					int currz=zz+k;
					if(currz>0 && currz<dimz){
						int currx=xx+i;int curry=yy+j;
						if(currx>=dimx){currx=currx-dimx;}else if(currx<0){currx=currx+dimx;}
						if(curry>=dimy){curry=curry-dimy;}else if(curry<0){curry=curry+dimy;}
						ECMconcfield[currx][curry][currz]=ECMconcfield[currx][curry][currz]+pushedvalue;
					}
				}
			}}}
			ECMconcfield[xx][yy][zz]=ECMconcfield[xx][yy][zz]*(1-currwpush);
		}
	}
	
	if(currentStep%ECMSavePeriod==0){
		int timestepwrite=currentStep;
		const char *outputname;
		string paddedtsstring, ecmoutputname_var,filetype;
		stringstream paddedts;
		paddedts <<setw(5) << setfill('0') <<timestepwrite;
		paddedtsstring=paddedts.str();
		filetype=".txt";
		ecmoutputname_var=ECMOutputFileName+paddedtsstring+filetype;
		cout<<ecmoutputname_var<<endl;
		outputname=ecmoutputname_var.c_str();
		ECMOutputFile.open(outputname,ofstream::trunc);
		if (!ECMOutputFile.is_open()){cerr<<"could not open "<<ecmoutputname_var<<"!!"<<endl;}
		ECMOutput();
		ECMOutputFile.close();
	}
}

void OrganotypicECMConcentrationUpdate::ECMOutput(){
	const int x=dimx;
	const int y=dimy;
	const int z=dimz;
	for (int i = 0; i< x; i++){
		for (int j = 0; j< y; j++){
			for (int k = 0; k< z; k++){
				if (ECMconcfield[i][j][k]!=1.0){
					if(ECMconcfield[i][j][k]!=0.0){
						ECMOutputFile<<setw(5) <<i<<" "<<setw(5) <<j<<" "<<setw(5) <<k<<"	" <<ECMconcfield[i][j][k]<<endl;
					}
					else if(ECMconcfield[i][j][k]==0.0){	
						ECMOutputFile<<setw(5) <<i<<" "<<setw(5) <<j<<" "<<setw(5) <<k<<"	" <<showpoint<<setprecision(2)<<0.0<<endl;
					}
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void OrganotypicECMConcentrationUpdate::finish(){	
}


void OrganotypicECMConcentrationUpdate::update(CC3DXMLElement *_xmlData, bool _fullInitFlag){

     if(_xmlData->findElement("VCAFtype"))
     VCAFtype=_xmlData->getFirstElement("VCAFtype")->getDouble();
     if(_xmlData->findElement("SCCtype")){SCCtype=_xmlData->getFirstElement("SCCtype")->getDouble();}
     else{SCCtype=1;cerr<<"SCCtype not specified, defaul t value set to: "<<SCCtype<<endl;}
     if(_xmlData->findElement("dimx")){dimx=_xmlData->getFirstElement("dimx")->getDouble();}
     if(_xmlData->findElement("dimy")){dimy=_xmlData->getFirstElement("dimy")->getDouble();}
     if(_xmlData->findElement("dimz")){dimz=_xmlData->getFirstElement("dimz")->getDouble();}
     if(_xmlData->findElement("SCCLambdaTaxis_VCAFequiv")){lambdaSCC_VCAFequiv=_xmlData->getFirstElement("SCCLambdaTaxis_VCAFequiv")->getDouble();}
     if(_xmlData->findElement("VCAFLambdaTaxis_min")){lambdamin=_xmlData->getFirstElement("VCAFLambdaTaxis_min")->getDouble();}
     if(_xmlData->findElement("VCAFLambdaTaxis_max")){lambdamax=_xmlData->getFirstElement("VCAFLambdaTaxis_max")->getDouble();}
     if(_xmlData->findElement("ECMDegradationRate_max_vcaf")){ecmDegradationRate_max_vcaf=_xmlData->getFirstElement("ECMDegradationRate_max_vcaf")->getDouble();}
     if(_xmlData->findElement("ECMPushingRate_max_vcaf")){ecmPushingRate_max_vcaf=_xmlData->getFirstElement("ECMPushingRate_max_vcaf")->getDouble();}
     if(_xmlData->findElement("ECMDegradationRate_max_scc")){ecmDegradationRate_max_scc=_xmlData->getFirstElement("ECMDegradationRate_max_scc")->getDouble();}
     if(_xmlData->findElement("ECMPushingRate_max_scc")){ecmPushingRate_max_scc=_xmlData->getFirstElement("ECMPushingRate_max_scc")->getDouble();}
     if(_xmlData->findElement("VCAFLambdaECMpenetration")){VCAFLambdaECMpenetration=_xmlData->getFirstElement("VCAFLambdaECMpenetration")->getDouble();}
     if(_xmlData->findElement("SCCLambdaECMpenetration")){SCCLambdaECMpenetration=_xmlData->getFirstElement("SCCLambdaECMpenetration")->getDouble();}
     if(_xmlData->findElement("ECMOutputFileName")){ECMOutputFileName=_xmlData->getFirstElement("ECMOutputFileName")->getText();}
     if(_xmlData->findElement("ECMOutputSavePeriod")){ECMOutputSavePeriod=_xmlData->getFirstElement("ECMOutputSavePeriod")->getDouble();}
     if(_xmlData->findElement("ECMInputFileName")){ECMInputFileName=_xmlData->getFirstElement("ECMInputFileName")->getText();}
     if(_xmlData->findElement("ECMSpaceBoundary")){SpaceBoundary=_xmlData->getFirstElement("ECMSpaceBoundary")->getDouble();}
     else{SpaceBoundary=10000;cerr<<"ECMSpaceBoundary not specified, default value set to: "<<SpaceBoundary<<endl;}
     if(_xmlData->findElement("CellECMConcentration")){CellECMConc=_xmlData->getFirstElement("CellECMConcentration")->getDouble();}
     else{CellECMConc=0.0;cerr<<"CellECMConcentration not specified, default value set to: "<<CellECMConc<<endl;}
}

std::string OrganotypicECMConcentrationUpdate::toString(){
   return "OrganotypicECMConcentrationUpdate";
}


std::string OrganotypicECMConcentrationUpdate::steerableName(){
   return toString();
}


