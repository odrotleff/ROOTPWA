///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////

// This Class' Header ------------------
#include "pwaPlotter.h"

// C/C++ Headers ----------------------
#include <iostream>
#include <limits>

// Collaborating Class Headers --------
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "fitResult.h"

// Class Member definitions -----------

using namespace std;
using namespace rpwa;

TH2D* drawDensity(TGraphErrors* g, TH2D* h, double weight){
  
  unsigned int ybins=h->GetNbinsY();
  unsigned int xpoints=g->GetN();
  TAxis* ax=h->GetYaxis();
  double ymin=ax->GetXmin();
  double ymax=ax->GetXmax();

  TF1 gaus("gaus","gausn(0)",ymin,ymax);
  TString tit=g->GetTitle();
  
  for(unsigned int ip=0; ip<xpoints; ++ip){
    double x=g->GetX()[ip];
    double err=g->GetEY()[ip];
    double val=g->GetY()[ip];
  
 //    if(tit.Contains("1-1++0+sigma_01_a11269") && val<100 && x>1.6 && x<2.0){
//       //cerr << "x="<<x
//       //	   << "  err="<< err 
//       //	   << "  val="<< val << endl;
//       //return h;
//     }
      
    gaus.SetParameters(1,val,err);
    for(unsigned int ibin=1; ibin<ybins;++ibin){
      double y=ax->GetBinCenter(ibin);
      //if(fabs(y-val)<5.*err){
	double w=gaus.Eval(ax->GetBinCenter(ibin))*weight;
	if(w==w)h->Fill(x,y,w);
	//}
    }
  }
  return h;
}



string
getIGJPCMEps(const std::string& wavename){
  return wavename.substr(0,7);
}



ClassImp(pwaPlotter);

pwaPlotter::pwaPlotter()
  : mMinEvidence(0)
{
  mLogLikelihood= new TMultiGraph();
  mLogLikelihood->SetTitle("LogLikelihood");
  mLogLikelihood->SetName("LogLikelihood");
  mLogLikelihoodPerEvent=new TMultiGraph();
  mLogLikelihoodPerEvent->SetTitle("LogLikelihood/Event");
  mLogLikelihoodPerEvent->SetName("LogLikelihoodPerEvent");
  mEvidence= new TMultiGraph();
  mEvidence->SetTitle("Evidence");
  mEvidence->SetName("Evidence");
  mEvidencePerEvent=new TMultiGraph();
  mEvidencePerEvent->SetTitle("Evidence/Event");
  mEvidencePerEvent->SetName("EvidencePerEvent");
  
}


pwaPlotter::~pwaPlotter(){
  mWavenames.clear();
  
}

void 
pwaPlotter::addFit(const std::string& filename,
		   const std::string& title,
		   const unsigned int colour,
		   const std::string& treename,
		   const std::string& branchname){
  
  // Open and test file and tree
  TFile* infile = TFile::Open(filename.c_str(),"READ");
  if(infile==NULL || infile->IsZombie()){
    cerr << "Input file "<<filename<<" is not a valid file!" << endl;
    return;
  }
  TTree* intree=(TTree*)infile->FindObjectAny(treename.c_str());
  if(intree==NULL || intree->IsZombie()){
    cerr << "Tree "<<treename<<" not found in file "<<filename<< endl;
    return;
  }
  fitResult* result=0;
  if(intree->FindBranch(branchname.c_str())==NULL){
    cerr << "Invalid branch "<<treename<<"."<<branchname<<" in file "
	 <<filename<<endl;
    return;
  }
  
  unsigned int ifit=mResultMetaInfo.size();
  cerr << "Adding file "<< filename << endl;
  
  intree->SetBranchAddress(branchname.c_str(),&result);
  unsigned int nbins=intree->GetEntries();
  // extract info for this fit
  // loop through bins
  // -> getRange in Mass bins
  // -> collect all used waves
  // -> integrate loglikelihood and evidence
  double mass_min=1E6;
  double mass_max=0;
  double logli=0;
  double logliperevt=0;
  double evi=0;
  double eviperevt=0;
  
  set<string> wavesinthisfit;
  
  for(unsigned int i=0;i<nbins;++i){
    intree->GetEntry(i);
    
    double massBinCenter=result->massBinCenter()*0.001;
    if(massBinCenter>mass_max)mass_max=massBinCenter;
    if(massBinCenter<mass_min)mass_min=massBinCenter;
    
    registerWave(".*"); // Total intensity
    wavesinthisfit.insert(".*");
    registerWave("^.....0"); // Total M=0
    wavesinthisfit.insert("^.....0");
    registerWave("^.....1"); // Total M=1
    wavesinthisfit.insert("^.....1");

    registerWave         ("^......\\+"); // Total Eps=+
    wavesinthisfit.insert("^......\\+");
    registerWave("^......-"); // Total Eps=-
    wavesinthisfit.insert("^......-");

    // check fitResult for used waves
    // if not already registered -> register wave (will create TMultiGraph)
    const vector<string>& waveNames=result->waveNames();
    unsigned int nwaves=waveNames.size();
    for(unsigned int iw=0;iw<nwaves;++iw){
      registerWave(waveNames[iw]);
      wavesinthisfit.insert(waveNames[iw]);
      // spin totals...
      registerWave(getIGJPCMEps(waveNames[iw]));
      wavesinthisfit.insert(getIGJPCMEps(waveNames[iw]));
    }
    
    // get loglikelihoods
    logli+=result->logLikelihood();
    evi+=result->evidence();
    logliperevt+=result->logLikelihood()/result->nmbEvents();
    eviperevt+=result->evidence()/result->nmbEvents();
  }
  double binwidth=(mass_max-mass_min)/(double)(nbins-1);
  cerr << "Number of bins: " << nbins 
       << "   Width: " << binwidth << endl;
  
  
  // create intensity plots ----------------------------------------------
  // We have registered all graphs in the step before...
  // This has to be done in a separate step! Try not to merge the following
  // with the loop above! You will loose generality!!! You have been warned!
  
  //cout << "creating graphs" << endl;
 

  // create graphs for this fit
  set<string>::iterator it=wavesinthisfit.begin();
  while(it!=wavesinthisfit.end()){
    TPwaFitGraphErrors* g = new TPwaFitGraphErrors(nbins,ifit);
    stringstream graphName;
    if(*it==".*") graphName << "g" << title << "_total";
    else if(*it=="^.....0") graphName << "g" << title << "_M0";
    else if(*it=="^.....1") graphName << "g" << title << "_M1";
    else if(*it=="^......\\+") graphName << "g" << title << "_E+";
    else if(*it=="^......-") graphName << "g" << title << "_E-";
    else graphName << "g" << title << "_" << *it;

    //cout << "creating graph   " << graphName.str() << endl;

    g->SetName (graphName.str().c_str());
    g->SetTitle(graphName.str().c_str());
    g->SetMarkerStyle(21);
    g->SetMarkerSize(0.5);
    g->SetMarkerColor(colour);
    g->SetLineColor(colour);
    mIntensities[*it]->Add(g,"p");
    ++it;
  }

  //cout << "building Likelihood graphs" << endl;

  // evidence and likelihood
  TGraph* gLikeli=new TGraph(nbins);
  stringstream graphName;
  graphName << "g" << title << "_LogLikelihood";
  gLikeli->SetName (graphName.str().c_str());
  gLikeli->SetTitle(graphName.str().c_str());
  gLikeli->SetMarkerStyle(21);
  gLikeli->SetMarkerSize(0.5);
  gLikeli->SetMarkerColor(colour);
  gLikeli->SetLineColor(colour);
  mLogLikelihood->Add(gLikeli,"p");
  TGraph* gLikeliPE=new TGraph(nbins);
  graphName.clear();
  graphName << "g" << title << "_LogLikelihoodPerEvent";
  gLikeliPE->SetName (graphName.str().c_str());
  gLikeliPE->SetTitle(graphName.str().c_str());
  gLikeliPE->SetMarkerStyle(21);
  gLikeliPE->SetMarkerSize(0.5);
  gLikeliPE->SetMarkerColor(colour);
  gLikeliPE->SetLineColor(colour);
  mLogLikelihoodPerEvent->Add(gLikeliPE,"p");
  TGraph* gEvidence=new TGraph(nbins);
  graphName.clear();
  graphName << "g" << title << "_Evidence";
  gEvidence->SetName (graphName.str().c_str());
  gEvidence->SetTitle(graphName.str().c_str());
  gEvidence->SetMarkerStyle(21);
  gEvidence->SetMarkerSize(0.5);
  gEvidence->SetMarkerColor(colour);
  gEvidence->SetLineColor(colour);
  mEvidence->Add(gEvidence,"p");
  TGraph* gEvidencePE=new TGraph(nbins);
  graphName.clear();
  graphName << "g" << title << "_EvidencePerEvent";
  gEvidencePE->SetName (graphName.str().c_str());
  gEvidencePE->SetTitle(graphName.str().c_str());
  gEvidencePE->SetMarkerStyle(21);
  gEvidencePE->SetMarkerSize(0.5);
  gEvidencePE->SetMarkerColor(colour);
  gEvidencePE->SetLineColor(colour);
  mEvidencePerEvent->Add(gEvidencePE,"p");

  //cout << "filling data" << endl;


  // loop again over fitResults and extract all info simultaneously
  for(unsigned int i=0;i<nbins;++i){
    intree->GetEntry(i);
    // loop through waves
    it=wavesinthisfit.begin();
    while(it!=wavesinthisfit.end()){
      TMultiGraph* mg=mIntensities[*it];
      TGraphErrors* g=dynamic_cast<TGraphErrors*>(mg->GetListOfGraphs()->Last());
      g->SetPoint(i,
		  result->massBinCenter()*0.001,
		  result->intensity(it->c_str()));
      g->SetPointError(i,
		       binwidth*0.5,
		       result->intensityErr(it->c_str()));
      ++it;
    }
    gLikeli->SetPoint(i,
		      result->massBinCenter()*0.001,
		      result->logLikelihood());
    gLikeliPE->SetPoint(i,
			result->massBinCenter()*0.001,
			result->logLikelihood()/result->nmbEvents());
    gEvidence->SetPoint(i,
			result->massBinCenter()*0.001,
			result->evidence());
    gEvidencePE->SetPoint(i,
			result->massBinCenter()*0.001,
			result->evidence()/result->nmbEvents());
  }
  
  //cout << "writing meta" << endl;
  // write MetaInfo
  fitResultMetaInfo meta(filename,
			 title,
			 colour,treename,
			 branchname);
  meta.setLikelihoods(logli,logliperevt,evi,eviperevt);
  meta.setBinRange(mass_min-binwidth*0.5,mass_max+binwidth*0.5,nbins);
  mResultMetaInfo.push_back(meta);
  mMinEvidence=(mMinEvidence*ifit+evi)/(ifit+1);
    
  // cleanup
  infile->Close();
  cerr << endl;
  
  
}



bool 
pwaPlotter::registerWave(const std::string& wavename){
  pair<set<string>::iterator,bool> inserted=mWavenames.insert(wavename);
  if(inserted.second){ // we had a true insterion
    cerr << "New wave ("<<mWavenames.size()<<"): " << wavename << endl;
    // create intensity graph:
    mIntensities[wavename]=new TMultiGraph();
    mIntensities[wavename]->SetTitle(wavename.c_str());
    mIntensities[wavename]->SetName(wavename.c_str());
  }
  
  return inserted.second;
}


void 
pwaPlotter::produceDensityPlots(){
 map<string,TMultiGraph*>::iterator it=mIntensities.begin();
  while(it!=mIntensities.end()){
    // get ranges
    TList* graphs=it->second->GetListOfGraphs();
    unsigned int ng=graphs->GetEntries();
    double xmin=1E9;double ymin=1E9;
    double xmax=0; double ymax=-1E9;
    unsigned int nbins=0;
    for(unsigned int ig=0;ig<ng;++ig){
      TPwaFitGraphErrors* g=dynamic_cast<TPwaFitGraphErrors*>(graphs->At(ig));
      double xmin1,xmax1,ymin1,ymax1;
      g->ComputeRange(xmin1,ymin1,xmax1,ymax1);

      if(ymin > ymin1)ymin=ymin1;
      if(ymax < ymax1)ymax=ymax1;
      unsigned int ifit=g->fitindex;
      if(xmin > mResultMetaInfo[ifit].mMin)xmin=mResultMetaInfo[ifit].mMin;
      if(xmax < mResultMetaInfo[ifit].mMax)xmax=mResultMetaInfo[ifit].mMax;
      if(nbins< mResultMetaInfo[ifit].mNumBins)
	nbins=mResultMetaInfo[ifit].mNumBins;
    }
    double r=fabs(ymax-ymin)*0.1;
    // create 2D Histogram:
    string name="d";name.append(it->first);
    
    //cerr << ymin << " .. " << ymax << endl;
    TH2D* h=new TH2D(name.c_str(),name.c_str(),
		     nbins,xmin,xmax,
		     400,ymin-r,ymax+r);
    
    mIntensityDensityPlots[it->first]=h;
    
    // fill histo
    for(unsigned int ig=0;ig<ng;++ig){
      TPwaFitGraphErrors* g=dynamic_cast<TPwaFitGraphErrors*>(graphs->At(ig));
      unsigned int ifit=g->fitindex;
      double w=(mResultMetaInfo[ifit].mTotalPerEventEvidence)/(double)mResultMetaInfo[ifit].mNumBins;
      //double likeli=0;
      //cout << "weight: "<<TMath::Exp(w)<<endl;
      if(w==w)h=drawDensity(g,h,TMath::Exp(w));
    }
    ++it;
    // rescale each x-bin
    for(unsigned int ibin=1; ibin<=nbins; ++ibin){
      // get maximum bin in y for this x-bin
      double max=0;
      for(unsigned int iy=0;iy<400;++iy){
	unsigned int bin = h->GetBin(ibin,iy);
	double val=h->GetBinContent(bin);
	if(val>max)max=val;
      }
      if(max!=0 && max==max){
	for(unsigned int iy=0;iy<400;++iy){
	  unsigned int bin = h->GetBin(ibin,iy);
	  h->SetBinContent(bin,h->GetBinContent(bin)/max);
	}
      }
    }
    
    
  }// end loop over waves

}



void 
pwaPlotter::writeAll(std::string filename){
  TFile* outfile=TFile::Open(filename.c_str(),"RECREATE");
  if(outfile!=0 && !outfile->IsZombie()){
    writeAll(outfile);
    outfile->Close();
  }
  else{
    cerr << "Error opening file " << filename << endl;
  }
}

void 
pwaPlotter::writeAll(TFile* outfile){
   outfile->cd();
   // write evidence and loglikelihoods
   mLogLikelihood->Write();
   mLogLikelihoodPerEvent->Write();
   mEvidence->Write();
   mEvidencePerEvent->Write();
   writeAllIntensities(outfile);
}


void 
pwaPlotter::writeAllIntensities(std::string filename){
  TFile* outfile=TFile::Open(filename.c_str(),"RECREATE");
  if(outfile!=0 && !outfile->IsZombie()){
    writeAllIntensities(outfile);
    outfile->Close();
  }
  else{
    cerr << "Error opening file " << filename << endl;
  }
}



void 
pwaPlotter::writeAllIntensities(TFile* outfile){
  outfile->cd();
  map<string,TMultiGraph*>::iterator it=mIntensities.begin();
  while(it!=mIntensities.end()){
    it->second->Write();
    ++it;
  }
  map<string,TH2D*>::iterator itd=mIntensityDensityPlots.begin();
  while(itd!=mIntensityDensityPlots.end()){
    itd->second->Write();
    ++itd;
  } 
}






// void
// pwaPlotter::plotIntensity(const std::string& wavename, TTree* tr){
//   stringstream drawExpr;
//   drawExpr << branchName << ".intensity(\"" << waveName << "\"):"
// 	   << branchName << ".intensityErr(\"" << waveName << "\"):"
// 	   << branchName << ".massBinCenter() >> h" << waveName << "_" << i;
//   cout << "    running TTree::Draw() expression '" << drawExpr.str() << "' "
//        << "on tree '" << trees[i]->GetName() << "', '" << trees[i]->GetTitle() << "'" << endl;
  
//   cerr << "Drawing" << endl;
  
//   try{
//     trees[i]->Draw(drawExpr.str().c_str(), selectExpr.c_str(), "goff");
//   }
//   catch(std::exception&){
//     cerr << "Cought Exception" << endl;
//     continue;
//   }
  
  
  
//   // extract data from TTree::Draw() result and build graph
//   const int nmbBins = trees[i]->GetSelectedRows();
//   vector<double> x(nmbBins), xErr(nmbBins);
//   vector<double> y(nmbBins), yErr(nmbBins);
//   for (int j = 0; j < nmbBins; ++j) {
//     x   [j] = trees[i]->GetV3()[j] * 0.001;  // convert mass to GeV
//     xErr[j] = 0;
//     y   [j] = trees[i]->GetV1()[j] * normalization;  // scale intensities
//     yErr[j] = trees[i]->GetV2()[j] * normalization;  // scale intensity errors
//   }
  
  
  
//   TGraphErrors* g = new TGraphErrors(nmbBins,
// 				     &(*(x.begin())),      // mass
// 				     &(*(y.begin())),      // intensity
// 				     &(*(xErr.begin())),   // mass error
// 				     &(*(yErr.begin())));  // intensity error
//   {
//     stringstream graphName;
//     graphName << ((graphTitle == "") ? waveName : graphTitle) << "_" << i;
//     g->SetName (graphName.str().c_str());
//     g->SetTitle(graphName.str().c_str());
//   }
//   g->SetMarkerStyle(21);
//   g->SetMarkerSize(0.5);
//   if (graphColors) {
//     g->SetMarkerColor(graphColors[i]);
//     g->SetLineColor  (graphColors[i]);
//   }
//   graph->Add(g);
  
//   // compute maximum for y-axis
//   for (int j = 0; j < nmbBins; ++j)
//     if (maxYVal < (y[j] + yErr[j]))
//       maxYVal = y[j] + yErr[j];
//   const double yMean = g->GetMean(2);
//   if (maxYMean < yMean)
//     maxYMean = yMean;
// }