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
//-----------------------------------------------------------
//
// Description:
//      Draws intensity graph for single wave from tree.
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <iostream>
#include <sstream>

#include <TAxis.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TTree.h>

#include "fitResult.h"
#include "reportingUtils.hpp"


using namespace std;
using namespace rpwa;


// signature with wave name
TMultiGraph*
plotIntensity(const unsigned int nmbTrees,                        // number of fitResult trees
              TTree**            trees,                           // array of fitResult trees
              const string&      waveName,                        // wave name
              const bool         saveEps       = false,           // if set, EPS file with name wave ID is created
              const int*         graphColors   = NULL,            // array of colors for graph line and marker
              const double*      graphScales   = NULL,            // array of scales for graphgroups
              const bool         drawLegend    = true,            // if set legend is drawn
              const string&      graphTitle    = "",              // name and title of graph (default is wave ID)
              const char*        drawOption    = "AP",            // draw option for graph
              const double       normalization = 1,               // scale factor for intensities
              const double       yAxisRangeMax = 0,               // if != 0; range of y-axis is limited to this value
              const string&      selectExpr    = "",              // TTree::Draw() selection expression
              const string&      branchName    = "fitResult_v2")  // fitResult branch name
{

	for (unsigned int i = 0; i < nmbTrees; ++i)
		if (!trees[i]) {
			printErr << "null pointer to tree[" << i << "]. aborting." << endl;
			return 0;
		}
	printInfo << "plotting wave intensity for wave '" << waveName;
	if (selectExpr != "")
		cout << " using selection criterion '" << selectExpr << "'";
	cout << endl;

	// create multiGraph
	TMultiGraph* graph = new TMultiGraph();
	{
		if (graphTitle != "") {
			graph->SetName (graphTitle.c_str());
			graph->SetTitle(graphTitle.c_str());
		} else {
			stringstream suffix;
			//suffix << " [" << "no waveIndex" << "]";
			graph->SetName (waveName.c_str());
			graph->SetTitle((waveName + suffix.str()).c_str());
		}
	}

	// fill multiGraph
	double maxY = 0;
	for (unsigned int i = 0; i < nmbTrees; ++i) {

		// get wave index for this tree
		fitResult* massBin = new fitResult();
		trees[i]->SetBranchAddress(branchName.c_str(), &massBin);
		vector<int> waveIndexThisTree; // vector of wave indices per bin
		bool found(false);
		for (int imassBin = 0; imassBin < trees[i]->GetEntries(); imassBin++){
			trees[i]->GetEntry(imassBin);
			waveIndexThisTree.push_back(massBin->waveIndex(waveName));
			if (waveIndexThisTree[imassBin] >= 0){
				found = true;
			}
		}
		if (!found) {
			printInfo << "cannot find wave '" << waveName << "' in tree '" << trees[i]->GetTitle() << "'. "
			          << "skipping." << endl;
			continue;
		}

		// running bin by bin since wave name does necessarily correspond to a fixed wave number

		// build and run TTree::Draw() expression
		stringstream drawExpr;
		drawExpr << branchName << ".intensity(\"" << waveName << "\"):"
		         << branchName << ".intensityErr(\"" << waveName << "\"):"
		         << branchName << ".massBinCenter() >> h" << waveName << "_" << i;
		cout << "    running TTree::Draw() expression '" << drawExpr.str() << "' "
		     << "on tree '" << trees[i]->GetName() << "', '" << trees[i]->GetTitle() << "'" << endl;
		trees[i]->Draw(drawExpr.str().c_str(), selectExpr.c_str(), "goff");

		// extract data from TTree::Draw() result and build graph
		const int nmbBins = trees[i]->GetSelectedRows();
		vector<double> x, xErr;
		vector<double> y, yErr;

		double scale = normalization;
		if(graphScales) {
		  scale = scale * graphScales[i];
		}

		for (int j = 0; j < nmbBins; ++j) {
//			if(trees[i]->GetV3()[j] * 0.001 <= 1.35) continue;
//			if(trees[i]->GetV3()[j] * 0.001 >= 4.05) continue;

//			if(trees[i]->GetV3()[j] * 0.001 <= 1.25) continue;
//			if(trees[i]->GetV3()[j] * 0.001 >= 3.15) continue;
			y	.push_back(trees[i]->GetV1()[j] * scale);  // scale intensities
			// if (y[y.size()-1] == 0){ // remove 0 entries
			// 	y.pop_back();
			// 	continue;
			// }
			x	.push_back(trees[i]->GetV3()[j] * 0.001);  // convert mass to GeV
//			xErr.push_back(0.015);
			xErr.push_back(0);
			yErr.push_back(trees[i]->GetV2()[j] * scale);  // scale intensity errors
		}
		TGraphErrors* g = new TGraphErrors(x.size(),
		                                   &(*(x.begin())),      // mass
		                                   &(*(y.begin())),      // intensity
		                                   &(*(xErr.begin())),   // mass error
		                                   &(*(yErr.begin())));  // intensity error

		// beautify graph
		stringstream graphName;
		graphName << graph->GetName() << "_" << i;
		g->SetName (graphName.str().c_str());
		g->SetTitle(graphName.str().c_str());
		g->SetMarkerStyle(21);
		g->SetMarkerSize(0.5);
		if (graphColors) {
			g->SetMarkerColor(graphColors[i]);
			g->SetLineColor  (graphColors[i]);
		}
		graph->Add(g);

		// compute maximum for y-axis
		for (unsigned int j = 0; j < x.size(); ++j) {
			const double val = y[j] + yErr[j];
			if (maxY < val)
				maxY = val;
		}
	}
	cout << "    maximum intensity for graph " << graph->GetName() << " is " << maxY << endl;

	if ((yAxisRangeMax > 0) && (maxY > yAxisRangeMax))
		maxY = yAxisRangeMax;
	graph->SetMinimum(-maxY * 0.1);
	graph->SetMaximum( maxY * 1.1);
	// draw graph
	graph->Draw(drawOption);
	graph->GetXaxis()->SetTitle("Mass [GeV]");
//	graph->GetXaxis()->SetRangeUser(1.3, 3.1);
	graph->GetYaxis()->SetTitle("Intensity");
	TLine line;
	line.SetLineStyle(3);
	line.DrawLine(graph->GetXaxis()->GetXmin(), 0, graph->GetXaxis()->GetXmax(), 0);
	gPad->Update();

	// add legend
	if (drawLegend && (nmbTrees > 1)) {
		TLegend* legend = new TLegend(0.65,0.80,0.99,0.99);
		legend->SetFillColor(10);
		legend->SetBorderSize(1);
		legend->SetMargin(0.2);
		for (unsigned int i = 0; i < nmbTrees; ++i) {
			TGraph* g = static_cast<TGraph*>(graph->GetListOfGraphs()->At(i));
			legend->AddEntry(g, trees[i]->GetTitle(), "LPE");
		}
		legend->Draw();
	}

	// create EPS file
	if (saveEps)
		gPad->SaveAs(((string)graph->GetName() + ".eps").c_str());

	return graph;
}

// signature with wave index

TMultiGraph*
plotIntensity(const unsigned int nmbTrees,                       // number of fitResult trees
              TTree**            trees,                          // array of fitResult trees
              const int          waveIndex,                      // wave index
              const bool         saveEps       = false,          // if set, EPS file with name wave ID is created
              const int*         graphColors   = NULL,           // array of colors for graph line and marker
              const double*      graphScales   = NULL,           // array of scales for graphgroups
              const bool         drawLegend    = true,           // if set legend is drawn
              const string&      graphTitle    = "",             // name and title of graph (default is wave ID)
              const char*        drawOption    = "AP",           // draw option for graph
              const double       normalization = 1,              // scale factor for intensities
              const double       yAxisRangeMax = 0,              // if != 0; range of y-axis is limited to this value
              const string&      selectExpr    = "",             // TTree::Draw() selection expression
              const string&      branchName    = "fitResult_v2") // fitResult branch name
{
	for (unsigned int i = 0; i < nmbTrees; ++i)
		if (!trees[i]) {
			printErr << "null pointer to tree[" << i << "]. aborting." << endl;
			return 0;
		}
	string waveName;
	{
		// get wave name (assume same wave set in all bins)
		fitResult* massBin = new fitResult();
		trees[0]->SetBranchAddress(branchName.c_str(), &massBin);
		trees[0]->GetEntry(0);
		waveName = massBin->waveName(waveIndex).Data();
	}
	printInfo << "plotting wave intensity for wave '" << waveName << "' [" << waveIndex << "]";
	if (selectExpr != "")
		cout << " using selection criterion '" << selectExpr << "'";
	cout << endl;

	// create multiGraph
	TMultiGraph* graph = new TMultiGraph();
	{
		if (graphTitle != "") {
			graph->SetName (graphTitle.c_str());
			graph->SetTitle(graphTitle.c_str());
		} else {
			stringstream suffix;
			suffix << " [" << waveIndex << "]";
			graph->SetName (waveName.c_str());
			graph->SetTitle((waveName + suffix.str()).c_str());
		}
	}

	// fill multiGraph
	double maxY = 0;
	for (unsigned int i = 0; i < nmbTrees; ++i) {

		// get wave index for this tree (assume same wave set in all bins)
		fitResult* massBin = new fitResult();
		trees[i]->SetBranchAddress(branchName.c_str(), &massBin);
		trees[i]->GetEntry(0);
		const int waveIndexThisTree = massBin->waveIndex(waveName);
		if (waveIndexThisTree < 0) {
			printInfo << "cannot find wave '" << waveName << "' in tree '" << trees[i]->GetTitle() << "'. "
			          << "skipping." << endl;
			continue;
		}

		// build and run TTree::Draw() expression
		stringstream drawExpr;
		drawExpr << branchName << ".intensity(" << waveIndexThisTree << "):"
		         << branchName << ".intensityErr(" << waveIndexThisTree << "):"
		         << branchName << ".massBinCenter() >> h" << waveName << "_" << i;
		cout << "    running TTree::Draw() expression '" << drawExpr.str() << "' "
		     << "on tree '" << trees[i]->GetName() << "', '" << trees[i]->GetTitle() << "'" << endl;
		trees[i]->Draw(drawExpr.str().c_str(), selectExpr.c_str(), "goff");

		// extract data from TTree::Draw() result and build graph
		const int nmbBins = trees[i]->GetSelectedRows();
		vector<double> x(nmbBins), xErr(nmbBins);
		vector<double> y(nmbBins), yErr(nmbBins);

		double scale = normalization;
    if (graphScales) {
      scale = scale * graphScales[i];
    }

		for (int j = 0; j < nmbBins; ++j) {
			x   [j] = trees[i]->GetV3()[j] * 0.001;  // convert mass to GeV
			xErr[j] = 0;
			y   [j] = trees[i]->GetV1()[j] * scale;  // scale intensities
			yErr[j] = trees[i]->GetV2()[j] * scale;  // scale intensity errors
		}
		TGraphErrors* g = new TGraphErrors(nmbBins,
		                                   &(*(x.begin())),      // mass
		                                   &(*(y.begin())),      // intensity
		                                   &(*(xErr.begin())),   // mass error
		                                   &(*(yErr.begin())));  // intensity error

		// beautify graph
		stringstream graphName;
		graphName << graph->GetName() << "_" << i;
		g->SetName (graphName.str().c_str());
		g->SetTitle(graphName.str().c_str());
		g->SetMarkerStyle(21);
		g->SetMarkerSize(0.5);
		if (graphColors) {
			g->SetMarkerColor(graphColors[i]);
			g->SetLineColor  (graphColors[i]);
		}
		graph->Add(g);

		// compute maximum for y-axis
		for (int j = 0; j < nmbBins; ++j) {
			const double val = y[j] + yErr[j];
			if (maxY < val)
				maxY = val;
		}
	}
	cout << "    maximum intensity for graph " << graph->GetName() << " is " << maxY << endl;

	if ((yAxisRangeMax > 0) && (maxY > yAxisRangeMax))
		maxY = yAxisRangeMax;
	graph->SetMinimum(-maxY * 0.1);
	graph->SetMaximum( maxY * 1.1);
	// draw graph
	graph->Draw(drawOption);
	graph->GetXaxis()->SetTitle("Mass [GeV]");
	graph->GetYaxis()->SetTitle("Intensity");
	TLine line;
	line.SetLineStyle(3);
	line.DrawLine(graph->GetXaxis()->GetXmin(), 0, graph->GetXaxis()->GetXmax(), 0);
	gPad->Update();

	// add legend
	if (drawLegend && (nmbTrees > 1)) {
		TLegend* legend = new TLegend(0.65,0.80,0.99,0.99);
		legend->SetFillColor(10);
		legend->SetBorderSize(1);
		legend->SetMargin(0.2);
		for (unsigned int i = 0; i < nmbTrees; ++i) {
			TGraph* g = static_cast<TGraph*>(graph->GetListOfGraphs()->At(i));
			legend->AddEntry(g, trees[i]->GetTitle(), "LPE");
		}
		legend->Draw();
	}

	// create EPS file
	if (saveEps)
		gPad->SaveAs(((string)graph->GetName() + ".eps").c_str());

	return graph;
}


TMultiGraph*
plotIntensity(vector<TTree*>&             trees,                             // array of fitResult trees
              const int                   waveIndex,                         // wave index
              const bool                  saveEps       = false,             // if set, EPS file with name wave ID is created
              const vector<int>&          graphColors   = vector<int>(),     // array of colors for graph line and marker
              const vector<double>&       graphScales   = vector<double>(),  // array of scales for graphgroups
              const bool                  drawLegend    = true,              // if set legend is drawn
              const string&               graphTitle    = "",                // name and title of graph (default is wave ID)
              const char*                 drawOption    = "AP",              // draw option for graph
              const double                normalization = 1,                 // scale factor for intensities
              const double                yAxisRangeMax = 0,                 // if != 0; range of y-axis is limited to this value
              const string&               selectExpr    = "",                // TTree::Draw() selection expression
              const string&               branchName    = "fitResult_v2")    // fitResult branch name
{
	return plotIntensity(trees.size(), &(*(trees.begin())), waveIndex, saveEps,
	                     &(*(graphColors.begin())), &(*(graphScales.begin())), drawLegend, graphTitle, drawOption,
	                     normalization, yAxisRangeMax, selectExpr, branchName);
}


TMultiGraph*
plotIntensity(TTree*             tree,                            // fitResult tree
              const int          waveIndex,                       // wave index
              const bool         saveEps       = false,           // if set, EPS file with name wave ID is created
              const int          graphColor    = kBlack,          // color of line and marker
              const bool         drawLegend    = false,           // if set legend is drawn
              const string&      graphTitle    = "",              // name and title of graph (default is wave ID)
              const char*        drawOption    = "AP",            // draw option for graph
              const double       normalization = 1,               // scale factor for intensities
              const double       yAxisRangeMax = 0,               // if != 0; range of y-axis is limited to this value
              const string&      selectExpr    = "",              // TTree::Draw() selection expression
              const string&      branchName    = "fitResult_v2")  // fitResult branch name
{
	return plotIntensity(1, &tree, waveIndex, saveEps, &graphColor, NULL, drawLegend, graphTitle,
	                     drawOption, normalization,  yAxisRangeMax, selectExpr, branchName);
}


TMultiGraph*
plotIntensity(vector<TTree*>&            trees,                             // array of fitResult trees
              const string&              waveName,                          // wave name
              const bool                 saveEps       = false,             // if set, EPS file with name wave ID is created
              const vector<int>&         graphColors   = vector<int>(),     // array of colors for graph line and marker
              const vector<double>&      graphScales   = vector<double>(),  // array of scales for graphgroups
              const bool                 drawLegend    = true,              // if set legend is drawn
              const string&              graphTitle    = "",                // name and title of graph (default is wave ID)
              const char*                drawOption    = "AP",              // draw option for graph
              const double               normalization = 1,                 // scale factor for intensities
              const double               yAxisRangeMax = 0,                 // if != 0; range of y-axis is limited to this value
              const string&              selectExpr    = "",                // TTree::Draw() selection expression
              const string&              branchName    = "fitResult_v2")    // fitResult branch name
{
	return plotIntensity(trees.size(), &(*(trees.begin())), waveName, saveEps,
	                     &(*(graphColors.begin())), &(*(graphScales.begin())), drawLegend, graphTitle, drawOption,
	                     normalization, yAxisRangeMax, selectExpr, branchName);
}


TMultiGraph*
plotIntensity(TTree*             tree,                            // fitResult tree
              const string&      waveName,                        // wave name
              const bool         saveEps       = false,           // if set, EPS file with name wave ID is created
              const int          graphColor    = kBlack,          // color of line and marker
              const bool         drawLegend    = true,            // if set legend is drawn
              const string&      graphTitle    = "",              // name and title of graph (default is wave ID)
              const char*        drawOption    = "AP",            // draw option for graph
              const double       normalization = 1,               // scale factor for intensities
              const double       yAxisRangeMax = 0,               // if != 0; range of y-axis is limited to this value
              const string&      selectExpr    = "",              // TTree::Draw() selection expression
              const string&      branchName    = "fitResult_v2")  // fitResult branch name
{
	return plotIntensity(1, &tree, waveName, saveEps, &graphColor, NULL, drawLegend, graphTitle,
	                     drawOption, normalization, yAxisRangeMax, selectExpr, branchName);
}
