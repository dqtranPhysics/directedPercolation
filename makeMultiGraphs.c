#include <iostream>
#include <fstream>
using namespace std;

void makeMultiGraphs(){
   auto c45 = new TCanvas("c45","c45",200,10,600,400);
   const int N = 300;
   double x[N] = {0.};
   double y1[N], y2[N], y3[N], y4[N] = {0.};
   char filename1[50] = "numberOfActiveSites-L5-p0.600.txt";
   char filename2[50] = "numberOfActiveSites-L5-p0.640.txt";
   char filename3[50] = "numberOfActiveSites-L5-p0.650.txt";
   char filename4[50] = "numberOfActiveSites-L5-p0.700.txt";

   double value1,value2,value3, value4 = 0.;
   double error1,error2,error3,error4 = 0.;
   double t = 0.;
   FILE *inFile1 = fopen(filename1, "r");
   FILE *inFile2 = fopen(filename2, "r");
   FILE *inFile3 = fopen(filename3, "r");
   FILE *inFile4 = fopen(filename4, "r");

   for (int i = 0; i < N; i++) {
      fscanf(inFile1, "%lf %lf %lf", &t, &value1, &error1);
	   fscanf(inFile2, "%lf %lf %lf", &t, &value2, &error2);
	   fscanf(inFile3, "%lf %lf %lf", &t, &value3, &error3);
	   fscanf(inFile4, "%lf %lf %lf", &t,&value4, &error4);
      x[i] = t;
      y1[i] = value1;
      y2[i] = value2;
      y3[i] = value3;
      y4[i] = value4;
   }
	
	
 
   auto mg = new TMultiGraph();
 
   auto gr1 = new TGraph(N,x,y1);
   gr1->SetName("gr1");
   gr1->SetTitle("p = 0.6");
   gr1->SetLineColor(kBlue);

   auto gr2 = new TGraph(N,x,y2);
   gr2->SetName("gr2");
   gr2->SetTitle("p = 0.64");
   gr2->SetLineColor(kRed);

   auto gr3 = new TGraph(N,x,y3);
   gr3->SetName("gr3");
   gr3->SetTitle("p = 0.65");
   gr3->SetLineColor(kGreen);

   auto gr4 = new TGraph(N,x,y4);
   gr4->SetName("gr4");
   gr4->SetTitle("p = 0.7");
   gr4->SetLineColor(kBlack);

   //mg -> Add(gr1);
   mg -> Add(gr2);
   mg -> Add(gr3);
   mg -> Add(gr4);
   
   auto c0 = new TCanvas("c1","number of active sites for different percolation probabilities",200,10,700,500);
   gPad->SetLogx();
   gPad->SetLogy();
   mg->GetXaxis()->SetTitle("time - t [mcs]");
   mg->GetYaxis()->SetTitle("Number of Active Sites - <N_{a}>");
   mg->Draw("ALP");
   c0->BuildLegend();
   c0->Update();

}

