#include <iostream>
#include <fstream>
using namespace std;

void dataCollapse(){
   auto c45 = new TCanvas("c45","c45",200,10,600,400);
   const int N = 50;
   double theta = -0.314;
   double nu = 3.73;
   double x1[N], x2[N], x3[N], x4[N] = {0.};
   double y1[N], y2[N], y3[N], y4[N] = {0.};
   char filename1[50] = "numberOfActiveSites-L50-p0.63.txt";
   char filename2[50] = "numberOfActiveSites-L50-p0.635.txt";
   char filename3[50] = "numberOfActiveSites-L50-p0.645.txt";
   char filename4[50] = "numberOfActiveSites-L50-p0.65.txt";

   double value1,value2,value3, value4 = 0.;
   double error1, error2, error3, error4 = 0.;
   double t = 0.;
   FILE *inFile1 = fopen(filename1, "r");
   FILE *inFile2 = fopen(filename2, "r");
   FILE *inFile3 = fopen(filename3, "r");
   FILE *inFile4 = fopen(filename4, "r");

   for (int i = 0; i < N; i++) {
        fscanf(inFile1, "%lf %lf %lf", &t, &value1, &error1);
	fscanf(inFile2, "%lf %lf %lf", &t, &value2, &error2);
	fscanf(inFile3, "%lf %lf %lf", &t, &value3, &error3);
	fscanf(inFile4, "%lf %lf %lf", &t,&value4, &error3);
	x1[i] = pow(t*(0.64 - 0.63),nu);
	x2[i] = pow(t*(0.64 - 0.635),nu);
	x3[i] = pow(t*(0.645 - 0.64),nu);
	x4[i] = pow(t*(0.65 - 0.64),nu);
	y1[i] = pow(value1*t,theta);
	y2[i] = pow(value2*t,theta);
	y3[i] = pow(value3*t,theta);
	y4[i] = pow(value4*t,theta);
   }
	
	
 
   auto mg = new TMultiGraph();
 
   auto gr1 = new TGraph(N,x1,y1);
   gr1->SetName("gr1");
   gr1->SetTitle("p = 0.63");
   gr1->SetLineColor(kBlue);

   auto gr2 = new TGraph(N,x2,y2);
   gr2->SetName("gr1");
   gr2->SetTitle("p = 0.635");
   gr2->SetLineColor(kRed);

   auto gr3 = new TGraph(N,x3,y3);
   gr3->SetName("gr3");
   gr3->SetTitle("p = 0.645");
   gr3->SetLineColor(kGreen);

   auto gr4 = new TGraph(N,x4,y4);
   gr4->SetName("gr4");
   gr4->SetTitle("p = 0.65");
   gr4->SetLineColor(kBlack);

   mg -> Add(gr1);
   mg -> Add(gr2);
   mg -> Add(gr3);
   mg -> Add(gr4);
   
   auto c0 = new TCanvas("c1","N_{a} curves collapse for #nu_{||} = 3.73",200,10,700,500);
   gPad->SetLogx();
   gPad->SetLogy();
   mg -> SetTitle("N_{a} curves collapse for #nu_{||} = 3.73");
   mg->GetXaxis()->SetTitle("time - t(p - p_{c})^{#nu_{||}} [mcs]");
   mg->GetYaxis()->SetTitle("Number of Active Sites - N_{a}(t)*t^{-#theta}");
   mg->Draw("ALP");
   c0->BuildLegend();
   c0->Update();

}

