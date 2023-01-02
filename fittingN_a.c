//This is a ROOT program
//written by Dat Tran, 11/2022
//performs linear fitting of DP's observables'data (log-log scale) and outputs the slopes which correspond to critical exponents.
//Input: Text file with measurements of an observable
//Output: Plot on log-log scale with linear fitting and fit parameters displayed

#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

void fittingN_a(){
   const int N = 50;
   double x[N] = {0.};
   double xe[N] = {0.};
   double y[N] = {0.};
   double ye[N] = {0.};
   double time[N-2] = {0.};
   double N_a[N-2] = {0.};
   double err[N-2] = {0.};
   char filename[50] = "numberOfActiveSites-L50-p0.64.txt";

   double value = 0.;
   double t = 0.;
   double error = 0.;
   FILE *inFile1 = fopen(filename, "r");

   for (int i = 0; i < N; i++) {
        fscanf(inFile1, "%lf %lf %lf", &t, &value, &error);
        x[i] = log(t);
        y[i] = log(value);
	ye[i] = log(error);
   }

   for (int i = 0; i < N - 2; i++) {
      time[i] = x[i+2];
      N_a[i] = y[i+2];
      err[i] = ye[i+2];
   }
   auto gr1 = new TGraphErrors(N-2,time,N_a,xe,err);
   gr1->SetName("gr1");
   gr1->SetLineColor(kBlue);

   auto c0 = new TCanvas("c1","average N_{a} at critical percolation probability over 10,000 runs",200,10,700,500);
   gr1->SetTitle("average N_{a} at critical percolation probability over 10,000 runs");
   gr1->GetXaxis()->SetTitle("time - t [mcs]");
   gr1->GetYaxis()->SetTitle("Number of Active Sites - <N_{a}>");
   gr1->Draw("AP");
   c0->Update();
}

