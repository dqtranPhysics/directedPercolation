//This is a ROOT program
//written by Dat Tran, 11/2022
//performs linear fitting of DP's observables'data (log-log scale) and outputs the slopes which correspond to critical exponents.
//Input: Text file with measurements of an observable
//Output: Plot on log-log scale with linear fitting and fit parameters displayed

#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

void fitting(){
   const int N = 50;
   double x[N] = {0.};
   double y[N] = {0.};
   double time[N-2] = {0.};
   double N_a[N-2] = {0.};
   char filename[50] = "numberOfActiveSites-L50-p0.64.txt";

   double value = 0.;
   double t = 0.;
   FILE *inFile1 = fopen(filename, "r");

   for (int i = 0; i < N; i++) {
        fscanf(inFile1, "%lf %lf", &t, &value);
        x[i] = log2(t);
        y[i] = log2(value);
	cout << x[i] << " " << y[i] << endl;
   }

   for (int i = 0; i < N - 2; i++) {
      time[i] = x[i+2];
      N_a[i] = y[i+2];
   }
   auto gr1 = new TGraph(N-2,time,N_a);
   gr1->SetName("gr1");
   gr1->SetLineColor(kBlue);

   auto c0 = new TCanvas("c1","number of active sites at critical percolation probability",200,10,700,500);
   gr1->GetXaxis()->SetTitle("time - t [mcs]");
   gr1->GetYaxis()->SetTitle("Number of Active Sites - <N_{a}>");
   gr1->Draw("AP");
   c0->Update();
}

