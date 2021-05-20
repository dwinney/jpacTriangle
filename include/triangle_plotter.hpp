// A custom plotter object
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _TRIPLOT_
#define _TRIPLOT_

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include "TLine.h"
#include "constants.hpp"

namespace jpacTriangle
{
    class triangle_plotter : public jpacGraph1D
    {
        public:
        triangle_plotter()
        : jpacGraph1D()
        {
            line_width = 7;
            SetLegendOffset(0.3, 0.1);
            jpacStyle->SetTitleOffset(1.15, "x");
        };

        double th1, th2, th3;
        void inline SetVerticals(double t, double p, double s)
        {
            th1 = t;
            th2 = p;
            th3 = s;
        };

        inline void AddEntry(std::vector<double> x, std::vector<std::complex<double>> * fx)
        {
            std::vector<double> real0 = vec_real(fx[0]);
            std::vector<double> imag0 = vec_imag(fx[0]);

            jpacGraph1D::AddEntry(x, real0, "Real");
            jpacGraph1D::AddEntry(x, imag0, "Imaginary");

            std::vector<double> real1 = vec_real(fx[1]);
            std::vector<double> imag1 = vec_imag(fx[1]);
            gR = new TGraph(x.size(), &(x[0]), &(real1[0]));
            gI = new TGraph(x.size(), &(x[0]), &(imag1[0]));
        };

        private:
        TGraph* gR = NULL;
        TGraph* gI = NULL;

        inline void AddExtra()
        {
            DrawVerticals();
            DrawCompare();
            canvas->Modified();
            canvas->Draw();
        };

        inline void DrawCompare()
        {
            gR->SetLineWidth(line_width);
            gR->SetLineColor(1);
            gR->SetLineStyle(4);
            gR->Draw("same");

            gI->SetLineWidth(line_width);
            gI->SetLineColor(1);
            gI->SetLineStyle(4);
            gI->Draw("same");
        };

        inline void DrawVerticals()
        {
            double yMax, yMin;
            if (yCustom == true)
            {
                yMax = yhigh;
                yMin = ylow;
            }
            else
            {
                yMax = yAxis->GetXmax();
                yMin = yAxis->GetXmin();
            }

            TLine* V1 = new TLine(th1, yMax, th1, yMin);
            V1->SetLineStyle(3);
            V1->Draw("same");

            TLine* V2 = new TLine(th2, yMin, th2, yMax);
            V2->SetLineStyle(3);
            V2->Draw("Same");

            TLine* V3 = new TLine(th3, yMin, th3, yMax);
            V3->SetLineStyle(3);
            V3->Draw("Same");
        };
    };
};

#endif