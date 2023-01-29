#define _ITERATOR_DEBUG_LEVEL 0
#include <string>

#include <matplot/matplot.h>
using namespace matplot;


void PlotData(float *x, float *y, int length, std::string title, int figurenum, int usex) {
    static std::map<int, figure_handle> figures;

    std::vector<double> xd;
    std::vector<double> yd;

    xd.resize(length);
    yd.resize(length);

    for (int i = 0; i < length; i++) {
        if (usex) {
            xd[i] = (double)(x[i]);
        }
        else {
            xd[i] = (double)i;
        }
        yd[i] = (double) (y[i]);
    }
    if (figures[figurenum] == 0) {
        figures[figurenum] = figure(true);
    }

    figures[figurenum]->position(figurenum * 600, 0);
    figures[figurenum]->title(title);
 
    figure(figures[figurenum]);
    plot(xd, yd, "-o");
    figures[figurenum]->draw();
    
}

void PlotData(std::vector<int> y, std::string title, int figure = 1) {
    float x;
    float* yf = (float *) malloc(y.size() * sizeof(float));
    for (int i = 0; i < y.size(); i++) {
        yf[i] = (float)y[i];
    }
    PlotData(&x, yf, y.size(), title, figure, 0);
    free(yf);
}

void PlotEverything(float * tempdatabuffer, int tempdatabufferlen, std::vector<int> neighbors, std::vector<int> maxindices, std::vector<int> distances, std::vector<int> cleandistances, float bpm) {


    static int figurenum = 0;
    int i;

    if (figurenum > 5) {
        figurenum = 0;
    }


    static std::map<int, figure_handle> figures;
    //allocate new figure if necessary
    if (figures[figurenum] == 0) {
        figures[figurenum] = figure(true);
    }

    figure_handle f = figures[figurenum];
    //plot sound data
    std::vector<double> xd;
    std::vector<double> yd;

    xd.resize(tempdatabufferlen);
    yd.resize(tempdatabufferlen);

    for ( i = 0; i < tempdatabufferlen; i++) {
        xd[i] = (double)i;
        yd[i] = (double)(tempdatabuffer[i]);
    }

    f->position(figurenum * 600, 600);
    std::string title1 = "BPM ";
    std::string title2;
    title2 = std::to_string(bpm);

    std::string title = title1 + title2;
    f->title(title);

    figure(f);
    hold(false);

    plot(xd, yd);
    hold(true);

    //draw neighbors
    for (i = 0; i < neighbors.size(); i++) {
        xd[i] = (double)neighbors[i];
        yd[i] = (double)(tempdatabuffer[i]);
    }

    //plot(xd, yd);

    //draw maxindices
    for (i = 0; i < maxindices.size(); i++) {
        xd[i] = (double)maxindices[i];
        yd[i] = (double)(tempdatabuffer[i]);
        
        std::string legendstr;
        if (distances.size() > i) {
            std::string distancestr = std::to_string(distances[i]/1000);
            distancestr = distancestr + "k ";
            legendstr = legendstr + distancestr;
        }
        auto [t, a] = textarrow(xd[i], -0.05, xd[i], yd[i], legendstr);
    }
    ylim({ -0.1, inf });





    figures[figurenum]->draw();


    figurenum++;
}
