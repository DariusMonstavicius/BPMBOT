#ifndef _PLOT_H_
#define _PLOT_H_

#include <vector>

void PlotData(float *x, float *y, int length, std::string title, int figure = 1, int usex = 0);
void PlotData(std::vector<int> y, std::string title, int figure = 1);
void PlotEverything(float *tempdatabuffer, int tempdatabufferlen, std::vector<int> neighbors, std::vector<int> maxindices, std::vector<int> distances, std::vector<int> cleandistances, float bpm);

#endif
