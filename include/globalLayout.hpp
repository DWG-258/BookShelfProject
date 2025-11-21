#include "placedata.hpp"
#include "objects.hpp"
#include <algorithm>





void computeBinDensity(
    std::vector<std::vector<Bin>>& bins,
    const vector<shared_ptr<Module>> &modules,
    int m, double binWidth, double binHeight,
    double targetDensity)
{
    // 1. 遍历单元，累加网格贡献
    for (auto& node : modules) {    
        int x_start = std::max(0, int(node->center.x/ binWidth));
        int x_end   = std::min(m-1, int((node->center.x + node->width) / binWidth));
        int y_start = std::max(0, int(node->center.y / binHeight));
        int y_end   = std::min(m-1, int((node->center.y + node->height) / binHeight));

        for (int i = x_start; i <= x_end; ++i) {
            for (int j = y_start; j <= y_end; ++j) {
                double overlapArea = binWidth * binHeight; // 可计算真实重叠
                double scaleX = node->width/ binWidth;
                double scaleY = node->height / binHeight;
                double contribution = overlapArea * scaleX * scaleY;

                if (!node->isFixed&&!node->isFiller&&!node->isMacro)
                    bins[i][j].nodeDensity += contribution;
                else if (node->isMacro)
                    bins[i][j].nodeDensity += contribution * targetDensity;
                else if (node->isFiller)
                    bins[i][j].fillerDensity += contribution * targetDensity;
            }
        }
    }

    // 2. 计算总密度
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            double totalArea = bins[i][j].nodeDensity
                             + bins[i][j].fillerDensity
                             + bins[i][j].darkDensity
                             + bins[i][j].terminalDensity;
            bins[i][j].eDensity = totalArea / (binWidth * binHeight);
        }
    }
}

std::vector<std::vector<double>> computeDCT(
    const std::vector<std::vector<Bin>>& bins, int m)
{
    std::vector<std::vector<double>> a(m, std::vector<double>(m, 0.0));
    for (int j = 0; j < m; ++j) {
        for (int k = 0; k < m; ++k) {
            double sum = 0.0;
            for (int x = 0; x < m; ++x) {
                for (int y = 0; y < m; ++y) {
                    sum += bins[x][y].eDensity
                         * cos(M_PI*j*(x+0.5)/m)
                         * cos(M_PI*k*(y+0.5)/m);
                }
            }
            a[j][k] = 2.0 * sum / m; // DCT-II
        }
    }
    return a;
}