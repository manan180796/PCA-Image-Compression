#include <omp.h>
#include <chrono>
#include <iostream>
#include "Image.hpp"
#include "eigen.hpp"
#include "pca.hpp"
using namespace std;

int main() {
    // cv::Mat image = readImage("lights.jpg");
    // cv::Mat image = readImage("office-365-logo-small-1.png");
    // DenseMatrix<double> data = Mat2DenseMatrix(image);
    // saveRawImage(data, "lights.jpg.raw");
    // SaveCompressedImage(data, "lights.jpg.compressed.400", 400);
    // auto start = std::chrono::high_resolution_clock::now();
    // DenseMatrix<double> V = PCA(data);
    // int k = 50;
    // DenseMatrix<double> V_ = V.GetSubMatrix(0, V.Columns() - k, V.Rows(), k);
    DenseMatrix<double> newData =
        ReadCompressedImage("lights.jpg.compressed.100");
    cv::Mat image = DenseMatrix2Mat(newData);
    cv::namedWindow("Display window",
                    cv::WINDOW_AUTOSIZE);  // Create a window for display.
    cv::imshow("Display window", image);
    cv::waitKey(0);
    return 0;
}