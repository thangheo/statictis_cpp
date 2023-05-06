#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

int main()
{
    Vector3d v1(1, 2, 3);
    Vector3d v2(4, 5, 6);

    double dot_product = v1.dot(v2);

    std::cout << "Dot product of v1 and v2: " << dot_product << std::endl;

    return 0;
}
