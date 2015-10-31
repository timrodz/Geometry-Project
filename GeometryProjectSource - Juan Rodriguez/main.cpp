#include <iostream>

//#include "geometry.h"
#include "test.h"

int main() {

	testEquals();
	testAdd();
	testSubstract();
	testScaleVector();
	testMagnitude();
	testDotProduct();
	testCrossProduct();
	testNormalise();
	testProjection();
	testComputeAngleBetween2D();
	testComputeAngleBetween3D();
	testComputeDistancePointToLine();
	testComputeDistancePointToPlane();
	testComputeDistancePointToSphere();
	testComputeDistanceCircleToCircle();
	testComputeDistanceCircleToTriangle();
	testComputeLineSphereIntersection();

	std::cout << "\nPress ENTER/RETURN to exit";
	std::cin.get();

	return 0;
}