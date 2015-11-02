#include <iostream>
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
	// Not yet done
	testComputeAngleBetween2D();
	testComputeAngleBetween3D();
	testComputeDistancePointToLine();
	testComputeDistancePointToPlane();
	testComputeDistancePointToSphere();
	testComputeDistanceCircleToCircle();
	testComputeDistanceCircleToTriangle();
	// Not yet done
	testComputeLineSphereIntersection();
	testIsLinePlaneIntersection();
	testIsIntersection();

	std::cout << "\nPress ENTER/RETURN to exit";
	std::cin.get();

	return 0;
}