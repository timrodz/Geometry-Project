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
	testComputeAngleBetween2D();
	testComputeAngleBetween3D();
	testComputeDistancePointToLine();
	testComputeDistancePointToPlane();
	testComputeDistancePointToSphere();
	testComputeDistanceCircleToCircle();
	testComputeDistanceCircleToTriangle();
	testComputeLineSphereIntersection();
	testIsLinePlaneIntersection();
	testIsIntersection();
	testComputeIntersectionBetweenLines();
	testIsInFieldOfView();
	testIsSurfaceLit();
	testFindTriangleNormal();
	testRotateTriangleAroundPoint();

	std::cout << "Press ENTER/RETURN to exit";
	std::cin.get();

	return 0;

}