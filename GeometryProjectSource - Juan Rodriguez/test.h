#define _USE_MATH_DEFINES
#include <iostream>
#include "geometry.h"
#include "windows.h"

using std::cout;
using std::endl;

// Colors, just for fun!
enum Color {
	BLACK = 0,
	DARKBLUE = 1,
	DARKGREEN = 2,
	DARKTEAL = 3,
	DARKRED = 4,
	DARKPINK = 5,
	DARKYELLOW = 6,
	GRAY = 7,
	DARKGRAY = 8,
	BLUE = 9,
	GREEN = 10,
	TEAL = 11,
	RED = 12,
	PINK = 13,
	YELLOW = 14,
	WHITE = 15
};

void setColor(Color c);

void testEquals();

void testAdd();

void testSubstract();

void testScaleVector();

void testMagnitude();

void testDotProduct();

void testCrossProduct();

void testNormalise();

void testProjection();

void testComputeAngleBetween2D();

void testComputeAngleBetween3D();

void testComputeDistancePointToLine();

void testComputeDistancePointToPlane();

void testComputeDistancePointToSphere();

void testComputeDistanceCircleToCircle();

void testComputeDistanceCircleToTriangle();

void testComputeLineSphereIntersection();

void testIsLinePlaneIntersection();

void testIsIntersection();

void testComputeIntersectionBetweenLines();

void testIsInFieldOfView();

void testIsSurfaceLit();

void testFindTriangleNormal();

void testRotateTriangleAroundPoint();