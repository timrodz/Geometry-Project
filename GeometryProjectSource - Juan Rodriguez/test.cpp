#include "test.h"

void setColor(Color c) {
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), c);
}

void testEquals() {

	TVector3 V1;
	TVector3 V2;
	bool result;

	V1 = {
		4.000001f, 7, 3
	};

	V2 = {
		4, 7, 3.0000f
	};

	result = Equals(V1, V2);

	setColor(WHITE);
	cout << "> test: Equals: ";
	if (result == true) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: " << std::boolalpha << true << ", instead got: " << std::boolalpha << result << "\n";

	}

}

void testAdd() {

	TVector3 expectedResult;
	TVector3 result;
	TVector3 V1;
	TVector3 V2;

	V1 = {
		1, -5.0000000001534f, 8
	};

	V2 = {
		4, 1, 90.0f
	};

	result.m_fX = V1.m_fX + V2.m_fX;
	result.m_fY = V1.m_fY + V2.m_fY;
	result.m_fZ = V1.m_fZ + V2.m_fZ;

	expectedResult = Add(V1, V2, expectedResult);

	setColor(WHITE);
	cout << "> test: Add: ";
	if ((result.m_fX == expectedResult.m_fX) &&
		(result.m_fY == expectedResult.m_fY) &&
		(result.m_fZ == expectedResult.m_fZ)) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: x:" << expectedResult.m_fX << ", y: " << expectedResult.m_fY << ", z: " << expectedResult.m_fZ << ".\n";
		cout << "Instead got : x:" << result.m_fX << ", y: " << result.m_fY << ", z: " << result.m_fZ << ".\n";

	}

}

void testSubstract() {

	TVector3 expectedResult;
	TVector3 result;
	TVector3 V1;
	TVector3 V2;

	V1 = {
		3, 7, -4.15f
	};

	V2 = {
		5.18f, -2, 0.00000000000000f
	};

	result.m_fX = V1.m_fX - V2.m_fX;
	result.m_fY = V1.m_fY - V2.m_fY;
	result.m_fZ = V1.m_fZ - V2.m_fZ;

	expectedResult = Subtract(V1, V2, expectedResult);

	setColor(WHITE);
	cout << "> test: Subtract: ";
	if ((result.m_fX == expectedResult.m_fX) &&
		(result.m_fY == expectedResult.m_fY) &&
		(result.m_fZ == expectedResult.m_fZ)) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: x:" << expectedResult.m_fX << ", y: " << expectedResult.m_fY << ", z: " << expectedResult.m_fZ << ".\n";
		cout << "Instead got : x:" << result.m_fX << ", y: " << result.m_fY << ", z: " << result.m_fZ << ".\n";

	}

}

void testScaleVector() {

	TVector3 expectedResult;
	TVector3 result;
	TVector3 vector;
	float scalar;

	vector = {
		3.14f, 7, -43.0012542f
	};

	scalar = 5.1353400001f;

	result.m_fX = vector.m_fX * scalar;
	result.m_fY = vector.m_fY * scalar;
	result.m_fZ = vector.m_fZ * scalar;

	expectedResult = ScaleVector(vector, scalar, expectedResult);

	setColor(WHITE);
	cout << "> test: ScaleVector: ";
	if ((result.m_fX == expectedResult.m_fX) &&
		(result.m_fY == expectedResult.m_fY) &&
		(result.m_fZ == expectedResult.m_fZ)) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: x:" << expectedResult.m_fX << ", y: " << expectedResult.m_fY << ", z: " << expectedResult.m_fZ << ".\n";
		cout << "Instead got : x:" << result.m_fX << ", y: " << result.m_fY << ", z: " << result.m_fZ << ".\n";

	}

}

void testMagnitude() {

	float expectedResult;
	TVector3 vector;
	float result;

	vector = {
		9.04143f, 0.1012031f, 5.14f
	};

	result = sqrt(pow(vector.m_fX, 2) + pow(vector.m_fY, 2) + pow(vector.m_fZ, 2));
	expectedResult = Magnitude(vector);

	setColor(WHITE);
	cout << "> test: Magnitude: ";
	if (result == expectedResult) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected result -> " << expectedResult << "\n";
		cout << "Actual result   -> " << result << "\n";

	}

}

void testDotProduct() {

	float expectedResult;
	TVector3 V1;
	TVector3 V2;
	float result;

	V1 = {
		-3.14f, 7.751231023f, 8
	};

	V2 = {
		6, -2.500001534f, 0.000000001f
	};

	result = (V1.m_fX * V2.m_fX) + (V1.m_fY * V2.m_fY) + (V1.m_fZ * V2.m_fZ);
	expectedResult = DotProduct(V1, V2);

	setColor(WHITE);
	cout << "> test: DotProduct: ";
	if (result == expectedResult) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: " << expectedResult << ", instead got: " << result;

	}

}

void testCrossProduct() {

	TVector3 expectedResult;
	TVector3 result;
	TVector3 V1;
	TVector3 V2;

	V1 = {
		34.35534f, 9, 7.751231023f
	};

	V2 = {
		2.5731274f, 124.4930f, 5.00043434000001f
	};

	result.m_fX = (V1.m_fY * V2.m_fZ) - (V2.m_fY * V1.m_fZ);
	result.m_fY = ((V1.m_fX * V2.m_fZ) - (V2.m_fX * V1.m_fZ)) * -1;
	result.m_fZ = (V1.m_fX * V2.m_fY) - (V2.m_fX * V1.m_fY);

	expectedResult = CrossProduct(V1, V2, expectedResult);

	setColor(WHITE);
	cout << "> test: CrossProduct: ";
	if ((result.m_fX == expectedResult.m_fX) &&
		(result.m_fY == expectedResult.m_fY) &&
		(result.m_fZ == expectedResult.m_fZ)) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: x:" << expectedResult.m_fX << ", y: " << expectedResult.m_fY << ", z: " << expectedResult.m_fZ << ".\n";
		cout << "Instead got : x:" << result.m_fX << ", y: " << result.m_fY << ", z: " << result.m_fZ << ".\n";

	}

}

void testNormalise() {

	TVector3 expectedResult;
	float magnitude;
	TVector3 V1;
	TVector3 result;

	V1 = {
		3, 2.00041243f, -4.049278388f
	};

	magnitude = Magnitude(V1);

	setColor(WHITE);
	cout << "> test: Normalise: ";
	if (magnitude != 0.0f) {

		result.m_fX = V1.m_fX / magnitude;
		result.m_fY = V1.m_fY / magnitude;
		result.m_fZ = V1.m_fZ / magnitude;

		expectedResult = Normalise(V1, expectedResult); \

			if ((result.m_fX == expectedResult.m_fX) &&
				(result.m_fY == expectedResult.m_fY) &&
				(result.m_fZ == expectedResult.m_fZ)) {

				setColor(GREEN);
				cout << "SUCCESS\n";

			}
			else {

				setColor(RED);
				cout << "FAILURE\n";
				setColor(WHITE);
				cout << "Expected: x:" << expectedResult.m_fX << ", y: " << expectedResult.m_fY << ", z: " << expectedResult.m_fZ << ".\n";
				cout << "Instead got : x:" << result.m_fX << ", y: " << result.m_fY << ", z: " << result.m_fZ << ".\n";

			}

	}
	else {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}

}

void testProjection() {

	TVector3 expectedResult;
	TVector3 result;
	TVector3 V1;
	TVector3 V2;
	float magSq;

	V1 = {
		2.40004f, 6.24390f, 5
	};

	V2 = {
		4.314f, 2, -3.14159f
	};

	magSq = pow(Magnitude(V2), 2);

	setColor(WHITE);
	cout << "> test: Projection: ";
	if (magSq != 0.0f) {

		result.m_fX = V2.m_fX * (DotProduct(V1, V2)) / magSq;
		result.m_fY = V2.m_fY * (DotProduct(V1, V2)) / magSq;
		result.m_fZ = V2.m_fZ * (DotProduct(V1, V2)) / magSq;

		expectedResult = Projection(V1, V2, expectedResult);

		if ((result.m_fX == expectedResult.m_fX) &&
			(result.m_fY == expectedResult.m_fY) &&
			(result.m_fZ == expectedResult.m_fZ)) {

			setColor(GREEN);
			cout << "SUCCESS\n";

		}
		else {

			setColor(RED);
			cout << "FAILURE\n";
			setColor(WHITE);
			cout << "Expected: x:" << expectedResult.m_fX << ", y: " << expectedResult.m_fY << ", z: " << expectedResult.m_fZ << ".\n";
			cout << "Instead got : x:" << result.m_fX << ", y: " << result.m_fY << ", z: " << result.m_fZ << ".\n";

		}

	}
	else {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}

}

void testComputeAngleBetween2D() {

	float expectedResult;
	float magnitudeV1;
	float magnitudeV2;
	float dotProduct;
	float result;
	TVector2 V1;
	TVector2 V2;

	V1.m_fX = 5;
	V1.m_fY = -7;

	V2.m_fX = 3;
	V2.m_fY = 2;

	dotProduct = (V1.m_fX * V2.m_fX) + (V1.m_fY * V2.m_fY);
	magnitudeV1 = sqrt(pow(V1.m_fX, 2) + pow(V1.m_fY, 2));
	magnitudeV2 = sqrt(pow(V2.m_fX, 2) + pow(V2.m_fY, 2));

	result = acos(dotProduct / (magnitudeV1 * magnitudeV2));

	expectedResult = ComputeAngleBetween(V1, V2);

	setColor(WHITE);
	cout << "> test: CopmuteAngleBetween2D: ";
	if (result == expectedResult) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: " << expectedResult << ", instead got: " << result;

	}

}

void testComputeAngleBetween3D() {

	float expectedResult;
	float result;
	TVector3 V1;
	TVector3 V2;

	V1.m_fX = 5;
	V1.m_fY = -7;
	V1.m_fZ = 4;

	V2.m_fX = 3;
	V2.m_fY = 2;
	V2.m_fZ = 0.0001f;

	result = acosf(DotProduct(V1, V2) / (Magnitude(V1) * Magnitude(V2)));

	expectedResult = ComputeAngleBetween(V1, V2);

	setColor(WHITE);
	cout << "> test: CopmuteAngleBetween3D: ";
	if (result == expectedResult) {

		setColor(GREEN);
		cout << "SUCCESS\n";
	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: " << expectedResult << ", instead got: " << result;

	}

}

void testComputeDistancePointToLine() {

	float expectedResult;
	TVector3 vector;
	TVector3 crossP;
	TVector3 point;
	float result;
	T3DLine line;

	// Point
	point.m_fX = 5;
	point.m_fY = 3;
	point.m_fZ = 0;

	// Line
	line.m_v3q.m_fX = 0;
	line.m_v3q.m_fY = 2;
	line.m_v3q.m_fZ = 0;

	line.m_v3v.m_fX = 1;
	line.m_v3v.m_fY = 3;
	line.m_v3v.m_fZ = 1;

	vector = Subtract(point, line.m_v3q, vector);
	crossP = CrossProduct(vector, line.m_v3v, crossP);

	result = Magnitude(crossP) / Magnitude(line.m_v3v);

	expectedResult = ComputeDistancePointToLine(line, point);

	setColor(WHITE);
	cout << "> test: ComputeDistancePointToLine: ";
	if (result == expectedResult) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: " << expectedResult << ", instead got: " << result;

	}

}

void testComputeDistancePointToPlane() {

	float expectedResult;
	TVector3 vector;
	TVector3 point;
	TVector3 proj;
	TPlane plane;
	float result;

	// Point
	point = {
		5, 3, 2
	};

	// Plane
	plane.m_v3point.m_fX = 1;
	plane.m_v3point.m_fY = 3;
	plane.m_v3point.m_fZ = 5;

	plane.m_v3normal.m_fX = 4;
	plane.m_v3normal.m_fY = 1;
	plane.m_v3normal.m_fZ = 2;

	/*float d = ((plane.m_v3normal.m_fX * plane.m_v3point.m_fX) +
	(plane.m_v3normal.m_fY * plane.m_v3point.m_fY) +
	(plane.m_v3normal.m_fZ * plane.m_v3point.m_fZ));*/

	vector = Subtract(point, plane.m_v3point, vector);
	proj = Projection(vector, plane.m_v3normal, proj);

	result = Magnitude(proj);

	expectedResult = ComputeDistancePointToPlane(plane, point);

	/*float distance = (abs((plane.m_v3normal.m_fX * point.m_fX) +
	(plane.m_v3normal.m_fY * point.m_fY) +
	(plane.m_v3normal.m_fZ * point.m_fZ)) + d) /
	Magnitude(plane.m_v3normal);*/

	setColor(WHITE);
	cout << "> test: ComputeDistancePointToPlane: ";
	if (result == expectedResult) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: " << expectedResult << ", instead got: " << result;

	}

}

void testComputeDistancePointToSphere() {

	float expectedResult;
	TVector3 point;
	TSphere sphere;
	float result;
	float SQRT;

	// Sphere
	point = {
		3, 5, 8
	};

	// Sphere
	sphere.m_fRadius = 5;
	sphere.m_v3center.m_fX = 1;
	sphere.m_v3center.m_fY = 3;
	sphere.m_v3center.m_fZ = 0;

	// square root of all the components squared
	SQRT = sqrt(pow(sphere.m_v3center.m_fX - point.m_fX, 2) +
		pow(sphere.m_v3center.m_fY - point.m_fY, 2) +
		pow(sphere.m_v3center.m_fZ - point.m_fZ, 2));

	result = abs(SQRT - sphere.m_fRadius);

	expectedResult = ComputeDistancePointToSphere(sphere, point);

	setColor(WHITE);
	cout << "> test: ComputeDistancePointToSphere: ";
	if (result == expectedResult) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: " << expectedResult << ", instead got: " << result;

	}

}

void testComputeDistanceCircleToCircle() {

	float expectedResult;
	TCircle circle1;
	TCircle circle2;
	float result;
	float deltaR;
	float SQRT;

	circle1 = {
		3, 5, 6
	};

	circle2 = {
		1, 4, 3
	};

	SQRT = sqrt(pow(circle1.m_v2center.m_fX - circle2.m_v2center.m_fX, 2) + pow(circle1.m_v2center.m_fY - circle2.m_v2center.m_fY, 2));
	deltaR = circle2.m_fRadius - circle1.m_fRadius;

	result = SQRT - deltaR;

	expectedResult = ComputeDistanceCircleToCircle(circle1, circle2);

	setColor(WHITE);
	cout << "> test: ComputeDistanceCircleToCircle: ";
	if (result == expectedResult) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: " << expectedResult << ", instead got: " << result;

	}

}

void testComputeDistanceCircleToTriangle() {

	float triangleCenterX;
	float triangleCenterY;
	float expectedResult;
	TTriangle2 triangle;
	TCircle circle;
	float result;
	float SQRT;

	circle = {
		3, 5, 6
	};

	triangle = {
		3, 6
	};

	triangleCenterX = (triangle.m_v2p1.m_fX + triangle.m_v2p2.m_fX + triangle.m_v2p3.m_fX) / 3;
	triangleCenterY = (triangle.m_v2p1.m_fY + triangle.m_v2p2.m_fY + triangle.m_v2p3.m_fY) / 3;

	SQRT = sqrt(pow(circle.m_v2center.m_fX - triangleCenterX, 2) + pow(circle.m_v2center.m_fY - triangleCenterY, 2));

	result = abs(SQRT - circle.m_fRadius);

	expectedResult = ComputeDistanceCircleToTriangle(circle, triangle);

	setColor(WHITE);
	cout << "> test: ComputeDistanceCircleToTriangle: ";
	if (result == expectedResult) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: " << expectedResult << ", instead got: " << result;

	}

}

void testComputeLineSphereIntersection() {

	EIntersections result;
	TVector3 intersection1;
	TVector3 intersection2;
	TSphere sphere;
	T3DLine line;

	line.m_v3q.m_fX = 4;       // x0
	line.m_v3q.m_fY = 5;       // y0
	line.m_v3q.m_fZ = 6;       // z0

	line.m_v3v.m_fX = 3;       // x1
	line.m_v3v.m_fY = -1;      // y1
	line.m_v3v.m_fZ = 5;       // z1

	sphere.m_fRadius = 5;       // r
	sphere.m_v3center.m_fX = 2; // h
	sphere.m_v3center.m_fY = 3; // k
	sphere.m_v3center.m_fZ = 1; // j

	result = ComputeLineSphereIntersection(line, sphere, intersection1, intersection2);

	setColor(WHITE);
	cout << "> test: ComputeLineSphereIntersection: ";
	if (result == INTERSECTION_TWO) {

		setColor(GREEN);
		cout << "SUCCESS - TWO INTERSECTIONS\n";

	}
	else if (result == INTERSECTION_ONE) {

		setColor(GREEN);
		cout << "SUCCESS - ONE INTERSECTION\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: " << std::boolalpha << true << ", instead got: " << std::boolalpha << result << "\n";

	}

}

void testIsLinePlaneIntersection() {

	TVector3 intersection;
	TPlane plane;
	T3DLine line;
	bool result;

	plane.m_v3normal.m_fX = 2;
	plane.m_v3normal.m_fY = 1;
	plane.m_v3normal.m_fZ = 4;

	plane.m_v3point.m_fX = 2;
	plane.m_v3point.m_fY = 5;
	plane.m_v3point.m_fZ = 6;

	line.m_v3q.m_fX = 0;
	line.m_v3q.m_fY = 2;
	line.m_v3q.m_fZ = 0;

	line.m_v3v.m_fX = 1;
	line.m_v3v.m_fY = 3;
	line.m_v3v.m_fZ = 1;

	result = IsLinePlaneIntersection(line, plane, intersection);

	setColor(WHITE);
	cout << "> test: IsLinePlaneIntersection ";
	if (result == true) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: " << std::boolalpha << true << ", instead got: " << std::boolalpha << result << "\n";

	}

}

void testIsIntersection() {

	T3DLine line1;
	T3DLine line2;

	// Line 1
	line1.m_v3q.m_fX = 3;
	line1.m_v3q.m_fY = -1;
	line1.m_v3q.m_fZ = 1;

	line1.m_v3v.m_fX = 1;
	line1.m_v3v.m_fY = 2;
	line1.m_v3v.m_fZ = 3;

	// Line 2
	line2.m_v3q.m_fX = 2;
	line2.m_v3q.m_fY = 5;
	line2.m_v3q.m_fZ = 0;

	line2.m_v3v.m_fX = 1;
	line2.m_v3v.m_fY = -1;
	line2.m_v3v.m_fZ = 1;

	bool result = IsIntersection(line1, line2);

	setColor(WHITE);
	cout << "> test: IsIntersection: ";
	if (result == true) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {
		
		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: " << std::boolalpha << true << ", instead got: " << std::boolalpha << result << "\n";

	}

}

void testComputeIntersectionBetweenLines() {

	TVector3 expectedResult;
	TVector3 result;
	T3DLine line1;
	T3DLine line2;

	// Line 1
	line1.m_v3q.m_fX = 3;
	line1.m_v3q.m_fY = 1;
	line1.m_v3q.m_fZ = -1;

	line1.m_v3v.m_fX = 1;
	line1.m_v3v.m_fY = 2;
	line1.m_v3v.m_fZ = 3;

	// Line 2
	line2.m_v3q.m_fX = 2;
	line2.m_v3q.m_fY = 5;
	line2.m_v3q.m_fZ = 0;

	line2.m_v3v.m_fX = 1;
	line2.m_v3v.m_fY = -1;
	line2.m_v3v.m_fZ = 1;

	// The value 't' for the formula: r(t) = q + t(v)
	float t = line1.m_v3q.m_fX + line1.m_v3v.m_fX - (line2.m_v3q.m_fX + line2.m_v3v.m_fX);

	result.m_fX = line1.m_v3q.m_fX + (t * line1.m_v3v.m_fX);
	result.m_fY = line1.m_v3q.m_fY + (t * line1.m_v3v.m_fY);
	result.m_fZ = line1.m_v3q.m_fZ + (t * line1.m_v3v.m_fZ);

	expectedResult = ComputeIntersectionBetweenLines(line1, line2, expectedResult);

	setColor(WHITE);
	cout << "> test: ComputeIntersectionBetweenLines: ";
	if ((result.m_fX == expectedResult.m_fX) && (result.m_fY == expectedResult.m_fY) && (result.m_fZ == expectedResult.m_fZ)) {
		
		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: x: " << expectedResult.m_fX << ", y: " << expectedResult.m_fY << ", z: " << expectedResult.m_fZ << ", ";
		cout << "instead got : x: " << result.m_fX << ", y: " << result.m_fY << ", z: " << result.m_fZ << ".\n";

	}

}

void testIsInFieldOfView() {

	// Our initial vectors
	float fieldOfViewInRadians;
	TVector2 cameraDirection;
	TVector2 cameraPosition;
	TVector2 objectPosition;

	// Values to check
	TVector2 vectorObjCameraPosNormalised;
	TVector2 cameraDirectionNormalised;
	TVector2 vectorObjCameraPos;

	float vectorObjCameraPosMagnitude;
	float cameraDirectionMagnitude;

	float normalisedVectorDotProduct;
	float fieldOfViewInDegrees;
	float angle;

	fieldOfViewInDegrees = 160;
	fieldOfViewInRadians = fieldOfViewInDegrees * float(M_PI) / 180.0f;

	cameraDirection = {
		1, 1
	};

	cameraPosition = {
		1, 3
	};

	objectPosition = {
		3, 2
	};

	// Vector from the camera position to the object
	vectorObjCameraPos.m_fX = objectPosition.m_fX - cameraPosition.m_fX;
	vectorObjCameraPos.m_fY = objectPosition.m_fY - cameraPosition.m_fY;

	// Magnitude
	cameraDirectionMagnitude = sqrt(pow(cameraDirection.m_fX, 2) + pow(cameraDirection.m_fX, 2));
	vectorObjCameraPosMagnitude = sqrt(pow(vectorObjCameraPos.m_fX, 2) + pow(vectorObjCameraPos.m_fY, 2));

	// Normalised vectors
	cameraDirectionNormalised.m_fX = cameraDirection.m_fX / cameraDirectionMagnitude;
	cameraDirectionNormalised.m_fY = cameraDirection.m_fY / cameraDirectionMagnitude;

	vectorObjCameraPosNormalised.m_fX = vectorObjCameraPos.m_fX / vectorObjCameraPosMagnitude;
	vectorObjCameraPosNormalised.m_fY = vectorObjCameraPos.m_fY / vectorObjCameraPosMagnitude;

	// Dot product of the normalised vectors
	normalisedVectorDotProduct = (cameraDirectionNormalised.m_fX * vectorObjCameraPosNormalised.m_fX) +
		(cameraDirectionNormalised.m_fY * vectorObjCameraPosNormalised.m_fY);

	angle = acos(normalisedVectorDotProduct);
	
	bool expectedResult = IsInFieldOfView(cameraPosition, cameraDirection, fieldOfViewInRadians, objectPosition);
	bool result = (angle < fieldOfViewInRadians / 2);

	setColor(WHITE);
	cout << "> test: IsInFieldOfView: ";
	if (result == expectedResult) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: " << std::boolalpha << expectedResult << ", instead got: " << std::boolalpha << result << "\n";

	}

}

void testIsSurfaceLit() {

	// Defining our points
	TVector3 pointOnSurface;
	TVector3 lightSource;
	TTriangle3 surface;

	pointOnSurface = {
		3, 5, 2
	};

	lightSource = {
		2, 7, 3
	};
	
	// P1
	surface.m_v3p1 = {
		5, 2, 4
	};

	// P2
	surface.m_v3p2 = {
		8, 3, 7
	};

	// P3
	surface.m_v3p3 = {
		2, 9, 6
	};

	// Knowing how much percentage of light has hit our surface
	float fPercentageOfLightingInSurface;
	// For storing the direction between the light and the surface normal
	float fAngleBetweenLightAndSurface;
	// For easier calculations
	TVector3 tNormalisedLightDirection;
	TVector3 tNormalisedSurfaceNormal;
	// Direction of the light
	TVector3 tLightDirection;
	// The surface normal
	TVector3 tSurfaceNormal;

	// Determining the normal from our surface
	tSurfaceNormal = FindTriangleNormal(surface, tSurfaceNormal);

	// Making a vector from the light source to the point on the surface
	tLightDirection = Subtract(pointOnSurface, lightSource, tLightDirection);

	tNormalisedLightDirection = Normalise(tLightDirection, tLightDirection);
	tNormalisedSurfaceNormal = Normalise(tSurfaceNormal, tSurfaceNormal);

	// Angle between the light direction vector and the surface normal, both normalised so we don't use magnitudes
	fAngleBetweenLightAndSurface = acosf(DotProduct(tNormalisedLightDirection, tNormalisedSurfaceNormal));

	// Knowing how much percent of lighting has actually struck our surface
	fPercentageOfLightingInSurface = cosf(fAngleBetweenLightAndSurface);

	bool expectedResult = IsSurfaceLit(pointOnSurface, lightSource, surface);
	bool result = fPercentageOfLightingInSurface > 0;

	setColor(WHITE);
	cout << "> test: IsSurfaceLit: ";
	if (result == expectedResult) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: " << std::boolalpha << expectedResult << ", instead got: " << std::boolalpha << result << "\n";

	}

}

void testFindTriangleNormal() {

	TVector3 expectedResult;
	TTriangle3 triangle;
	TVector3 result;
	TVector3 V; // P2 - P1
	TVector3 W; // P3 - P1

	// P1
	triangle.m_v3p1 = {
		8, 6.18384f, 3
	};
	// P2
	triangle.m_v3p2 = {
		7, 4, 2
	};
	// P3
	triangle.m_v3p3 = {
		8, 6 / 0.243f, 3
	};

	// V = P2 - P1
	V = Subtract(triangle.m_v3p2, triangle.m_v3p1, V);

	// W = P3 - P1
	W = Subtract(triangle.m_v3p3, triangle.m_v3p1, W);

	// N = our normal
	result = CrossProduct(V, W, result);

	expectedResult = FindTriangleNormal(triangle, expectedResult);

	setColor(WHITE);
	cout << "> test: FindTriangleNormal: ";
	if ((result.m_fX == expectedResult.m_fX) &&
		(result.m_fY == expectedResult.m_fY) &&
		(result.m_fZ == expectedResult.m_fZ)) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected: x:" << expectedResult.m_fX << ", y: " << expectedResult.m_fY << ", z: " << expectedResult.m_fZ << ".\n";
		cout << "Instead got : x:" << result.m_fX << ", y: " << result.m_fY << ", z: " << result.m_fZ << ".\n";

	}

}

void testRotateTriangleAroundPoint() {

	float angleToRotateInRadians;
	TTriangle2 expectedResult;
	TVector2 rotationPoint;
	TTriangle2 triangle;
	TTriangle2 result;

	// Angle to rotate by
	angleToRotateInRadians = 60 * float(M_PI) / 180.0f;

	// First point
	triangle.m_v2p1.m_fX = 2;
	triangle.m_v2p1.m_fY = 2;

	// Second point
	triangle.m_v2p2.m_fX = 6;
	triangle.m_v2p2.m_fY = 4;

	// Third point
	triangle.m_v2p3.m_fX = 6;
	triangle.m_v2p3.m_fY = 2;

	// Point to rotate
	rotationPoint.m_fX = 4;
	rotationPoint.m_fY = 3;

	// First point on the triangle
	float x1 = triangle.m_v2p1.m_fX;
	float y1 = triangle.m_v2p1.m_fY;

	// Second point on the triangle
	float x2 = triangle.m_v2p2.m_fX;
	float y2 = triangle.m_v2p2.m_fY;

	// Third point on the triangle
	float x3 = triangle.m_v2p3.m_fX;
	float y3 = triangle.m_v2p3.m_fY;

	// Rotation points
	float rx = rotationPoint.m_fX;
	float ry = rotationPoint.m_fY;

	// Translate every point to the center by subtracting our rotation points
	// X
	x1 -= rx;
	x2 -= rx;
	x3 -= rx;
	// Y
	y1 -= ry;
	y2 -= ry;
	y3 -= ry;

	// Rotating the first point
	rx = x1 * cosf(angleToRotateInRadians / 180.0f * float(M_PI)) - y1 * sinf(angleToRotateInRadians / 180.0f * float(M_PI));
	ry = y1 * cosf(angleToRotateInRadians / 180.0f * float(M_PI)) + x1 * sinf(angleToRotateInRadians / 180.0f * float(M_PI));
	x1 += rx;
	y1 += ry;

	// Rotating the second point
	rx = x2 * cosf(angleToRotateInRadians / 180.0f * float(M_PI)) - y2 * sinf(angleToRotateInRadians / 180.0f * float(M_PI));
	ry = y2 * cosf(angleToRotateInRadians / 180.0f * float(M_PI)) + x2 * sinf(angleToRotateInRadians / 180.0f * float(M_PI));
	x2 += rx;
	y2 += ry;

	// Rotating the third point
	rx = x3 * cosf(angleToRotateInRadians / 180.0f * float(M_PI)) - y3 * sinf(angleToRotateInRadians / 180.0f * float(M_PI));
	ry = y3 * cosf(angleToRotateInRadians / 180.0f * float(M_PI)) + x3 * sinf(angleToRotateInRadians / 180.0f * float(M_PI));
	x3 += rx;
	y3 += ry;

	// Asigning our rotated points to our final triangle
	// First point
	result.m_v2p1.m_fX = x1;
	result.m_v2p1.m_fY = y1;

	// Second point
	result.m_v2p2.m_fX = x2;
	result.m_v2p2.m_fY = y2;

	// Third point
	result.m_v2p3.m_fX = x3;
	result.m_v2p3.m_fY = y3;

	expectedResult = RotateTriangleAroundPoint(triangle, angleToRotateInRadians, rotationPoint, result);

	// Since operations to check each point are too long, we'll compare each point individually
	bool isFirstPointRotated = (result.m_v2p1.m_fX == expectedResult.m_v2p1.m_fX) &&
		(result.m_v2p1.m_fY == expectedResult.m_v2p1.m_fY);

	bool isSecondPointRotated = (result.m_v2p2.m_fX == expectedResult.m_v2p2.m_fX) &&
		(result.m_v2p2.m_fY == expectedResult.m_v2p2.m_fY);

	bool isThirdPointRotated = (result.m_v2p3.m_fX == expectedResult.m_v2p3.m_fX) &&
		(result.m_v2p3.m_fY == expectedResult.m_v2p3.m_fY);

	setColor(WHITE);
	cout << "> test: RotateTriangleAroundPoint: ";
	if ((isFirstPointRotated) && (isSecondPointRotated) && (isThirdPointRotated)) {

		setColor(GREEN);
		cout << "SUCCESS\n";

	}
	else {

		setColor(RED);
		cout << "FAILURE\n";
		setColor(WHITE);
		cout << "Expected:\t\t\tInstead got:\n";
		// P1
		cout << "P1: x: " << expectedResult.m_v2p1.m_fX << ", y: " << expectedResult.m_v2p1.m_fY << "\t";
		cout << "P1: x: " << result.m_v2p1.m_fX << ", y: " << result.m_v2p1.m_fY << "\n";
		// P2
		cout << "P2: x: " << expectedResult.m_v2p2.m_fX << ",  y: " << expectedResult.m_v2p2.m_fY << "\t";
		cout << "P2: x: " << result.m_v2p2.m_fX << ",  y: " << result.m_v2p2.m_fY << "\n";
		// P3
		cout << "P3: x: " << expectedResult.m_v2p3.m_fX << ",  y: " << expectedResult.m_v2p3.m_fY << "\t";
		cout << "P3: x: " << result.m_v2p3.m_fX << ",  y: " << result.m_v2p3.m_fY << "\n";

	}

}