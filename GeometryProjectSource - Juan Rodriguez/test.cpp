#include "test.h"

using std::cout;

void testEquals() {

	TVector3 vec1;
	TVector3 vec2;

	vec1 = {
		4.00000001f, 7, 3
	};

	vec2 = {
		4, 7, 3.0000f
	};

	cout << "> testEquals ";
	if (Equals(vec1, vec2) == true) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testAdd() {

	TVector3 resultant;
	TVector3 vec1;
	TVector3 vec2;

	vec1 = {
		1, -5.0000000001534f, 8
	};

	vec2 = {
		4, 1, 90.0f
	};

	resultant.m_fX = vec1.m_fX + vec2.m_fX;
	resultant.m_fY = vec1.m_fY + vec2.m_fY;
	resultant.m_fZ = vec1.m_fZ + vec2.m_fZ;

	cout << "> testAdd ";
	if ((resultant.m_fX == Add(vec1, vec2, resultant).m_fX) &&
		(resultant.m_fY == Add(vec1, vec2, resultant).m_fY) &&
		(resultant.m_fZ == Add(vec1, vec2, resultant).m_fZ)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testSubstract() {

	TVector3 resultant;
	TVector3 vec1;
	TVector3 vec2;

	vec1 = {
		3, 7, -4.15f
	};

	vec2 = {
		5.18f, -2, 0.00000000000000f
	};

	resultant.m_fX = vec1.m_fX - vec2.m_fX;
	resultant.m_fY = vec1.m_fY - vec2.m_fY;
	resultant.m_fZ = vec1.m_fZ - vec2.m_fZ;

	cout << "> testSubtract ";
	if ((resultant.m_fX == Subtract(vec1, vec2, resultant).m_fX) &&
		(resultant.m_fY == Subtract(vec1, vec2, resultant).m_fY) &&
		(resultant.m_fZ == Subtract(vec1, vec2, resultant).m_fZ)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testScaleVector() {

	TVector3 resultant;
	TVector3 vector;
	float scalar;

	vector = {
		3.14f, 7, -43.0012542f
	};

	scalar = 5.1353400001f;

	resultant.m_fX = vector.m_fX * scalar;
	resultant.m_fY = vector.m_fY * scalar;
	resultant.m_fZ = vector.m_fZ * scalar;

	cout << "> testScaleVector ";
	if ((resultant.m_fX == ScaleVector(vector, scalar, resultant).m_fX) &&
		(resultant.m_fY == ScaleVector(vector, scalar, resultant).m_fY) &&
		(resultant.m_fZ == ScaleVector(vector, scalar, resultant).m_fZ)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testMagnitude() {

	TVector3 vector;
	float magnitude;

	vector = {
		9.04143f, 0.1012031f, 5.14f
	};

	magnitude = sqrt(pow(vector.m_fX, 2) + pow(vector.m_fY, 2) + pow(vector.m_fZ, 2));

	cout << "> testMagnitude ";
	if (magnitude == Magnitude(vector)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testDotProduct() {

	float dotProduct;
	TVector3 vec1;
	TVector3 vec2;

	vec1 = {
		-3.14f, 7.751231023f, 8
	};

	vec2 = {
		6, -2.500001534f, 0.000000001f
	};

	dotProduct = (vec1.m_fX * vec2.m_fX) + (vec1.m_fY * vec2.m_fY) + (vec1.m_fZ * vec2.m_fZ);

	cout << "> testDotProduct ";
	if (dotProduct == DotProduct(vec1, vec2)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testCrossProduct() {

	TVector3 resultant;
	TVector3 vec1;
	TVector3 vec2;

	vec1 = {
		34.35534f, 9, 7.751231023f
	};

	vec2 = {
		2.5731274f, 124.4930f, 5.00043434000001f
	};

	resultant.m_fX = (vec1.m_fY * vec2.m_fZ) - (vec2.m_fY * vec1.m_fZ);
	resultant.m_fY = ((vec1.m_fX * vec2.m_fZ) - (vec2.m_fX * vec1.m_fZ)) * -1;
	resultant.m_fX = (vec1.m_fX * vec2.m_fY) - (vec2.m_fX * vec1.m_fY);

	cout << "> testCrossProduct ";
	if ((resultant.m_fX == CrossProduct(vec1, vec2, resultant).m_fX) &&
		(resultant.m_fY == CrossProduct(vec1, vec2, resultant).m_fY) &&
		(resultant.m_fZ == CrossProduct(vec1, vec2, resultant).m_fZ)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testNormalise() {

	float magnitude;
	TVector3 vec1;
	TVector3 res;

	vec1 = {
		3, 2.00041243f, -4.049278388f
	};

	magnitude = Magnitude(vec1);

	if (magnitude != 0.0f) {

		res.m_fX = vec1.m_fX / magnitude;
		res.m_fY = vec1.m_fY / magnitude;
		res.m_fZ = vec1.m_fZ / magnitude;

		cout << "> testNormalise ";
		if ((res.m_fX == Normalise(vec1, res).m_fX) &&
			(res.m_fY == Normalise(vec1, res).m_fY) &&
			(res.m_fZ == Normalise(vec1, res).m_fZ)) {
			cout << "SUCCESS\n";
		}
		else {
			cout << "FAILURE\n";
		}

	}
	else {
		cout << "> testNormalise SUCCESS\n";
	}

}

void testProjection() {

	TVector3 resultant;
	TVector3 vec1;
	TVector3 vec2;
	float magSq;

	vec1 = {
		2.40004f, 6.24390f, 5
	};

	vec2 = {
		4.314f, 2, -3.14159f
	};

	magSq = pow(Magnitude(vec2), 2);

	if (magSq != 0.0f) {

		resultant.m_fX = vec2.m_fX * (DotProduct(vec1, vec2)) / magSq;
		resultant.m_fY = vec2.m_fY * (DotProduct(vec1, vec2)) / magSq;
		resultant.m_fZ = vec2.m_fZ * (DotProduct(vec1, vec2)) / magSq;

		cout << "> testProjection ";
		if ((resultant.m_fX == Projection(vec1, vec2, resultant).m_fX) &&
			(resultant.m_fY == Projection(vec1, vec2, resultant).m_fY) &&
			(resultant.m_fZ == Projection(vec1, vec2, resultant).m_fZ)) {
			cout << "SUCCESS\n";
		}
		else {
			cout << "FAILURE\n";
		}

	}
	else {
		cout << "> testProjection SUCCESS\n";
	}

}

void testComputeAngleBetween2D() {

	TVector2 vec1;
	TVector2 vec2;

	vec1.m_fX = 5;
	vec1.m_fY = -7;

	vec2.m_fX = 3;
	vec2.m_fY = 2;

	float DP = (vec1.m_fX * vec2.m_fX) + (vec1.m_fY * vec2.m_fY);
	float magV1 = sqrt(pow(vec1.m_fX, 2) + pow(vec1.m_fY, 2));
	float magV2 = sqrt(pow(vec2.m_fX, 2) + pow(vec2.m_fY, 2));

	cout << "> testCopmuteAngleBetween2D ";
	if (acos(DP / (magV1 * magV2) == (ComputeAngleBetween(vec1, vec2)))) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testComputeAngleBetween3D() {

	TVector3 vec1;
	TVector3 vec2;

	vec1.m_fX = 5;
	vec1.m_fY = -7;
	vec1.m_fZ = 4;

	vec2.m_fX = 3;
	vec2.m_fY = 2;
	vec2.m_fZ = 0;

	float DP = (vec1.m_fX * vec2.m_fX) + (vec1.m_fY * vec2.m_fY) + (vec1.m_fZ + vec2.m_fZ);
	float magV1 = sqrt(pow(vec1.m_fX, 2) + pow(vec1.m_fY, 2) + pow(vec1.m_fZ, 2));
	float magV2 = sqrt(pow(vec2.m_fX, 2) + pow(vec2.m_fY, 2) + pow(vec2.m_fZ, 2));

	cout << "> testCopmuteAngleBetween2D ";
	if (acos(DP / (magV1 * magV2) == (ComputeAngleBetween(vec1, vec2)))) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testComputeDistancePointToLine() {

	TVector3 vector;
	TVector3 crossP;
	TVector3 point;
	float distance;
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
	distance = Magnitude(crossP) / Magnitude(line.m_v3v);

	cout << "> testComputeDistancePointToLine ";
	if (distance == ComputeDistancePointToLine(line, point)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testComputeDistancePointToPlane() {

	TVector3 point;
	TPlane plane;

	TVector3 vector;
	TVector3 proj;

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

	/*float distance = (abs((plane.m_v3normal.m_fX * point.m_fX) + 
						 (plane.m_v3normal.m_fY * point.m_fY) +
						 (plane.m_v3normal.m_fZ * point.m_fZ)) + d) / 
					  Magnitude(plane.m_v3normal);*/

	cout << "> testComputeDistancePointToPlane ";
	if (ComputeDistancePointToPlane(plane, point) == Magnitude(proj)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testComputeDistancePointToSphere() {

	TVector3 point;
	TSphere sphere;

	// Sphere
	point = {
		3, 5, 8
	};

	// Sphere
	sphere.m_fRadius = 5;
	sphere.m_v3center.m_fX = 1;
	sphere.m_v3center.m_fY = 3;
	sphere.m_v3center.m_fZ = 0;

	float SQRT = sqrt(pow(sphere.m_v3center.m_fX - point.m_fX, 2) +
		pow(sphere.m_v3center.m_fY - point.m_fY, 2) +
		pow(sphere.m_v3center.m_fZ - point.m_fZ, 2));

	float distance = abs(SQRT - sphere.m_fRadius);

	cout << "> testComputeDistancePointToSphere ";
	if (distance == ComputeDistancePointToSphere(sphere, point)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testComputeDistanceCircleToCircle() {

	TCircle circle1;
	TCircle circle2;

	circle1 = {
		3, 5, 6
	};

	circle2 = {
		1, 4, 3
	};

	float SQRT = sqrt(pow(circle1.m_v2center.m_fX - circle2.m_v2center.m_fX, 2) + pow(circle1.m_v2center.m_fY - circle2.m_v2center.m_fY, 2));
	float deltaR = circle2.m_fRadius - circle1.m_fRadius;
	float distance = SQRT - deltaR;

	cout << "> testComputeDistanceCircleToCircle ";
	if (distance == ComputeDistanceCircleToCircle(circle1, circle2)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testComputeDistanceCircleToTriangle() {

	TCircle circle;
	TTriangle2 triangle;

	circle = {
		3, 5, 6
	};

	triangle = {
		3, 6
	};

	float triangleCenterX = (triangle.m_v2p1.m_fX + triangle.m_v2p2.m_fX + triangle.m_v2p3.m_fX) / 3;
	float triangleCenterY = (triangle.m_v2p1.m_fY + triangle.m_v2p2.m_fY + triangle.m_v2p3.m_fY) / 3;

	float SQRT = sqrt(pow(circle.m_v2center.m_fX - triangleCenterX, 2) + pow(circle.m_v2center.m_fY - triangleCenterY, 2));
	float distance = abs(SQRT - circle.m_fRadius);

	cout << "> testComputeDistanceCircleToTriangle ";
	if (distance == ComputeDistanceCircleToTriangle(circle, triangle)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testComputeLineSphereIntersection() {

	T3DLine line;
	TSphere sphere;
	TVector3 intersection1;
	TVector3 intersection2;

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

	cout << "> testComputeLineSphereIntersection ";
	if (ComputeLineSphereIntersection(line, sphere, intersection1, intersection2) == INTERSECTION_TWO) {
		cout << "SUCCESS - TWO INTERSECTIONS\n";
	}
	else if (ComputeLineSphereIntersection(line, sphere, intersection1, intersection2) == INTERSECTION_ONE) {
		cout << "SUCCESS - ONE INTERSECTION\n";
	}
	else {
		cout << "FAILURE - NO INTERSECTIONS\n";
	}

}

void testIsLinePlaneIntersection() {

	TPlane plane;
	T3DLine line;
	TVector3 intersection;

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

	cout << "> testIsLinePlaneIntersection ";
	if (IsLinePlaneIntersection(line, plane, intersection) == true) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
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

	cout << "> testIsIntersection ";
	if ((IsIntersection(line1, line2) || ((Equals(line1.m_v3q, line2.m_v3q)) && (Equals(line1.m_v3v, line2.m_v3v)))) == true) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testComputeIntersectionBetweenLines() {

	T3DLine line1;
	T3DLine line2;
	TVector3 intersection;

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

	intersection = ComputeIntersectionBetweenLines(line1, line2, intersection);

	cout << "> testComputeIntersectionBetweenLines ";
	if (intersection.m_fX && intersection.m_fY && intersection.m_fZ) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
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

	fieldOfViewInDegrees = 120;
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

	cout << "> testIsInFieldOfView ";
	if (IsInFieldOfView(cameraPosition, cameraDirection, fieldOfViewInRadians, objectPosition) == true) {
		cout << "SUCCESS - IS IN FIELD OF VIEW\n";
	}
	else {
		cout << "FAILURE - NOT IN FIELD OF VIEW\n";
	}

}

void testIsSurfaceLit() {



}

void testFindTriangleNormal() {

	// Triangle
	TTriangle3 triangle;
	// Point 1
	triangle.m_v3p1.m_fX = 8;
	triangle.m_v3p1.m_fY = 6.18384f;
	triangle.m_v3p1.m_fZ = 3;
	// Point 2
	triangle.m_v3p2.m_fX = 7.193840f;
	triangle.m_v3p2.m_fY = 4.34f;
	triangle.m_v3p2.m_fZ = 2;
	// Point 3
	triangle.m_v3p3.m_fX = 8.43f;
	triangle.m_v3p3.m_fY = 6/0.4242f;
	triangle.m_v3p3.m_fZ = 3;

	TVector3 N; // Normal
	TVector3 V; // P2 - P1
	TVector3 W; // P3 - P1

	// First point 
	TVector3 P1;
	P1.m_fX = triangle.m_v3p1.m_fX;
	P1.m_fY = triangle.m_v3p1.m_fY;
	P1.m_fZ = triangle.m_v3p1.m_fZ;

	// Second point
	TVector3 P2;
	P2.m_fX = triangle.m_v3p2.m_fX;
	P2.m_fY = triangle.m_v3p2.m_fY;
	P2.m_fZ = triangle.m_v3p2.m_fZ;

	// Third point
	TVector3 P3;
	P3.m_fX = triangle.m_v3p3.m_fX;
	P3.m_fY = triangle.m_v3p3.m_fY;
	P3.m_fZ = triangle.m_v3p3.m_fZ;

	// V = P2 - P1
	V.m_fX = Subtract(P2, P1, V).m_fX;
	V.m_fY = Subtract(P2, P1, V).m_fY;
	V.m_fZ = Subtract(P2, P1, V).m_fZ;

	// W = P3 - P1
	W.m_fX = Subtract(P3, P1, W).m_fX;
	W.m_fY = Subtract(P3, P1, W).m_fY;
	W.m_fZ = Subtract(P3, P1, W).m_fZ;

	// N = our normal
	N = CrossProduct(V, W, N);

	cout << "> testFindTriangleNormal ";
	if ((N.m_fX == FindTriangleNormal(triangle, N).m_fX) &&
		(N.m_fY == FindTriangleNormal(triangle, N).m_fY) && 
		(N.m_fZ == FindTriangleNormal(triangle, N).m_fZ)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testRotateTriangleAroundPoint() {

	float angleToRotateInRadians;
	TTriangle2 resultingTriangle;
	TVector2 rotationPoint;
	TTriangle2 triangle;

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
	resultingTriangle.m_v2p1.m_fX = x1;
	resultingTriangle.m_v2p1.m_fY = y1;

	// Second point
	resultingTriangle.m_v2p2.m_fX = x2;
	resultingTriangle.m_v2p2.m_fY = y2;

	// Third point
	resultingTriangle.m_v2p3.m_fX = x3;
	resultingTriangle.m_v2p3.m_fY = y3;

	// Since operations to check each point are too long, we'll compare each point individually
	bool isFirstPointRotated = 
		(resultingTriangle.m_v2p1.m_fX == RotateTriangleAroundPoint(triangle, angleToRotateInRadians, rotationPoint, resultingTriangle).m_v2p1.m_fX) &&
		(resultingTriangle.m_v2p1.m_fY == RotateTriangleAroundPoint(triangle, angleToRotateInRadians, rotationPoint, resultingTriangle).m_v2p1.m_fY);

	bool isSecondPointRotated =
		(resultingTriangle.m_v2p2.m_fX == RotateTriangleAroundPoint(triangle, angleToRotateInRadians, rotationPoint, resultingTriangle).m_v2p2.m_fX) &&
		(resultingTriangle.m_v2p2.m_fY == RotateTriangleAroundPoint(triangle, angleToRotateInRadians, rotationPoint, resultingTriangle).m_v2p2.m_fY);

	bool isThirdPointRotated =
		(resultingTriangle.m_v2p3.m_fX == RotateTriangleAroundPoint(triangle, angleToRotateInRadians, rotationPoint, resultingTriangle).m_v2p3.m_fX) &&
		(resultingTriangle.m_v2p3.m_fY == RotateTriangleAroundPoint(triangle, angleToRotateInRadians, rotationPoint, resultingTriangle).m_v2p3.m_fY);

	cout << "> testRotateTriangleAroundPoint ";
	if (isFirstPointRotated && isSecondPointRotated && isThirdPointRotated) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}