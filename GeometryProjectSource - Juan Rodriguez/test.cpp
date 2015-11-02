#include "test.h"

using std::cout;

void testEquals() {

	TVector3 vec1;
	TVector3 vec2;

	vec1.m_fX = 4;
	vec1.m_fY = 7;
	vec1.m_fZ = 3;

	vec2.m_fX = 4;
	vec2.m_fY = 7;
	vec2.m_fZ = 3;

	cout << "testEquals ";
	if (Equals(vec1, vec2) == true) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testAdd() {

	TVector3 vec1;
	TVector3 vec2;
	TVector3 res;

	vec1.m_fX = 1;
	vec1.m_fY = -5;
	vec1.m_fZ = 8;

	vec2.m_fX = 4;
	vec2.m_fY = 1;
	vec2.m_fZ = 90.0f;

	res.m_fX = vec1.m_fX + vec2.m_fX;
	res.m_fY = vec1.m_fY + vec2.m_fY;
	res.m_fZ = vec1.m_fZ + vec2.m_fZ;

	cout << "testAdd ";
	if ((res.m_fX == Add(vec1, vec2, res).m_fX) &&
		(res.m_fY == Add(vec1, vec2, res).m_fY) &&
		(res.m_fZ == Add(vec1, vec2, res).m_fZ)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testSubstract() {

	TVector3 vec1;
	TVector3 vec2;
	TVector3 res;

	vec1.m_fX = 3;
	vec1.m_fY = 7;
	vec1.m_fZ = -4.15f;

	vec2.m_fX = 5.18f;
	vec2.m_fY = -2;
	vec2.m_fZ = 0;

	res.m_fX = vec1.m_fX - vec2.m_fX;
	res.m_fY = vec1.m_fY - vec2.m_fY;
	res.m_fZ = vec1.m_fZ - vec2.m_fZ;

	cout << "testSubtract ";
	if ((res.m_fX == Subtract(vec1, vec2, res).m_fX) &&
		(res.m_fY == Subtract(vec1, vec2, res).m_fY) &&
		(res.m_fZ == Subtract(vec1, vec2, res).m_fZ)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testScaleVector() {

	TVector3 vec1;
	float scalar;
	TVector3 res;

	vec1.m_fX = 3.14f;
	vec1.m_fY = 7.00f;
	vec1.m_fZ = -43;

	scalar = 5.13f;

	res.m_fX = vec1.m_fX * scalar;
	res.m_fY = vec1.m_fY * scalar;
	res.m_fZ = vec1.m_fZ * scalar;

	cout << "testScaleVector ";
	if ((res.m_fX == ScaleVector(vec1, scalar, res).m_fX) &&
		(res.m_fY == ScaleVector(vec1, scalar, res).m_fY) &&
		(res.m_fZ == ScaleVector(vec1, scalar, res).m_fZ)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testMagnitude() {

	TVector3 vec1;
	float mag;

	vec1.m_fX = 9;
	vec1.m_fY = 0;
	vec1.m_fZ = 5.14f;

	mag = sqrt(pow(vec1.m_fX, 2) + pow(vec1.m_fY, 2) + pow(vec1.m_fZ, 2));

	cout << "testMagnitude ";
	if (mag == Magnitude(vec1)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testDotProduct() {

	TVector3 vec1;
	TVector3 vec2;
	float dotProduct;

	vec1.m_fX = -3.14f;
	vec1.m_fY = 7.75f;
	vec1.m_fZ = 8;

	vec2.m_fX = 6;
	vec2.m_fY = -2.50f;
	vec2.m_fZ = 0.0f;

	dotProduct = (vec1.m_fX * vec2.m_fX) + (vec1.m_fY * vec2.m_fY) + (vec1.m_fZ * vec2.m_fZ);

	cout << "testDotProduct ";
	if (dotProduct == DotProduct(vec1, vec2)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testCrossProduct() {

	TVector3 vec1;
	TVector3 vec2;
	TVector3 res;

	vec1.m_fX = -3.14f;
	vec1.m_fY = 7.75f;
	vec1.m_fZ = 8;

	vec2.m_fX = 6;
	vec2.m_fY = -2.50f;
	vec2.m_fZ = 0.0f;

	res.m_fX = (vec1.m_fY * vec2.m_fZ) - (vec2.m_fY * vec1.m_fZ);
	res.m_fY = ((vec1.m_fX * vec2.m_fZ) - (vec2.m_fX * vec1.m_fZ)) * -1;
	res.m_fX = (vec1.m_fX * vec2.m_fY) - (vec2.m_fX * vec1.m_fY);

	cout << "testCrossProduct ";
	if ((res.m_fX == CrossProduct(vec1, vec2, res).m_fX) &&
		(res.m_fY == CrossProduct(vec1, vec2, res).m_fY) &&
		(res.m_fZ == CrossProduct(vec1, vec2, res).m_fZ)) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testNormalise() {

	TVector3 vec1;
	TVector3 res;
	float mag;

	vec1.m_fX = 0.15f;
	vec1.m_fY = 6;
	vec1.m_fZ = 3.14f;

	mag = Magnitude(vec1);

	if (mag != 0.0f) {

		res.m_fX = vec1.m_fX / mag;
		res.m_fY = vec1.m_fY / mag;
		res.m_fZ = vec1.m_fZ / mag;

		cout << "testNormalise ";
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
		cout << "testNormalise SUCCESS\n";
	}

}

void testProjection() {

	TVector3 vec1;
	TVector3 vec2;
	TVector3 res;
	float magSq;

	vec1.m_fX = 2;
	vec1.m_fY = 4;
	vec1.m_fZ = 5;

	magSq = pow(Magnitude(vec2), 2);

	if (magSq != 0.0f) {

		res.m_fX = vec2.m_fX * (DotProduct(vec1, vec2)) / magSq;
		res.m_fY = vec2.m_fY * (DotProduct(vec1, vec2)) / magSq;
		res.m_fZ = vec2.m_fZ * (DotProduct(vec1, vec2)) / magSq;

		cout << "testProjection ";
		if ((res.m_fX == Projection(vec1, vec2, res).m_fX) &&
			(res.m_fY == Projection(vec1, vec2, res).m_fY) &&
			(res.m_fZ == Projection(vec1, vec2, res).m_fZ)) {
			cout << "SUCCESS\n";
		}
		else {
			cout << "FAILURE\n";
		}

	}
	else {
		cout << "testProjection SUCCESS\n";
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

	cout << "testCopmuteAngleBetween2D ";
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

	cout << "testCopmuteAngleBetween2D ";
	if (acos(DP / (magV1 * magV2) == (ComputeAngleBetween(vec1, vec2)))) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testComputeDistancePointToLine() {



}

void testComputeDistancePointToPlane() {



}

void testComputeDistancePointToSphere() {



}

void testComputeDistanceCircleToCircle() {



}

void testComputeDistanceCircleToTriangle() {



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

	cout << "testComputeLineSphereIntersection ";
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

	cout << "testIsLinePlaneIntersection ";
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

	cout << "testIsIntersection ";
	if ((IsIntersection(line1, line2) || ((Equals(line1.m_v3q, line2.m_v3q)) && (Equals(line1.m_v3v, line2.m_v3v)))) == true) {
		cout << "SUCCESS\n";
	}
	else {
		cout << "FAILURE\n";
	}

}

void testComputeIntersectionBetweenLines() {



}

void testIsInFieldOfView() {



}

void testIsSurfaceLit() {



}

void testFindTriangleNormal() {



}

void testRotateTriangleAroundPoint() {



}