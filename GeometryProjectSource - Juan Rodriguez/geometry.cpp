#define _USE_MATH_DEFINES
#include <cmath>
#include "geometry.h"

/***********************
* name of the function: Equals
* @author: Juan Rodriguez
* @parameter: _krA the first vector to test
* @parameter: _krB the second vector to test
* @return: boolean
********************/
bool Equals(const TVector3& _krA, const TVector3& _krB) {

	// Checks if every component of the vector is equal, if all XYZ are, return true
	if ((_krA.m_fX == _krB.m_fX) && (_krA.m_fY == _krB.m_fY) && (_krA.m_fZ == _krB.m_fZ)) {
		return true;
	}
	else {
		return false;
	}

}

/***********************
* name of the function: Add
* @author: Juan Rodriguez
* @parameter: _krA the first vector to add
* @parameter: _krB the second vector to add
* @parameter: _rResultant the resulting vector
* @return: _rResultant
********************/
TVector3& Add(const TVector3& _krA,
	const TVector3& _krB,
	TVector3& _rResultant) {

	_rResultant.m_fX = _krA.m_fX + _krB.m_fX;
	_rResultant.m_fY = _krA.m_fY + _krB.m_fY;
	_rResultant.m_fZ = _krA.m_fZ + _krB.m_fZ;

	return _rResultant;

}

/***********************
* name of the function: Subtract
* @author: Juan Rodriguez
* @parameter: _krA the first vector to subtract
* @parameter: _krB the second vector to subtract
* @parameter: _rResultant the resulting vector
* @return: _rResultant
********************/
TVector3& Subtract(const TVector3& _krA,
	const TVector3& _krB,
	TVector3& _rResultant) {

	_rResultant.m_fX = _krA.m_fX - _krB.m_fX;
	_rResultant.m_fY = _krA.m_fY - _krB.m_fY;
	_rResultant.m_fZ = _krA.m_fZ - _krB.m_fZ;

	return _rResultant;

}

/***********************
* name of the function: ScaleVector
* @author: Juan Rodriguez
* @parameter: _krA the vector
* @parameter: _kfScalar the scalar that multiplies _krA
* @parameter: _rResultant the resulting vector
* @return: _rResultant
********************/
TVector3& ScaleVector(const TVector3& _krA,
	const float _kfScalar,
	TVector3& _rResultant) {

	_rResultant.m_fX = _krA.m_fX * _kfScalar;
	_rResultant.m_fY = _krA.m_fY * _kfScalar;
	_rResultant.m_fZ = _krA.m_fZ * _kfScalar;

	return _rResultant;

}

/***********************
* name of the function: Magnitude
* @author: Juan Rodriguez
* @parameter: _krA the vector to compute its magnitude
* @return: square root of the vector's XYZ components squared
********************/
float Magnitude(const TVector3& _krA) {

	float magnitude = sqrt(pow(_krA.m_fX, 2) + pow(_krA.m_fY, 2) + pow(_krA.m_fZ, 2));
		
	return (magnitude);

}

/***********************
* name of the function: DotProdcut
* @author: Juan Rodriguez
* @parameter: _krA the first vector to multiply
* @parameter: _krB the second vector to multiply
* @return: The product of the multiplication of both vectors' XYZ components
********************/
float DotProduct(const TVector3& _krA, const TVector3& _krB) {

	float dotProduct = (_krA.m_fX * _krB.m_fX) + (_krA.m_fY * _krB.m_fY) + (_krA.m_fZ * _krB.m_fZ);

	return (dotProduct);

}

/***********************
* name of the function: CrossProduct
* @author: Juan Rodriguez
* @parameter: _krA the first vector
* @parameter: _krA the second vector
* @parameter: _rResultant the resulting vector
* @return: _rResultant
********************/
TVector3& CrossProduct(const TVector3& _krA,
	const TVector3& _krB,
	TVector3& _rResultant) {

	_rResultant.m_fX = (_krA.m_fY * _krB.m_fZ) - (_krB.m_fY * _krA.m_fZ);
	_rResultant.m_fY = ((_krA.m_fX * _krB.m_fZ) - (_krB.m_fX * _krA.m_fZ)) * -1.0f;
	_rResultant.m_fZ = (_krA.m_fX * _krB.m_fY) - (_krB.m_fX * _krA.m_fY);

	return _rResultant;
}

/***********************
* name of the function: Normalise
* @author: Juan Rodriguez
* @parameter: _krA the vector that will be normalised
* @parameter: _rResultant the normalised vector
* @return: _rResultant
********************/
TVector3& Normalise(const TVector3& _krA, TVector3& _rResultant) {

	_rResultant.m_fX = _krA.m_fX / Magnitude(_krA);
	_rResultant.m_fY = _krA.m_fY / Magnitude(_krA);
	_rResultant.m_fZ = _krA.m_fZ / Magnitude(_krA);

	return _rResultant;

}

/***********************
* name of the function: Projection
* @author: Juan Rodriguez
* @parameter: _krA the vector that will be projected
* @parameter: _krB the vector where _krA will be projected
* @parameter: _rResultant the resulting projected vector
* @return: _rResultant
********************/
TVector3& Projection(const TVector3& _krA,
	const TVector3& _krB,
	TVector3& _rResultant) {

	_rResultant.m_fX = (_krB.m_fX * DotProduct(_krA, _krB)) / pow(Magnitude(_krB), 2);
	_rResultant.m_fY = (_krB.m_fY * DotProduct(_krA, _krB)) / pow(Magnitude(_krB), 2);
	_rResultant.m_fZ = (_krB.m_fZ * DotProduct(_krA, _krB)) / pow(Magnitude(_krB), 2);

	return _rResultant;

}

/***********************
* name of the function: ComputeAngleBetween
* @author: Juan Rodriguez
* @parameter: _krA the first vector
* @parameter: _krB the second vector
* @return: angle between the two 2D Vectors
********************/
float ComputeAngleBetween(const TVector2& _krA,
	const TVector2& _krB) {

	float dotProduct = (_krA.m_fX * _krB.m_fX) + (_krA.m_fY * _krB.m_fY);
	float magnitude_krA = sqrt(pow(_krA.m_fX, 2) + pow(_krA.m_fY, 2));
	float magnitude_krB = sqrt(pow(_krB.m_fX, 2) + pow(_krB.m_fY, 2));

	return (
			acos(dotProduct / (magnitude_krA * magnitude_krB))
		);

}

/***********************
* name of the function: ComputeAngleBetween
* @author: Juan Rodriguez
* @parameter: _krA the first 3D vector
* @parameter: _krB the second 3D vector
* @return: angle between the two 3D Vectors
********************/
float ComputeAngleBetween(const TVector3& _krA,
	const TVector3& _krB) {

	// Since we've two 3D vectors, we recall the DotProduct and Magnitude functions that were previously created
	return (
			acos(DotProduct(_krA, _krB) / (Magnitude(_krA) * Magnitude(_krB)))
		);

}

/***********************
* name of the function: ComputeDistancePointToLine
* @author: Juan Rodriguez
* @parameter: _krLine the 3D line
* @parameter: _krPoint the 3D point
* @return: magnitude of the cross product between _krLine and _krPoint divided by the direction vector of _krLine
********************/
float ComputeDistancePointToLine(const T3DLine& _krLine,
	const TVector3& _krPoint) {

	TVector3 vector;
	TVector3 crossP;

	vector = Subtract(_krPoint, _krLine.m_v3q, vector);
	crossP = CrossProduct(vector, _krLine.m_v3v, crossP);

	return (
			Magnitude(crossP) / Magnitude(_krLine.m_v3v)
		);

}

/***********************
* name of the function: ComputeDistancePointToPlane
* @author: Juan Rodriguez
* @parameter: _krPlane the 3D plane
* @parameter: _krPoint the 3D point
* @return: distance
********************/
float ComputeDistancePointToPlane(const TPlane& _krPlane,
	const TVector3& _krPoint) {

	TVector3 vector;
	TVector3 proj;

	vector = Subtract(_krPoint, _krPlane.m_v3point, vector);
	proj = Projection(vector, _krPlane.m_v3normal, proj);

	float distance = Magnitude(proj);

	return (distance);

}

/***********************
* name of the function: ComputeDistancePointToSphere
* @author: Juan Rodriguez
* @parameter: _krSphere the 3D sphere
* @parameter: _krPoint the 3D point
* @return: the distance between the sphere minus its radius and the point
********************/
//Distance between point and center of the spheres
float ComputeDistancePointToSphere(const TSphere& _krSphere,
	const TVector3& _krPoint) {

	// Square root of the components except the radius
	float SQRT = sqrt(pow(_krSphere.m_v3center.m_fX - _krPoint.m_fX, 2) +
		pow(_krSphere.m_v3center.m_fY - _krPoint.m_fY, 2) +
		pow(_krSphere.m_v3center.m_fZ - _krPoint.m_fZ, 2));

	// We subtract the radius and, if it's a negative value, make it positive
	float distance = abs(SQRT - _krSphere.m_fRadius);

	return (distance);

}

/***********************
* name of the function: ComputeDistanceCircleToCircle
* @author: Juan Rodriguez
* @parameter: _krCircle1 the first circle
* @parameter: _krCircle2 the second circle
* @return: the distance between the two circles
********************/
float ComputeDistanceCircleToCircle(const TCircle& _krCircle1,
	const TCircle& _krCircle2) {

	// Square root of the components except the radius
	float distance = sqrt(pow(_krCircle1.m_v2center.m_fX - _krCircle2.m_v2center.m_fX, 2) +
		pow(_krCircle1.m_v2center.m_fY - _krCircle2.m_v2center.m_fY, 2));

	// The difference between the two radiuses
	float deltaR = _krCircle2.m_fRadius - _krCircle1.m_fRadius;

	return (distance - deltaR);

}

/***********************
* name of the function: ComputeDistanceCircleToTriangle
* @author: Juan Rodriguez
* @parameter: _krCircle the circle
* @parameter: _krTriangle the 2D triangle
* @return: distance between the circle minus the radius and the triangle
********************/
float ComputeDistanceCircleToTriangle(const TCircle& _krCircle,
	const TTriangle2& _krTriangle) {

	// Horizontal center of the triangle
	float triangleCentroidX = (_krTriangle.m_v2p1.m_fX + _krTriangle.m_v2p2.m_fX + _krTriangle.m_v2p3.m_fX) / 3;

	// Vertical center of the triangle
	float triangleCentroidY = (_krTriangle.m_v2p1.m_fY + _krTriangle.m_v2p2.m_fY + _krTriangle.m_v2p3.m_fY) / 3;

	// Square root of the center of the circle minus the center of the triangle
	float SQRT = sqrt(pow(_krCircle.m_v2center.m_fX - triangleCentroidX, 2) + 
				      pow(_krCircle.m_v2center.m_fY - triangleCentroidY, 2));

	float distance = abs(SQRT - _krCircle.m_fRadius);

	return (distance);

}

/***********************
* name of the function: ComputeDistanceCircleToTriangle
* @author: Juan Rodriguez
* @parameter: _krLine the 3D line
* @parameter: _krSphere the 3D sphere
* @parameter: _rv3IntersectionPoint1 the first point of intersection
* @parameter: _rv3IntersectionPoint2 the second point of intersection
* @return: EIntersections (two, one, none)
********************/
EIntersections ComputeLineSphereIntersection(const T3DLine& _krLine,
	const TSphere& _krSphere,
	TVector3& _rv3IntersectionPoint1,
	TVector3& _rv3IntersectionPoint2) {

	// Sphere
	float h = _krSphere.m_v3center.m_fX;
	float k = _krSphere.m_v3center.m_fY;
	float j = _krSphere.m_v3center.m_fZ;
	float r = _krSphere.m_fRadius;

	// Line
	float x0 = _krLine.m_v3q.m_fX;
	float y0 = _krLine.m_v3q.m_fY;
	float z0 = _krLine.m_v3q.m_fZ;
	
	// Saving time for our x1 - x0 calculations
	float x1_minus_x0 = x0 + (_krLine.m_v3v.m_fX - x0);
	float y1_minus_y0 = y0 + (_krLine.m_v3v.m_fY - y0);
	float z1_minus_z0 = z0 + (_krLine.m_v3v.m_fZ - z0);

	//    a = (x1 - x0) ^ 2       + (y1 - y0) ^ 2       + (z1 - z0) ^ 2
	float a = pow(x1_minus_x0, 2) + pow(y1_minus_y0, 2) + pow(z1_minus_z0, 2);

	float b = ((-2 * h)*(x1_minus_x0)+(2 * x0)*(x1_minus_x0)) -
		((-2 * k)*(y1_minus_y0)+(2 * y0)*(y1_minus_y0)) -
		((-2 * j)*(z1_minus_z0)+(2 * z0)*(z1_minus_z0));

	float c = pow(h, 2) + pow(k, 2) + pow(j, 2) +
		(-2 * h * x0) + pow(x0, 2) +
		(-2 * k * y0) + pow(y0, 2) +
		(-2 * j * z0) + pow(z0, 2) -
		pow(r, 2);

	float discriminant = sqrt(pow(b, 2) - (4.0f * a * c));

	if (discriminant > 0) {

		float t0 = (((-1.0f * b) + discriminant) / (2.0f * a));
		float t1 = (((-1.0f * b) - discriminant) / (2.0f * a));

		_rv3IntersectionPoint1.m_fX = t0 * (x1_minus_x0)+x0;
		_rv3IntersectionPoint1.m_fY = t0 * (y1_minus_y0)+y0;
		_rv3IntersectionPoint1.m_fZ = t0 * (z1_minus_z0)+z0;

		_rv3IntersectionPoint2.m_fX = t1 * (x1_minus_x0)+x0;
		_rv3IntersectionPoint2.m_fY = t1 * (y1_minus_y0)+y0;
		_rv3IntersectionPoint2.m_fZ = t1 * (z1_minus_z0)+z0;

		return INTERSECTION_TWO;

	}
	else if (discriminant == 0) {

		float t = (-1.0f * b) / (2.0f * a);

		_rv3IntersectionPoint1.m_fX = t * (x1_minus_x0)+x0;
		_rv3IntersectionPoint1.m_fY = t * (y1_minus_y0)+y0;
		_rv3IntersectionPoint1.m_fZ = t * (z1_minus_z0)+z0;

		return INTERSECTION_ONE;

	}
	else {

		return INTERSECTION_NONE;

	}

}

/***********************
* name of the function: IsLinePlaneIntersection
* @author: Juan Rodriguez
* @parameter: _krLine the 3D line
* @parameter: _krPlane the 3D plane
* @parameter: _rv3IntersectionPoint the point where the plane and the line intersect
* @return: boolean
********************/
bool IsLinePlaneIntersection(const T3DLine& _krLine,
	const TPlane& _krPlane,
	TVector3& _rv3IntersectionPoint) {

	// Normal = <a, b, c>
	float a = _krPlane.m_v3normal.m_fX;
	float b = _krPlane.m_v3normal.m_fY;
	float c = _krPlane.m_v3normal.m_fZ;

	// Plane = <x, y, z>
	float x = _krPlane.m_v3point.m_fX;
	float y = _krPlane.m_v3point.m_fY;
	float z = _krPlane.m_v3point.m_fZ;

	// Equation for getting d out of the normal:
	// ax + by + cz = -d 
	// We don't use -d because the change of the sign will be applied later
	float d = ((a * x) + (b * y) + (c * z));

	// Line
	float x0 = _krLine.m_v3q.m_fX;
	float y0 = _krLine.m_v3q.m_fY;
	float z0 = _krLine.m_v3q.m_fZ;
	float x1_minus_x0 = x0 + (_krLine.m_v3v.m_fX - x0);
	float y1_minus_y0 = y0 + (_krLine.m_v3v.m_fY - y0);
	float z1_minus_z0 = z0 + (_krLine.m_v3v.m_fZ - z0);

	//            ax    +    by    +    cz    + d
	float t = -((a * x0) + (b * y0) + (c * z0) + d) / (DotProduct(_krPlane.m_v3normal, _krLine.m_v3v));

	if (t != 0) {

		_rv3IntersectionPoint.m_fX = x0 + (t * x1_minus_x0);
		_rv3IntersectionPoint.m_fY = y0 + (t * y1_minus_y0);
		_rv3IntersectionPoint.m_fZ = z0 + (t * z1_minus_z0);

		return true;

	}
	else {

		return false;

	}

}

/***********************
* name of the function: IsIntersection
* @author: Juan Rodriguez
* @parameter: _krLine1 the first 3D line
* @parameter: _krLine2 the second 3D line
* @return: boolean
********************/
bool IsIntersection(const T3DLine& _krLine1,
	const T3DLine& _krLine2) {

	// Line 1
	float l1x0 = _krLine1.m_v3q.m_fX;
	float l1y0 = _krLine1.m_v3q.m_fY;
	float l1z0 = _krLine1.m_v3q.m_fZ;

	float l1x1 = _krLine1.m_v3v.m_fX;
	float l1y1 = _krLine1.m_v3v.m_fY;
	float l1z1 = _krLine1.m_v3v.m_fZ;

	// Line 2
	float l2x0 = _krLine2.m_v3q.m_fX;
	float l2y0 = _krLine2.m_v3q.m_fY;
	float l2z0 = _krLine2.m_v3q.m_fZ;

	float l2x1 = _krLine2.m_v3v.m_fX;
	float l2y1 = _krLine2.m_v3v.m_fY;
	float l2z1 = _krLine2.m_v3v.m_fZ;

	float d = l1x0 + l1x1 - (l2x0 + l2x1);

	if (d > 0) {

		return true;

	}
	else {

		return false;

	}

}

/***********************
* name of the function: ComputeIntersectionBetweenLines
* @author: Juan Rodriguez
* @parameter: _krLine1 the first 3D line
* @parameter: _krLine2 the second 3D line
* @parameter: _rIntersectionPoint the point of intersection
* @return: _rIntersectionPoint
********************/
TVector3& ComputeIntersectionBetweenLines(const T3DLine& _krLine1,
	const T3DLine& _krLine2,
	TVector3& _rIntersectionPoint) {

	// Line 1
	float l1x0 = _krLine1.m_v3q.m_fX;
	float l1y0 = _krLine1.m_v3q.m_fY;
	float l1z0 = _krLine1.m_v3q.m_fZ;

	float l1x1 = _krLine1.m_v3v.m_fX;
	float l1y1 = _krLine1.m_v3v.m_fY;
	float l1z1 = _krLine1.m_v3v.m_fZ;

	// Line 2
	float l2x0 = _krLine2.m_v3q.m_fX;
	float l2y0 = _krLine2.m_v3q.m_fY;
	float l2z0 = _krLine2.m_v3q.m_fZ;

	float l2x1 = _krLine2.m_v3v.m_fX;
	float l2y1 = _krLine2.m_v3v.m_fY;
	float l2z1 = _krLine2.m_v3v.m_fZ;

	float d = l1x0 + l1x1 - (l2x0 + l2x1);

	_rIntersectionPoint.m_fX = l1x0 + (d * l1x1);
	_rIntersectionPoint.m_fY = l1y0 + (d * l1y1);
	_rIntersectionPoint.m_fZ = l1z0 + (d * l1z1);

	return _rIntersectionPoint;

}

/***********************
* name of the function: IsInFieldOfView
* @author: Juan Rodriguez
* @parameter: _krCameraPosition the position of the camera
* @parameter: _krCameraDirection the direction where the camera looks at
* @parameter: _kfFieldOfViewInRadians the field of view given in radians
* @parameter: _krObjectPosition the position of the object
* @return: boolean
********************/
bool IsInFieldOfView(const TVector2& _krCameraPosition,
	const TVector2& _krCameraDirection,
	const float _kfFieldOfViewInRadians,
	const TVector2& _krObjectPosition) {

	// Normalised vectors
	TVector2 cameraDirectionNormalised;
	TVector2 vectorObjCameraPosNormalised;

	// Magnitude of our vectors (not yet normalised)
	float cameraDirectionMagnitude;
	float vectorObjCameraPosMagnitude;

	// Dot product of our two normalised vectors
	float normalisedVectorDotProduct;

	// Vector from the object to the position of the camera
	TVector2 vectorCameraPosObject;

	float angle;
	
	// Vector from the camera position to the object
	vectorCameraPosObject.m_fX = _krObjectPosition.m_fX - _krCameraPosition.m_fX;
	vectorCameraPosObject.m_fY = _krObjectPosition.m_fY - _krCameraPosition.m_fY;

	// Magnitude of our direction and object to camera position vectors
	cameraDirectionMagnitude = sqrt(pow(_krCameraDirection.m_fX, 2) + pow(_krCameraDirection.m_fY, 2));
	vectorObjCameraPosMagnitude = sqrt(pow(vectorCameraPosObject.m_fX, 2) + pow(vectorCameraPosObject.m_fY, 2));

	// Our normalised vectors
	cameraDirectionNormalised.m_fX = _krCameraDirection.m_fX / cameraDirectionMagnitude;
	cameraDirectionNormalised.m_fY = _krCameraDirection.m_fY / cameraDirectionMagnitude;

	vectorObjCameraPosNormalised.m_fX = vectorCameraPosObject.m_fX / vectorObjCameraPosMagnitude;
	vectorObjCameraPosNormalised.m_fY = vectorCameraPosObject.m_fY / vectorObjCameraPosMagnitude;
	
	// Recalling from the formula of dot product, we proceed to get the dot product of our normalised vectors
	normalisedVectorDotProduct = (cameraDirectionNormalised.m_fX * vectorObjCameraPosNormalised.m_fX) +
		(cameraDirectionNormalised.m_fY * vectorObjCameraPosNormalised.m_fY);

	// then apply with the form of:
	// 0 = acos(D' * V')
	angle = acos(normalisedVectorDotProduct);

	// If our angle is greater than half of our field of view, the object isn't visible
	if (angle > _kfFieldOfViewInRadians / 2) {

		return false;

	}
	else {

		return true;

	}

}

/***********************
* name of the function: IsSurfaceLit
* @author: Juan Rodriguez
* @parameter: _krCameraPosition the position of the camera
* @parameter: _krCameraDirection the direction where the camera looks at
* @parameter: _kfFieldOfViewInRadians the field of view given in radians
* @parameter: _krObjectPosition the position of the object
* @return: boolean
********************/
bool IsSurfaceLit(const TVector3& _krPointOnSurface,
	const TVector3& _krLightSourcePosition,
	const TTriangle3& _krSurface) {
	
	/*
	Got the surface normal of the triangle. 
	Then got the vector from the light source to the point on the surface.
	Then compared the direction of the two to check if the light source is in front of or behind the triangle. 
	Not sure if it's right, it's the best I could think of.
	*/

	TVector3 N; // Normal
	N = FindTriangleNormal(_krSurface, N);

	TVector3 vectorLightPointOnSurface;

	// Making a vector from the light source to the point on the surface
	vectorLightPointOnSurface = Subtract(_krPointOnSurface, _krLightSourcePosition, vectorLightPointOnSurface);


	
	return false;

}

/***********************
* name of the function: FindTriangleNormal
* @author: Juan Rodriguez
* @parameter: _krCameraPosition the position of the camera
* @parameter: _krCameraDirection the direction where the camera looks at
* @parameter: _kfFieldOfViewInRadians the field of view given in radians
* @parameter: _krObjectPosition the position of the object
* @return: boolean
********************/
TVector3& FindTriangleNormal(const TTriangle3& _krTriangle,
	TVector3& _rNormal) {

	TVector3 V; // P2 - P1
	TVector3 W; // P3 - P1

	// First point 
	TVector3 P1;
	P1.m_fX = _krTriangle.m_v3p1.m_fX;
	P1.m_fY = _krTriangle.m_v3p1.m_fY;
	P1.m_fZ = _krTriangle.m_v3p1.m_fZ;

	// Second point
	TVector3 P2;
	P2.m_fX = _krTriangle.m_v3p2.m_fX;
	P2.m_fY = _krTriangle.m_v3p2.m_fY;
	P2.m_fZ = _krTriangle.m_v3p2.m_fZ;

	// Third point
	TVector3 P3;
	P3.m_fX = _krTriangle.m_v3p3.m_fX;
	P3.m_fY = _krTriangle.m_v3p3.m_fY;
	P3.m_fZ = _krTriangle.m_v3p3.m_fZ;

	// V = P2 - P1
	V.m_fX = Subtract(P2, P1, V).m_fX;
	V.m_fY = Subtract(P2, P1, V).m_fY;
	V.m_fZ = Subtract(P2, P1, V).m_fZ;
		
	// W = P3 - P1
	W.m_fX = Subtract(P3, P1, W).m_fX;
	W.m_fY = Subtract(P3, P1, W).m_fY;
	W.m_fZ = Subtract(P3, P1, W).m_fZ;

	// _rNormal = V x W
	_rNormal = CrossProduct(V, W, _rNormal);

	return _rNormal;

}

/***********************
* name of the function: RotateTriangleAroundPoint
* @author: Juan Rodriguez
* @parameter: _krCameraPosition the position of the camera
* @parameter: _krCameraDirection the direction where the camera looks at
* @parameter: _kfFieldOfViewInRadians the field of view given in radians
* @parameter: _krObjectPosition the position of the object
* @return: boolean
********************/
TTriangle2& RotateTriangleAroundPoint(const TTriangle2& _krTriangle,
	const float _kfRotAngleInRadians,
	const TVector2& _krRotAroundPoint,
	TTriangle2& _rRotatedTriangle) {

	// First point on the triangle
	float x1 = _krTriangle.m_v2p1.m_fX;
	float y1 = _krTriangle.m_v2p1.m_fY;

	// Second point on the triangle
	float x2 = _krTriangle.m_v2p2.m_fX;
	float y2 = _krTriangle.m_v2p2.m_fY;

	// Third point on the triangle
	float x3 = _krTriangle.m_v2p3.m_fX;
	float y3 = _krTriangle.m_v2p3.m_fY;
	
	// Rotation points
	float rx = _krRotAroundPoint.m_fX;
	float ry = _krRotAroundPoint.m_fY;

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
	rx = x1 * cosf(_kfRotAngleInRadians / 180.0f * float(M_PI)) - y1 * sinf(_kfRotAngleInRadians / 180.0f * float(M_PI));
	ry = y1 * cosf(_kfRotAngleInRadians / 180.0f * float(M_PI)) + x1 * sinf(_kfRotAngleInRadians / 180.0f * float(M_PI));
	x1 += rx;
	y1 += ry;

	// Rotating the second point
	rx = x2 * cosf(_kfRotAngleInRadians / 180.0f * float(M_PI)) - y2 * sinf(_kfRotAngleInRadians / 180.0f * float(M_PI));
	ry = y2 * cosf(_kfRotAngleInRadians / 180.0f * float(M_PI)) + x2 * sinf(_kfRotAngleInRadians / 180.0f * float(M_PI));
	x2 += rx;
	y2 += ry;

	// Rotating the third point
	rx = x3 * cosf(_kfRotAngleInRadians / 180.0f * float(M_PI)) - y3 * sinf(_kfRotAngleInRadians / 180.0f * float(M_PI));
	ry = y3 * cosf(_kfRotAngleInRadians / 180.0f * float(M_PI)) + x3 * sinf(_kfRotAngleInRadians / 180.0f * float(M_PI));
	x3 += rx;
	y3 += ry;

	// Asigning our rotated points to our final triangle
	// First point
	_rRotatedTriangle.m_v2p1.m_fX = x1;
	_rRotatedTriangle.m_v2p1.m_fY = y1;

	// Second point
	_rRotatedTriangle.m_v2p2.m_fX = x2;
	_rRotatedTriangle.m_v2p2.m_fY = y2;

	// Third point
	_rRotatedTriangle.m_v2p3.m_fX = x3;
	_rRotatedTriangle.m_v2p3.m_fY = y3;

	return _rRotatedTriangle;
	
}