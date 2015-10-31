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

	if ((_krA.m_fX == _krB.m_fX) && (_krA.m_fY == _krB.m_fY) && (_krA.m_fZ == _krB.m_fZ) ) {
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

	return (sqrt( 
			pow(_krA.m_fX, 2) + pow(_krA.m_fY, 2) + pow(_krA.m_fZ, 2)
		));

}

/***********************
* name of the function: DotProdcut
* @author: Juan Rodriguez
* @parameter: _krA the first vector to multiply
* @parameter: _krB the second vector to multiply
* @return: The product of the multiplication of both vectors' XYZ components
********************/
float DotProduct(const TVector3& _krA, const TVector3& _krB) {

	return (
			(_krA.m_fX * _krB.m_fX) + (_krA.m_fY * _krB.m_fY) + (_krA.m_fZ * _krB.m_fZ)
		);

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

	_rResultant.m_fX =       (_krA.m_fY * _krB.m_fZ) - (_krB.m_fY * _krA.m_fZ);
	_rResultant.m_fY = -1.0f * ((_krA.m_fX * _krB.m_fZ) - (_krB.m_fX * _krA.m_fZ));
	_rResultant.m_fZ =       (_krA.m_fX * _krB.m_fY) - (_krB.m_fX * _krA.m_fY);

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
			acos( dotProduct / (magnitude_krA * magnitude_krB) )
		);

}

/***********************
* name of the function: ComputeAngleBetween
* @author: Juan Rodriguez
* @parameter: _krA the first vector
* @parameter: _krB the second vector
* @return: angle between the two 3D Vectors
********************/
float ComputeAngleBetween(const TVector3& _krA,
	const TVector3& _krB) {

	return (
			acos( DotProduct(_krA, _krB) / (Magnitude(_krA) * Magnitude(_krB)) )
		);

}

/***********************
* name of the function: ComputeDistancePointToLine
* @author: Juan Rodriguez
* @parameter: 
* @parameter: 
* @return: 
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
* @parameter: 
* @parameter: 
* @return: 
********************/
float ComputeDistancePointToPlane(const TPlane& _krPlane,
	const TVector3& _krPoint) {

	TVector3 vector;
	TVector3 proj;

	vector = Subtract(_krPoint, _krPlane.m_v3point, vector);
	proj = Projection(vector, _krPlane.m_v3normal, proj);

	return (
			Magnitude(proj)
		);

}

/***********************
* name of the function: ComputeDistancePointToSphere
* @author: Juan Rodriguez
* @parameter: 
* @parameter: 
* @return: 
********************/
//Distance between point and center of the spheres
float ComputeDistancePointToSphere(const TSphere& _krSphere,
	const TVector3& _krPoint) {

	float SQRT = sqrt(pow(_krSphere.m_v3center.m_fX - _krPoint.m_fX, 2) +
		pow(_krSphere.m_v3center.m_fY - _krPoint.m_fY, 2) +
		pow(_krSphere.m_v3center.m_fZ - _krPoint.m_fZ, 2) );

	return (
			abs( SQRT - _krSphere.m_fRadius )
		);

}

/***********************
* name of the function: ComputeDistanceCircleToCircle
* @author: Juan Rodriguez
* @parameter: 
* @parameter: 
* @return: 
********************/
//Distance between center of the circles
float ComputeDistanceCircleToCircle(const TCircle& _krCircle1,
	const TCircle& _krCircle2) {

	return (
			sqrt(pow(_krCircle1.m_v2center.m_fX - _krCircle2.m_v2center.m_fX, 2) + 
				 pow(_krCircle1.m_v2center.m_fY - _krCircle2.m_v2center.m_fY, 2))
		);

}

/***********************
* name of the function: ComputeDistanceCircleToTriangle
* @author: Juan Rodriguez
* @parameter: 
* @parameter: 
* @return: 
********************/
//Distance between center of the circle and triangle
float ComputeDistanceCircleToTriangle(const TCircle& _krCircle,
	const TTriangle2& _krTriangle) {

	float triangleCentroidX = (_krTriangle.m_v2p1.m_fX + _krTriangle.m_v2p2.m_fX + _krTriangle.m_v2p3.m_fX) / 3;
	float triangleCentroidY = (_krTriangle.m_v2p1.m_fY + _krTriangle.m_v2p2.m_fY + _krTriangle.m_v2p3.m_fY) / 3;
	float SQRT = sqrt(pow(_krCircle.m_v2center.m_fX - triangleCentroidX, 2) + pow(_krCircle.m_v2center.m_fY - triangleCentroidY, 2));

	return (
			abs(SQRT - _krCircle.m_fRadius)
		);

}

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
	float x1_minus_x0 = _krLine.m_v3v.m_fX - _krLine.m_v3q.m_fX;
	float y1_minus_y0 = _krLine.m_v3v.m_fY - _krLine.m_v3q.m_fY;
	float z1_minus_z0 = _krLine.m_v3v.m_fZ - _krLine.m_v3q.m_fZ;

	float a = pow(x1_minus_x0, 2) + pow(y1_minus_y0, 2) + pow(z1_minus_z0, 2);

	float b = ((-2 * h)*(x1_minus_x0) + (2 * x0)*(x1_minus_x0)) -
			  ((-2 * k)*(y1_minus_y0) + (2 * y0)*(y1_minus_y0)) -
			  ((-2 * j)*(z1_minus_z0) + (2 * z0)*(z1_minus_z0));

	float c = pow(h, 2) + pow(k, 2) + pow(j, 2) -
			  ((-2 * h * x0) + pow(x0, 2)) - 
			  ((-2 * k * y0) + pow(y0, 2)) - 
			  ((-2 * j * z0) + pow(z0, 2)) - 
			  pow(r, 2);

	float discriminant = sqrt(pow(b, 2) - (4.0f * a * c));

	if (discriminant > 0) {

		float t0 = (((-1.0f * b) + discriminant) / (2.0f * a));
		float t1 = (((-1.0f * b) - discriminant) / (2.0f * a));

		_rv3IntersectionPoint1.m_fX = t0 * (x1_minus_x0) + x0;
		_rv3IntersectionPoint1.m_fY = t0 * (y1_minus_y0) + y0;
		_rv3IntersectionPoint1.m_fZ = t0 * (z1_minus_z0) + z0;

		_rv3IntersectionPoint2.m_fX = t1 * (x1_minus_x0) + x0;
		_rv3IntersectionPoint2.m_fY = t1 * (y1_minus_y0) + y0;
		_rv3IntersectionPoint2.m_fZ = t1 * (z1_minus_z0) + z0;

		return INTERSECTION_TWO;

	}
	else if (discriminant == 0) {

		float t = (-1.0f * b) / (2.0f * a);

		_rv3IntersectionPoint1.m_fX = t * (x1_minus_x0) + x0;
		_rv3IntersectionPoint1.m_fY = t * (y1_minus_y0) + y0;
		_rv3IntersectionPoint1.m_fZ = t * (z1_minus_z0) + z0;

		return INTERSECTION_ONE;

	}
	else {

		return INTERSECTION_NONE;

	}

	

}