//////////////////////////////////////////////////////////////////////////////////
//																																							//
// Copyright (C) 2010 Engin Turetken																						//
//																																							//
// This program is free software: you can redistribute it and/or modify         //
// it under the terms of the version 3 of the GNU General Public License        //
// as published by the Free Software Foundation.                                //
//                                                                              //
// This program is distributed in the hope that it will be useful, but          //
// WITHOUT ANY WARRANTY; without even the implied warranty of                   //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU             //
// General Public License for more details.                                     //
//                                                                              //
// You should have received a copy of the GNU General Public License            //
// along with this program. If not, see <http://www.gnu.org/licenses/>.         //
//                                                                              //
// Contact <engin.turetken@epfl.ch> for comments & bug reports                  //
//////////////////////////////////////////////////////////////////////////////////

#ifndef __itkPolyLineParametricPathExtended_h
#define __itkPolyLineParametricPathExtended_h

#include "itkParametricPath.h"
#include "itkVectorContainer.h"
#include "itkContinuousIndex.h"
#include "itkIndex.h"
#include "itkOffset.h"
#include "itkVector.h"
#include "itkArray.h"
#include "itkMath.h"
#include "itkImage.h"
#include "itkPathConstIterator.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkCircularShapedNeighborhoodIterCreator.h"
#include "itkLinearInterpolateImageFunction.h"

#include <algorithm>

#include <boost/config.hpp>
#include <boost/serialization/vector.hpp>

// TODO
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace itk
{
	/** \class PolyLineParametricPathExtended
	 * \brief  Represents a path of line segments through ND Space
	 *
	 */
	template <unsigned int VDimension>
	class PolyLineParametricPathExtended : public ParametricPath< VDimension >
	{
	public:
		/** Standard class typedefs. */
		typedef PolyLineParametricPathExtended      Self;
		typedef ParametricPath<VDimension>  Superclass;
		typedef SmartPointer<Self>          Pointer;
		typedef SmartPointer<const Self>    ConstPointer;
		
		typedef PolyLineParametricPathExtended<(VDimension > 0) ? VDimension-1 : VDimension> NMinus1DPathType;
		typedef PolyLineParametricPathExtended<(VDimension < 5) ? VDimension+1 : VDimension> NPlus1DPathType;
		
		/** Dimension constant. */
		itkStaticConstMacro(Dimension, unsigned int, VDimension);
		
		/** New() method for dynamic construction */
		itkNewMacro( Self );
		
		/** Run-time type information (and related methods). */
		itkTypeMacro(PolyLineParametricPathExtended, ParametricPath);
		
		/** Input type */
		typedef typename Superclass::InputType  InputType;
		
		/** Output type */
		typedef typename Superclass::OutputType OutputType;
		
		
		/** Basic data-structure types used */
		typedef ContinuousIndex<double,VDimension>    ContinuousIndexType;
		typedef Index<  VDimension >                  IndexType;
		typedef ImageRegion<  VDimension >            ImageRegionType;
		typedef Offset< VDimension >                  OffsetType;
		typedef Size< VDimension >										SizeType;
		typedef typename 
		ImageBase<VDimension>::SpacingType						SpacingType;
		typedef Point<double,VDimension>              PointType;
		typedef Vector<double,VDimension>             VectorType;
		typedef ContinuousIndexType                   VertexType;
		typedef VectorContainer<unsigned, VertexType> VertexListType;		// TODO use an std::vector instead!
		typedef typename VertexListType::Pointer      VertexListPointer;
		
		typedef double																RadiusType;
		typedef std::vector<RadiusType>								RadiusListType;
			
		virtual OutputType Evaluate( const InputType & input ) const;
		virtual RadiusType EvaluateRadius( const InputType & input ) const;
	
		inline virtual bool AddVertex( const VertexType & vertex, RadiusType radi = 0 )
		{
			// Do not add the vertex if it is too close to the end vertex.
			if( m_VertexList->Size() > 0 )
			{
				const VertexType& endVertex = m_VertexList->ElementAt(m_VertexList->Size() - 1);
				
				if( vertex.SquaredEuclideanDistanceTo(endVertex) <  m_Epsilon )
				{
					return false;
				}
			}
			
			m_VertexList->InsertElement( m_VertexList->Size(), vertex );
			m_RadiusList.push_back( radi );
			this->Modified();
			
			return true;
		}
		
		/** Creates a clone of the path and returns it. */
		virtual Pointer Clone() const;
		
		/** Copies the given path input to this path. */
		virtual void Copy(const Self* path);
		
		/** Appends a given path to the end of this one. */
		virtual void AppendEnd(const Self* path);
		
		/** Appends a given path to the end of this one. */
		virtual void Split(InputType splitPoint, Self* path1, Self* path2)  const;
		
		/** Reverses the direction of the path. */
		virtual void Reverse();
		
		/** Extracts a segment of this path. */
		virtual void SubPath(InputType startPoint, InputType endPoint, Self* path)  const;
		
		/** Rounds the continuous indices this path traverses. */
		virtual void RoundIndices(  bool removeRepeatedVertices = true );
		
		/** Resamples the path at the given step intervals in world coordinate system. */
		template<class TImage>
		void Resample(double stepInWorldCoords, const TImage* image);
		
		// TODO Add downsample and upsample functions as well.

		/** Smooths the vertex locations and the radius values. */
		template<class TImage>
		void SmoothVertexLocationsAndRadii(double avgWindowRadiusInWorldCoords,
																			 const TImage* image);
		
		/** Smooths the vertex locations. */
		template<class TImage>
		void SmoothVertexLocations(double avgWindowRadiusInWorldCoords,
															 const TImage* image);
		
		/** Smooths the radius values along the path. */
		template<class TImage>
		void SmoothRadiusValues(double avgWindowRadiusInWorldCoords,
														const TImage* image);
		
		/** 
		 * Computes bounding image region of the tubular path with 
		 * the radius information taken into account. 
		 */
		virtual ImageRegionType 
		ComputeBoundingImageRegionWithRadius(const SpacingType& spacing) const;
		/** 
		 * Computes bounding image region of the centerline of the 
		 * path without the radius information taken into account. 
		 */
		virtual ImageRegionType ComputeBoundingImageRegionWithoutRadius() const;

		/** Compute mean vertex location of the path. */
		virtual VertexType ComputeMeanVertex() const;
		
		/** 
		 * Compute centroid location of the path. The centroid location 
		 * is taken as the weighted mean of the linear segments. The weight
		 * of a segment is taken as the length of it in world coordinates. 
		 */
		template<class TImage>
		VertexType ComputeCentroidVertex(TImage* image) const;
		
		/** Compute minimum/maximum/mean/centroid radius along the path. */
		virtual RadiusType ComputeMinRadius() const;
		virtual RadiusType ComputeMaxRadius() const;
		virtual RadiusType ComputeMeanRadius() const;
		template<class TImage>
		RadiusType ComputeCentroidRadius(TImage* image) const;
		
		/** 
		 * Scale and shift the radius values of the path. 
		 * The new values becomes r_new = scale * r_old + shift.
		 */
		virtual void ScaleAndShiftRadii(double scaleFactor = 1,
																		double shift = 0);
		
		/** 
		 * Scale and shift the vertex coordinates of the path that are 
		 * defined in the continuous image domain. The new coordinate values 
		 * becomes x_new = scale * x_old + shift.
		 */
		virtual void ScaleAndShiftVertices(double scaleFactor = 1,
																			 double shift = 0);
		
		/** 
		 * Get the list of unique image indices that this path traverses.
		 * If the path visits an index multiple times, only the first one 
		 * is kept in the array and all the subsequent ones are removed.
		 * HACK: This should ideally be an option in the path iterator (e.g., 
		 * a flag that determines whether repeated indices will be visited 
		 * or not.) Alternatively, it can be a path filter.
		 */
		template<class TImage>
		void GetUniqueImageIndices(const TImage* image, 
															 std::vector<IndexType>& uniqueIndices) const;

		/** 
		 * Create a path from the list of unique image indices that this path 
		 * traverses. The image indices of the source and the target point is kept 
		 * the same. If there exist a loop that includes the end index, then the 
		 * image indices on the loop are removed.
		 */
		template<class TImage>
		Pointer CreatePathFromUniqueImageIndices(const TImage* image) const;

		/** 
		 * Compute the minimum Euclidean distance to path points in image coordinates.
		 * Optionally, it also returns the path input for which the distance is minimum.
		 */
		virtual std::pair<double,  InputType>
		GetEuclDistToPointInImageCoords(const VertexType& vertex)  const;
		
		/** 
		 * Compute the minimum Euclidean distance to path points in world coordinates.
		 * Optionally, it also returns the path input for which the distance is minimum.
		 */
		template<class TImage>
		std::pair<double,  InputType> 
		GetEuclDistToPointInWorldCoords(const VertexType& vertex, 
																		const TImage* image)  const;

		/** 
		 * Checks whether this path contains the given point. 
		 * Returns true if it contains the point, or false 
		 * otherwise.
		 */
		template<class TImage>
		bool DoesContainPoint(const VertexType& vertex,
													const TImage* image,
													bool excludeTails = true) const;

		/** 
		 * Extract the set of path points whose distances to the given point is 
		 * minimum in world coordinate system. The returned path points are unique. 
		 * An additional constraint imposed on the path points is that they must 
		 * include the given point, that is the radius value at the path point must 
		 * be greater than or equal to the distance to the point.
		 */
		template <class TImage>
		void
		ExtractClosestPathPointsThatContainsPointInWorldCoords(const VertexType& vertex, 
																													 const TImage* image,
																													 std::pair< double, std::vector<InputType> >& outputPairs,
																													 bool excludeTails = true) const;
		
		
		/** 
		 * Checks whether this path contains the centerline 
		 * vertices of the given path. Returns true if it contains 
		 * all the vertices, or false otherwise.
		 */
		template<class TImage>
		bool DoesContainPointsOfPath(const Self* path,
																 const TImage* image,
																 bool excludeTails = true) const;

		/** 
		 * Approximately computes the length of this path's connected 
		 * longest segment that resides inside the given path in world 
		 * coordinates and returns it. 
		 */
		template <class TImage>
		double ComputeMaxPathSegmLengthOutsidePath(const Self* path,
																							 const TImage* image,
																							 bool excludeTails = true) const;		
		
		/** 
		 * Approximately computes the total length of this path's 
		 * centerline inside the given path in world coordinates 
		 * and returns it. 
		 */
		template <class TImage>
		double ComputePathLengthInsidePath(const Self* path,
																			 const TImage* image,
																			 bool excludeTails = true) const;
		
		/** 
		 * Approximately computes the total length of this path's 
		 * centerline outside the given path in world coordinates 
		 * and returns it. 
		 */
		template <class TImage>
		double ComputePathLengthOutsidePath(const Self* path,
																			  const TImage* image,
																				bool excludeTails = true) const;
		
		/** 
		 * Approximately computes the total length of the path's 
		 * centerline inside the given mask image in world coordinates. 
		 * Mask image membership is check by casting the pixel value to a bool.
		 */
		template <class TImage>
		double ComputePathLengthInsideMask(const TImage* maskImage) const;
		
		/** 
		 * Approximately computes the total length of the path's 
		 * centerline outside the given mask image in world coordinates. 
		 * Mask image membership is check by casting the pixel value to a bool.
		 */
		template <class TImage>
		double ComputePathLengthOutsideMask(const TImage* maskImage) const;
		
		/** 
		 * Convert an N-D path to an (N-1)-D path by discarding the radius 
		 * value and treating the last dimension as the world coordinate value 
		 * along the radius dimension.
		 */
		virtual typename NMinus1DPathType::Pointer 
		ConvertToNMinus1DPath(double radiusOrigin, double radiusSpacing) const;
		
		/** 
		 * Convert an N-D path to an (N+1)-D path by treating the radius  
		 * as the world coordinate value of the last dimension.
		 */
		virtual typename NPlus1DPathType::Pointer 
		ConvertToNPlus1DPath(double radiusOrigin, double radiusSpacing) const;
		
		/** 
		 * Returns the number intersections between this path and a given one. To 
		 * check for intersections, the points on the paths are first rounded and then 
		 * compared. The number of intersections can be at most the number of point on
		 * this path. 
		 * 
		 * TODO This implementation computes the intersection of the centerline
		 * points only and discards the radius values.
		 */
		virtual unsigned long GetNoOfIntersectionPoints(const Self* path) const;
		
		/** 
		 * Computes the numberß of pixels in both the intersection and union 
		 * regions of this path and the given one.
		 */
		template<class TImage>
		void ComputeIntersectionAndUnionInImageCoords(const Self* path,
																									const TImage* image,
																									unsigned long& intersectionVolume,
																									unsigned long& unionVolume,
																									bool excludeTails = true) const;
		/** 
		 * Computes the volumes of the both intersection and union regions  of 
		 * this path and the given one in world coordinates.
		 */		
		template<class TImage>
		void ComputeIntersectionAndUnionInWorldCoords(const Self* path,
																									const TImage* image,
																									double& intersectionVolume,
																									double& unionVolume,
																									bool excludeTails = true) const;
		
		
		/** 
		 * Tests if this path is completely inside the given region.
		 */
		bool IsInside(const ImageRegionType& region) const
		{
			return region.IsInside(ComputeBoundingImageRegionWithoutRadius());
		}
			
		/** 
		 * Computes the sum of the pixels along the path. It is assumed that 
		 * the + operator for the type TSum is defined. When the 
		 * interpolatePixels flag is true, bilinear interpolation of the pixels
		 * values are used, otherwise nearest neighbor interpolation is used.
		 */
		template<class TSum, class TImage> 
		TSum PixelSum(const TImage* image,
									bool interpolatePixels = true)
		{
			typedef PathConstIterator<TImage, Self>					IterType;
			typedef	typename TImage::PixelType							PixelType;
			typedef TSum																		SumType;
			typedef LinearInterpolateImageFunction<TImage, 
			typename ContinuousIndexType::CoordRepType>			InterpolatorType;
			
			SumType pixelSum = NumericTraits<SumType>::ZeroValue();
			IterType iter(image, this);
			if( interpolatePixels )
			{
				typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
				
				interpolator->SetInputImage( image );
				
				for(iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
				{
					PixelType val = 
					interpolator->EvaluateAtContinuousIndex
					(this->Evaluate(iter.GetPathPosition()));
					pixelSum = pixelSum + static_cast<SumType>(val);
				}
			}
			else
			{
				for(iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
				{
					pixelSum = pixelSum + static_cast<SumType>(iter.Get());
				}
			}
			
			return pixelSum;
		}
		
		/** 
		 * Computes the sum of the pixels along the path similar to the PixelSum() 
		 * function, but weights pixel values by the Eucliden distance of their 
		 * connecting path segments in the world coordinate system.
		 */
		template<class TSum, class TImage> 
		TSum PixelSumEuclDistWeighted(const TImage* image,
																	bool interpolatePixels = true)
		{
			typedef PathConstIterator<TImage, Self>					IterType;
			typedef	typename TImage::PixelType							PixelType;
			typedef TSum																		SumType;
			typedef LinearInterpolateImageFunction<TImage, 
			typename ContinuousIndexType::CoordRepType>			InterpolatorType;
			
			SumType pixelSum = NumericTraits<SumType>::ZeroValue();
			IterType iter(image, this);
			typename InterpolatorType::Pointer interpolator;
			if( interpolatePixels )
			{
				interpolator = InterpolatorType::New();
				interpolator->SetInputImage( image );
			}
			
			double distInWorldCoord;
			InputType prevPos;
			InputType currentPos;
			PixelType prevVal;
			PixelType currentVal;
			bool isAtBegin = true;
			for(iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
			{
				currentPos = iter.GetPathPosition();
				if( interpolatePixels )
				{
					currentVal = 
					interpolator->EvaluateAtContinuousIndex(this->Evaluate(currentPos));
				}
				else
				{
					currentVal = iter.Get();
				}
								
				if( !isAtBegin )
				{
					// Compute the distance between the previous position and 
					// the current one in world coordinates.
					distInWorldCoord = this->GetDistanceInWorldCoords(prevPos,
																														currentPos,
																														image);
					
					// Accumulate the average of the pixel value pairs of successive points 
					// weighted by the euclidean length of the path segment between them.
					pixelSum = pixelSum + 
					static_cast<SumType>(0.5 * distInWorldCoord * (prevVal + currentVal));
				}
				
				prevPos = currentPos;
				prevVal = currentVal;
				isAtBegin = false;
			}
			
			return pixelSum;
		}
		
		/** Computes the length of the path in world coordinates. */
		template<class TImage>
		double GetLengthInWorldCoords(const TImage* image ) const
		{
			return GetDistanceInWorldCoords(this->StartOfInput(),
																			this->EndOfInput(),
																			image);
		}
		
		/** 
		 * Computes the Euclidean length between the given two points 
		 * on the path in the world coordinates.
		 */
		template<class TImage>
		double GetDistanceInWorldCoords(InputType startPoint, 
																		InputType endPoint, 
																		const TImage* image ) const
		{
			if( m_VertexList->Size() < 2 )
			{
				return 0.0;
			}
			
			if( endPoint > this->EndOfInput() )
			{
				itkExceptionMacro(<<"endPoint is larger than the end of the path.");
			}
			else if( endPoint < this->StartOfInput() )
			{
				itkExceptionMacro(<<"endPoint is smaller than the start of the path.");
			}
			if( startPoint > this->EndOfInput() )
			{
				itkExceptionMacro(<<"startPoint is larger than the end of the path.");
			}
			else if( startPoint < this->StartOfInput() )
			{
				itkExceptionMacro(<<"startPoint is smaller than the start of the path.");
			}
			
			if( endPoint < startPoint )
			{
				std::swap(startPoint, endPoint);
			}
			
			if( vcl_fabs(endPoint - startPoint) < m_Epsilon )
			{
				return 0.0;
			}
			
			double dist = 0;
			ContinuousIndexType currentContIndx = Evaluate(startPoint);
			PointType currentPointLoc;
			image->TransformContinuousIndexToPhysicalPoint(currentContIndx, currentPointLoc);
			unsigned long nextVertexPoint = ((unsigned long)(startPoint - this->StartOfInput())) + 1;
			ContinuousIndexType nextContIndx = m_VertexList->ElementAt( nextVertexPoint );
			PointType nextPointLoc;
			while((endPoint - ((InputType)nextVertexPoint) - this->StartOfInput()) > 
						m_Epsilon )
			{
				image->TransformContinuousIndexToPhysicalPoint(nextContIndx, nextPointLoc);
				dist += nextPointLoc.EuclideanDistanceTo(currentPointLoc);
				
				currentPointLoc = nextPointLoc;
				nextVertexPoint++;
				nextContIndx = m_VertexList->ElementAt( nextVertexPoint );
			}
			
			// Finally, add the distance from the current point to the end point.
			image->TransformContinuousIndexToPhysicalPoint(Evaluate(endPoint), nextPointLoc);
			dist += nextPointLoc.EuclideanDistanceTo(currentPointLoc);
			
			return dist;
		}
		
		/** Computes the length of the path in image coordinates. */
		double GetLengthInImageCoords() const
		{
			return GetDistanceInImageCoords(this->StartOfInput(),
																			this->EndOfInput());
		}
		
		/** 
		 * Computes the Euclidean length between the given two points 
		 * on the path in the image coordinates.
		 */
		double GetDistanceInImageCoords(InputType startPoint, 
																		InputType endPoint) const
		{
			if( m_VertexList->Size() < 2 )
			{
				return 0;
			}
			if( endPoint < startPoint )
			{
				std::swap(startPoint, endPoint);
			}
			
			if( endPoint > this->EndOfInput() )
			{
				endPoint = this->EndOfInput();
			}
			else if( endPoint < this->StartOfInput() )
			{
				endPoint = this->StartOfInput();
			}
			if( startPoint > this->EndOfInput() )
			{
				startPoint = this->EndOfInput();
			}
			else if( startPoint < this->StartOfInput() )
			{
				startPoint = this->StartOfInput();
			}
			
			if( vcl_fabs(endPoint - startPoint) < m_Epsilon )
			{
				return 0;
			}
			
			double dist = 0;
			ContinuousIndexType currentContIndx = Evaluate(startPoint);
			unsigned long nextVertexPoint = ((unsigned long)(startPoint - this->StartOfInput())) + 1;
			ContinuousIndexType nextContIndx = m_VertexList->ElementAt( nextVertexPoint );
			while((endPoint - ((InputType)nextVertexPoint) - this->StartOfInput()) > 
						m_Epsilon )
			{
				dist += nextContIndx.EuclideanDistanceTo(currentContIndx);
				
				currentContIndx = nextContIndx;
				nextVertexPoint++;
				nextContIndx = m_VertexList->ElementAt( nextVertexPoint );
			}
			
			// Finally, add the distance from the current point to the end point.
			dist += Evaluate(endPoint).EuclideanDistanceTo(currentContIndx);
			
			return dist;
		}
		
		
		// TODO Why the default step size is diameter and, for instance, not radius.
		// If stepSizeInWorldCoords <= 0, then path diameter is used as the step size.
		// This is the default behaviour, i.e., when no argument is not provided to 
		// the function.
		// If the path radius is zero, then 1.0 is used as the default step size.
		template<class TImage>
		VectorType EvaluateDirectionInWorldCoords(const InputType & input, 
																							const TImage* image,
																							const double & stepSizeInWorldCoords = 0.0) const
		{
			VectorType diffVector;
			if( m_VertexList->Size() < 2 )
			{
				diffVector.Fill(0.0);
				return diffVector;
			}
			
			InputType pathPoint = input;
			if( pathPoint < this->StartOfInput() )
			{
				pathPoint = this->StartOfInput();
			}
			if( pathPoint > this->EndOfInput() )
			{
				pathPoint = this->EndOfInput();
			}
			
			double stepSize = stepSizeInWorldCoords;
			if( stepSize < m_Epsilon )
			{
				stepSize = 2.0 * this->EvaluateRadius(pathPoint);
			}
			if( stepSize < m_Epsilon )
			{
				stepSize = 1.0;
			}

			// Find the previous and the next path points from which the direction 
			// vector will be extracted.
			VertexType currentIndex;
			VertexType prevIndex;
			VertexType nextIndex;
			InputType nextPathPoint;
			InputType prevPathPoint;
			double traversedDistance;
			do
			{
				this->TraverseDistanceInWorldCoords(pathPoint,
																						image,
																						stepSize,
																						false,
																						prevPathPoint,
																						traversedDistance);
				this->TraverseDistanceInWorldCoords(pathPoint,
																						image,
																						stepSize,
																						true,
																						nextPathPoint,
																						traversedDistance);
				
				currentIndex = this->Evaluate( pathPoint );
				prevIndex = this->Evaluate( prevPathPoint );
				nextIndex = this->Evaluate( nextPathPoint );
				
				stepSize = stepSize / 2.0;
			}
			while((prevIndex.SquaredEuclideanDistanceTo(nextIndex) < m_Epsilon) &&
						(prevIndex.SquaredEuclideanDistanceTo(currentIndex) < m_Epsilon) &&
						(stepSize > m_Epsilon));
			
			// Check if two points are the same.
			if(prevIndex.SquaredEuclideanDistanceTo(nextIndex) < m_Epsilon)
			{
				// If all points are at the same location, then 
				// throw an exception.
				if(prevIndex.SquaredEuclideanDistanceTo(currentIndex) < m_Epsilon)
				{
					itkExceptionMacro(<<"Path direction vector could not be evaluated!");
				}
				else // else take the previous point as the current one 
					// to compute the vector.
				{
					prevIndex = currentIndex;
				}
			}
			
			// Compute the unit-length path direction vector in world coordinates.
			PointType prevPoint;
			PointType nextPoint;
			image->TransformContinuousIndexToPhysicalPoint(prevIndex, prevPoint);
			image->TransformContinuousIndexToPhysicalPoint(nextIndex, nextPoint);
			
			diffVector = nextPoint - prevPoint;
			diffVector /= diffVector.GetNorm();
			
			return diffVector;
		}
		
		// Traverse each pixel along the path and sets its value
		// to the given one.
		template<class TImage>
		void OverlayCenterlines(TImage* image, 
														typename TImage::PixelType overlayVal) const
		{
			typedef TImage	ImageType;
			
			// Check if the path is completely inside the image buffered region.
			ImageRegionType region = image->GetBufferedRegion();
			if( !this->IsInside(region) )
			{
				itkExceptionMacro(<<"The path is not completely inside "
													<<"the image buffered region "
													<<region);
			}
			
			PathConstIterator<ImageType, Self> iter(image, this );
			for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
			{				
				image->SetPixel(iter.GetIndex(), overlayVal);
			}
		}
		
		
		// TODO Fix the bug in following function. The reasoning about the tail 
		// map image is wrong because a path can turn and end on itself. 
		// Use the DoesContainPoint method instead.
//		template<class TImage>
//		void Overlay(TImage* image, 
//								 typename TImage::PixelType overlayVal,
//								 double radiusMultipFactor,
//								 bool clipTails = false) const
//		{
//			typedef TImage																						ImageType;
//			typedef ShapedNeighborhoodIterator<ImageType>							NeighborhoodIteratorType;		
//			
//			typedef Image<bool, ImageType::ImageDimension>						TailMapImageType;
//			typedef ShapedNeighborhoodIterator<TailMapImageType>			TailMapNeighborhoodIteratorType;
//			
//			ImageRegionType region = image->GetBufferedRegion();
//			NeighborhoodIteratorType edgeNeighIter;
//			
//			// Creating a circular shaped neighborhood iterator creator.
//			CircularShapedNeighborhoodIterCreator<NeighborhoodIteratorType> iterCreator;
//			iterCreator.SetRadiusInImageCoordinates( false );		// radius values of the path are stored in world coords.
//			iterCreator.SetRegion( region );
//			iterCreator.SetImage( image );
//			
//			
//			typename TailMapImageType::Pointer tailMapImage = 0;
//			if( clipTails )
//			{
//				if( m_VertexList->Size() <= 1 )
//				{
//					itkWarningMacro(<<"Clip Tails flag is turned on but the path has less than two vertices. "
//													<<"For overlaying with clipping, at least two vertices are required "
//													<<"to be in the path so that path direction at the start and end of the "
//													<<"path can be determined. Skipping the tail clipping step ...");
//				}
//				else
//				{
//					// Compute the bounding box of the tubular path after 
//					// rescaling the radius values with the given factor.
//					Pointer clone = this->Clone();
//					clone->ScaleAndShiftRadii(radiusMultipFactor);
//					ImageRegionType boundingRegion = 
//					clone->ComputeBoundingImageRegionWithRadius(image->GetSpacing());
//					boundingRegion.Crop(region);
//					
//					// Create a boolean map image for the tail regions.
//					tailMapImage = TailMapImageType::New();
//					tailMapImage->SetLargestPossibleRegion( image->GetLargestPossibleRegion() );
//					tailMapImage->SetBufferedRegion( boundingRegion );
//					tailMapImage->SetRequestedRegion( boundingRegion );
//					tailMapImage->SetSpacing( image->GetSpacing() );
//					tailMapImage->SetOrigin( image->GetOrigin() );
//					tailMapImage->SetDirection(image->GetDirection() );
//					tailMapImage->SetNumberOfComponentsPerPixel(image->GetNumberOfComponentsPerPixel() );
//					tailMapImage->Allocate();
//					tailMapImage->FillBuffer( true );
//
//					
//					// Creating a circular shaped neighborhood iterator creator only for the tail regions.
//					CircularShapedNeighborhoodIterCreator<TailMapNeighborhoodIteratorType> tailMapIterCreator;
//					tailMapIterCreator.SetRadiusInImageCoordinates( false );		// radius values of the path are stored in world coords.
//					tailMapIterCreator.SetRegion( boundingRegion );
//					tailMapIterCreator.SetImage( tailMapImage );
//					tailMapIterCreator.SetGenerateHalfCircularNeigh( true );
//					
//					// Create a half-circular neighborhood iterator for the start pixel 
//					// of the path.
//					TailMapNeighborhoodIteratorType tailMapEdgeNeighIter;
//					VertexType firstVertex = m_VertexList->ElementAt(1);
//					VertexType secondVertex = m_VertexList->ElementAt(0);
//					
//					tailMapIterCreator.SetRadius(this->EvaluateRadius(this->StartOfInput()) * 
//																										 radiusMultipFactor);
//					tailMapIterCreator.SetHalfCircleDirection( secondVertex - firstVertex );
//					tailMapIterCreator.GenerateIterator(tailMapEdgeNeighIter);
//					
//					typename ImageType::IndexType startIndex;
//					startIndex.CopyWithRound( secondVertex );
//					if( !(region.IsInside( startIndex )) )
//					{
//						itkExceptionMacro(<<"Start index of the path " << startIndex 
//															<<" is found to be outside the image buffered region "
//															<<region);
//					}
//					tailMapEdgeNeighIter.SetLocation( startIndex );
//					
//					// Now, mark the half-circular region for the start pixel 
//					// of the path.
//					this->OverlayFromNeighIter(false, 
//																		 tailMapEdgeNeighIter);
//					
//					
//					// Create a half-circular neighborhood iterator for the end pixel 
//					// of the path.
//					firstVertex = m_VertexList->ElementAt(m_VertexList->Size() - 2);
//					secondVertex = m_VertexList->ElementAt(m_VertexList->Size() - 1);
//					
//					tailMapIterCreator.SetRadius(this->EvaluateRadius(this->EndOfInput()) * 
//																			 radiusMultipFactor);
//					tailMapIterCreator.SetHalfCircleDirection( secondVertex - firstVertex );
//					tailMapIterCreator.GenerateIterator(tailMapEdgeNeighIter);
//					
//					typename ImageType::IndexType endIndex;
//					endIndex.CopyWithRound( secondVertex );
//					if( !(region.IsInside( endIndex )) )
//					{
//						itkExceptionMacro(<<"End index of the path " << endIndex 
//															<<" is found to be outside the path image buffered region "
//															<<region);
//					}
//					tailMapEdgeNeighIter.SetLocation( endIndex );
//					
//					// Now, mark the half-circular region for the end pixel 
//					// of the path.
//					this->OverlayFromNeighIter(false, 
//																		 tailMapEdgeNeighIter);
//					
//					// Finally, be sure that the start and end indices 
//					// of the path centerline are excluded from the tail map.
//					tailMapImage->SetPixel( startIndex, true );
//					tailMapImage->SetPixel( endIndex, true );
//				}				
//			}
//
//			
//			// Traverse each pixel along the path and for 
//			// each such pixel traverse its circular neighborhood.
//			PathConstIterator<ImageType, Self> iter(image, this );
//			for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
//			{				
//				// Check if the index is inside the image region.
//				if( !(region.IsInside( iter.GetIndex() )) )
//				{
//					itkExceptionMacro(<<"A path point index " << iter.GetIndex() 
//														<<" is found to be outside the image buffered region "
//														<<region);
//				}
//				
//				// Create a new neighborhood iterator with the radius attribute taken 
//				// from the path point of the tubularity graph edge.
//				iterCreator.SetRadius(this->EvaluateRadius(iter.GetPathPosition()) * 
//															radiusMultipFactor);
//				iterCreator.GenerateIterator(edgeNeighIter);
//				
//				edgeNeighIter.SetLocation( iter.GetIndex() );
//				
//				this->OverlayFromNeighIter(overlayVal, 
//																	 edgeNeighIter, 
//																	 tailMapImage.GetPointer());
//			}
//		}
		
		
		// TODO Maybe create a path to path filter from this function!
		template<class TImage>
		InputType SplitByMinLinePairFittingError(Self* pathLine1,
																										 Self* pathLine2,
																										 const TImage* image,
																										 double& errorLine1 = 0, 
																										 double& errorLine2 = 0,
																										 double& cosOfAngleBtwLines = 0)
		{
			// Get the path length (number of points on the path).
			unsigned long length = m_VertexList->Size();
		
			if( length < 3 )
			{
				itkExceptionMacro("Path length must be greater than 2." 
													<<"That is, there must be at least three points in a path.");
			}
			
			
			if( pathLine1 == 0 || pathLine2 == 0 )
			{
				itkExceptionMacro("One of the input path objects is null.");
			}
			
			
			Array<double> cosArray(length);
			cosArray.Fill( 1.0 );
			PointType startPoint;
			PointType endPoint;
			PointType centerPoint;
			PointType point;
			VectorType line1;
			VectorType line2;
			double minError = NumericTraits<double>::max();
			VectorType pointVector;
			InputType startOfInput = this->StartOfInput();
			InputType endOfInput = this->EndOfInput();
			unsigned long centerPointIndex = 1;
			
			image->TransformContinuousIndexToPhysicalPoint(this->Evaluate(startOfInput), startPoint);
			image->TransformContinuousIndexToPhysicalPoint(this->Evaluate(endOfInput), endPoint);
			for(unsigned long i = 1; i < length - 1; i++)
			{
				image->TransformContinuousIndexToPhysicalPoint
				(this->Evaluate(static_cast<InputType>(i)+startOfInput), centerPoint);
				
				line1 = centerPoint - startPoint;
				line2 = endPoint - centerPoint;
				
				line1 /= line1.GetNorm();
				line2 /= line2.GetNorm();
				
				// dot product of unit vectors.
				cosArray[i] = line1 * line2;
			}
			
			for(unsigned long i = 1; i < length - 1; i++)
			{	
				// Takes as candidates only those points
				// that are local maxima of curvature.
				if((cosArray[i] < cosArray[i-1]) && 
					 (cosArray[i] < cosArray[i+1]))
				{
					// Compute the line fitting error.
					image->TransformContinuousIndexToPhysicalPoint
					(this->Evaluate(static_cast<InputType>(i)+startOfInput), centerPoint);
					
					line1 = centerPoint - startPoint;
					line2 = endPoint - centerPoint;
					
					line1 /= line1.GetNorm();
					line2 /= line2.GetNorm();
					
					double currentErrorLine1 = 0;
					double pointError;
					for(unsigned long j = 1; j <= i; j++)
					{
						image->TransformContinuousIndexToPhysicalPoint
						(this->Evaluate(static_cast<InputType>(j)+startOfInput), point);
						
						pointVector = point - startPoint;
						
						// || v - (v.u)xu || : magnitude of the component of the vector v perpendicular to the unit vector u. 
						pointError = (pointVector - ((pointVector * line1) * line1)).GetNorm();
						
						currentErrorLine1 += pointError * pointError;
					}
					double currentErrorLine2 = 0;
					for(unsigned long j = i; j < length-1; j++)
					{
						image->TransformContinuousIndexToPhysicalPoint
						(this->Evaluate(static_cast<InputType>(j)+startOfInput), point);
						
						pointVector = endPoint - point;
						
						// || v - (v.u)xu || : magnitude of the component of the vector v perpendicular to the unit vector u. 								
						pointError = (pointVector - ((pointVector * line2) * line2)).GetNorm();
						
						currentErrorLine2 += pointError * pointError;
					}
					
					if( minError > (currentErrorLine1 + currentErrorLine2) )
					{
						minError  = currentErrorLine1 + currentErrorLine2;
						errorLine1 = currentErrorLine1;
						errorLine2 = currentErrorLine2;
						cosOfAngleBtwLines = cosArray[i];
						
						centerPointIndex = i;
					}
				}				
			}
			
			// Create two paths, each corresponding to a line.
			InputType splitPoint = ((InputType)centerPointIndex) + this->StartOfInput();
			this->Split(splitPoint, 
									pathLine1, 
									pathLine2 );
			
			return splitPoint;
		}
		
		
		// TODO Maybe create a path to path filter from this function!
		// The input image pixel type must be scalar. 
		template<class TImage>
		InputType SplitByMaxAbsMeanPixelValDiff(Self* path1,
																										Self* path2,
																										const TImage* image,
																										double& mean1 = 0.0,
																										double& mean2 = 0.0) const
		{
			typedef  PathConstIterator<TImage, Self>				IterType;
			typedef  typename TImage::PixelType							PixelType;
			
			// Create the integral path from pixel values.
			std::vector<double> integralPelVals;
			std::vector<InputType> pathInputVals;
			IterType iter(image, this);
			iter.GoToBegin();
			integralPelVals.push_back( static_cast<double>(iter.Get()) );
			pathInputVals.push_back( iter.GetPathPosition() );
			++iter;
			for(unsigned long i = 0; !iter.IsAtEnd(); ++iter, i++)
			{
				integralPelVals.push_back(integralPelVals[i] + 
																	 static_cast<double>(iter.Get()) );
				
				pathInputVals.push_back( iter.GetPathPosition() );
			}
			
			// Find the maximum absolute difference of mean values.
			double maxAbsMeanDiff = -NumericTraits<double>::max();
			InputType maxAbsMeanDiffPoint = 0;
			unsigned long maxIndex = 0; 
			double pathMean1;
			double pathMean2;
			double pathSum =  integralPelVals.back();
			unsigned long noOfPels = integralPelVals.size();			
			if( noOfPels < 3 )
			{
				itkExceptionMacro("Number of pixels on the path must be greater than 2.");	// TODO Does this make sense? Maybe we should return a success flag.
			}
			for(unsigned long i = 1; i < noOfPels-1; i++)
			{
				pathMean1 = integralPelVals[i] / (double(i + 1));
				pathMean2 = (pathSum - integralPelVals[i]) / (double(noOfPels - i - 1));
				
				if( maxAbsMeanDiff < vcl_fabs(pathMean1 - pathMean2) )
				{
					maxAbsMeanDiff = vcl_fabs(pathMean1 - pathMean2);
					maxAbsMeanDiffPoint = pathInputVals[i];
					maxIndex = i;
				}
			}

			this->Split(maxAbsMeanDiffPoint,
									path1,
									path2 );
			
			mean1 = integralPelVals[maxIndex] / (double(maxIndex + 1));
			mean2 = (pathSum - integralPelVals[maxIndex]) / (double(noOfPels - maxIndex - 1));
			
			return maxAbsMeanDiffPoint;
		}
		
		
		// TODO Delete this function. Used only for debugging.
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		template<class TImage>
		void OverlayAndWrite(const TImage* image, 
												 typename TImage::PixelType pixelStartVal, 
												 typename TImage::PixelType pixelEndVal,
												 std::string fileName)
		{

			typedef itk::ImageFileWriter< TImage >																								WriterFilterType;
			typedef itk::ImageDuplicator< TImage > DuplicatorType;
			
			typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
			duplicator->SetInputImage( image );
			duplicator->Update();
			typename TImage::Pointer duplicatorOutput = duplicator->GetOutput();
			duplicatorOutput->DisconnectPipeline();
			
			PathConstIterator<TImage, Self> pathIter(duplicatorOutput, this );
			unsigned long pathLength = 0;
			for( pathIter.GoToBegin(); !pathIter.IsAtEnd(); ++pathIter )
			{	
				pathLength++;
			}
			typename TImage::PixelType increment = (pixelEndVal - pixelStartVal) / ((typename TImage::PixelType)pathLength);
			typename TImage::PixelType currentPixelVal = pixelStartVal;
			for( pathIter.GoToBegin(); !pathIter.IsAtEnd(); ++pathIter )
			{		
				duplicatorOutput->SetPixel( pathIter.GetIndex(), currentPixelVal );
				currentPixelVal = currentPixelVal + increment;
			}
			
			typename WriterFilterType::Pointer writer = WriterFilterType::New();
			writer->SetFileName( fileName );
			writer->SetInput( duplicatorOutput );
			try
			{
				writer->Update();
			}
			catch( itk::ExceptionObject & excp )
			{
				std::cerr << "Problem encountered while writing ";
				std::cerr << " image file : " << fileName << std::endl;
				std::cerr << excp << std::endl;
			}
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		
		// TODO Maybe create a path to path filter from this function!
		// The input image pixel type must be scalar. 
		template<class TImage>
		std::pair<InputType,InputType>  ExtractPathSegmentByMaxAbsMeanPixelValDiff(Self* path,
																																							 unsigned int minNoOfPixels,
																																							 const TImage* image,
																																							 double& mean1 = 0.0,
																																							 double& mean2 = 0.0) const
		{
			typedef  PathConstIterator<TImage, Self>				IterType;
			typedef  typename TImage::PixelType							PixelType;
			
			// Create the integral path from pixel values.
			std::vector<double> integralPelVals;
			std::vector<InputType> pathInputVals;
			IterType iter(image, this);
			iter.GoToBegin();
			integralPelVals.push_back( static_cast<double>(iter.Get()) );
			pathInputVals.push_back( iter.GetPathPosition() );
			++iter;
			for(unsigned long i = 0; !iter.IsAtEnd(); ++iter, i++)
			{
				integralPelVals.push_back(integralPelVals[i] + 
																	 static_cast<double>(iter.Get()) );
				
				pathInputVals.push_back( iter.GetPathPosition() );
			}
			
			if( integralPelVals.size() <= minNoOfPixels )
			{
				path->Copy(this);
				mean1 = integralPelVals.back() / (double(integralPelVals.size()));
				mean2 = mean1;
				
				return std::make_pair(pathInputVals.front(), pathInputVals.back());
			}
			
			// Find the maximum absolute difference of mean values.
			double maxAbsMeanDiff = -NumericTraits<double>::max();
			InputType maxAbsMeanDiffStartPoint = 0;
			InputType maxAbsMeanDiffEndPoint = 0;
			double pathSegmentSum;
			double pathMean1;
			double pathMean2;
			double pathSum =  integralPelVals.back();
			unsigned long noOfPels = integralPelVals.size();			
			if( noOfPels < 3 )
			{
				itkExceptionMacro("Number of pixels on the path must be greater than 2.");	// TODO Does this make sense? Maybe we should return a success flag.
			}
			
			unsigned long endPathLength;
			if( noOfPels > 2 * minNoOfPixels )
			{
				endPathLength = minNoOfPixels + (noOfPels - 2 * minNoOfPixels) + 1;
			}
			else
			{
				endPathLength = minNoOfPixels + 1;	// only try the minNoOfPixels as a path length!
			}
			
			for(unsigned long segmentLength = minNoOfPixels; 
					segmentLength < endPathLength; 
					segmentLength++)
			{
				for(unsigned long i = 0; i + segmentLength <= noOfPels; i++)
				{				
					if( i != 0 )
					{
						pathSegmentSum = (integralPelVals[i+segmentLength-1] - integralPelVals[i-1]);
					}
					else
					{
						pathSegmentSum = integralPelVals[i+segmentLength-1];
					}
					pathMean1 = pathSegmentSum / (double(segmentLength));					pathMean2 = (pathSum - pathSegmentSum) / (double(noOfPels - segmentLength));
					
					if( maxAbsMeanDiff < (pathMean2 - pathMean1) )
					{
						maxAbsMeanDiff = (pathMean2 - pathMean1);					// Store the path segment that has the smallest mean.
						maxAbsMeanDiffStartPoint = pathInputVals[i];
						maxAbsMeanDiffEndPoint = pathInputVals[i+segmentLength-1];
						mean1 = pathMean1;
						mean2 = pathMean2;
					}
				}
			}
			
			
			this->SubPath(maxAbsMeanDiffStartPoint,
										maxAbsMeanDiffEndPoint,
										path);
			
//			// TODO
//			IndexType startIndex;
//			startIndex[0] = 381;
//			startIndex[1] = 201;
//			IndexType endIndex;
//			endIndex[0] = 360;
//			endIndex[1] = 182;
//			bool istiyonMu = (this->EvaluateToIndex(this->StartOfInput()) == startIndex) && 
//											 (this->EvaluateToIndex(this->EndOfInput()) == endIndex);
//			if( istiyonMu )
//			{
//				this->OverlayAndWrite(image, 0, -5, "/Users/engin/Desktop/EsasYavru_1.nrrd");
//				path->OverlayAndWrite(image, 0, -5, "/Users/engin/Desktop/MinikYavru_1.nrrd");
//			}
			
//			startIndex[0] = 396;
//			startIndex[1] = 169;
//			endIndex[0] = 409;
//			endIndex[1] = 147;
//			istiyonMu = (this->EvaluateToIndex(this->StartOfInput()) == startIndex) && 
//			(this->EvaluateToIndex(this->EndOfInput()) == endIndex);
//			if( istiyonMu )
//			{
//				this->OverlayAndWrite(image, 0, -5, "/Users/engin/Desktop/EsasYavru_2.nrrd");
//				path->OverlayAndWrite(image, 0, -5, "/Users/engin/Desktop/MinikYavru_2.nrrd");
//			}
			
			
			return std::make_pair(maxAbsMeanDiffStartPoint, maxAbsMeanDiffEndPoint);
		}
		
		
		// TODO Create a path to path filter from this function!
		template<class TImage>
		void  ExtractPathSegmentsAsymmetric(std::vector<Pointer>& pathSegmentArray,
																				double pathSegmentLengthInWorldCoords,
																				double stepSizeInWorldCoords,
																				const TImage* image) const
		{
			pathSegmentArray.clear();
			
			double traversedDistance; // dummy
			
			InputType pathSegmentEndPoint;
			for(InputType i = this->StartOfInput(); 
					(this->EndOfInput() - i)  > m_Epsilon;)
			{
				// Get the end point of the current path segment.
				if( !TraverseDistanceInWorldCoords(i,
																					 image,
																					 pathSegmentLengthInWorldCoords,
																					 true,
																					 pathSegmentEndPoint,
																					 traversedDistance) )
				{
					// At this point, we reach end of the path!
					
					// Check if we have extracted any segment before. 
					// This happens, when the total path 
					// length is smaller than the segment length.
					// If not, add the whole path to the list.
					if( pathSegmentArray.empty() )
					{
						Pointer pathSegment = Self::New();
						SubPath(i, pathSegmentEndPoint, pathSegment);
						pathSegmentArray.push_back(pathSegment);
					}
					else
					{					
						// If segment length is larger than the step size, try to span the 
						// tail of the path so that we cover the path entirely, while 
						// keeping the length of the last segment at 
						// pathSegmentLengthInWorldCoords.
						if( pathSegmentLengthInWorldCoords - stepSizeInWorldCoords > -m_Epsilon )
						{
							// Traverse in the backward direction starting from the end
							// of the path.
							TraverseDistanceInWorldCoords(this->EndOfInput(),
																						image,
																						pathSegmentLengthInWorldCoords,
																						false,
																						pathSegmentEndPoint,
																						traversedDistance);
							
							Pointer pathSegment = Self::New();
							SubPath(pathSegmentEndPoint, this->EndOfInput(), pathSegment);
							pathSegmentArray.push_back(pathSegment);
						}
					}
					
					break;
				}
				
				// Create the path segment and add it to the list.
				Pointer pathSegment = Self::New();
				SubPath(i, pathSegmentEndPoint, pathSegment);
				pathSegmentArray.push_back(pathSegment);
				
				// Get the start point for the next segment by traversing step 
				// size along the path.
				if( !TraverseDistanceInWorldCoords(i,
																					 image,
																					 stepSizeInWorldCoords,
																					 true,
																					 i,
																					 traversedDistance) )
				{
					break;
				}
			}
			
		}
		
		
		
		// TODO Create a path to path filter from this function!
		template<class TImage>
		void  ExtractPathSegmentsSymmetric(std::vector<Pointer>& pathSegmentArray,
																			 double pathSegmentLengthInWorldCoords,
																			 double stepSizeInWorldCoords,
																			 const TImage* image) const
		{
			pathSegmentArray.clear();
			
			double pathLength = this->GetLengthInWorldCoords(image);
			
			// Check if the total path length is smaller than the 
			// segment length. If so, add the whole path to the list.
			if( pathLength - pathSegmentLengthInWorldCoords < m_Epsilon )
			{
				Pointer pathSegment = this->Clone();
				pathSegmentArray.push_back(pathSegment);
				return;
			}
			
			double traversedDistance; // dummy
			
			// Determine the start point for starting from the middle and going in the 
			// backward direction. The first segment will be centralized at the middle 
			// of the path.
			InputType i = 0.5*(this->EndOfInput() - this->StartOfInput());
			TraverseDistanceInWorldCoords(i,
																		image,
																		0.5 * pathSegmentLengthInWorldCoords,
																		true,
																		i,
																		traversedDistance);
			
			InputType pathSegmentStartPoint;
			InputType pathSegmentEndPoint;
			while((i - this->StartOfInput())  > m_Epsilon)
			{
				// Get the start point of the current path segment.
				if( !TraverseDistanceInWorldCoords(i,
																					 image,
																					 pathSegmentLengthInWorldCoords,
																					 false,
																					 pathSegmentStartPoint,
																					 traversedDistance) )
				{
					// At this point, we reach end of the path!
								
					// If segment length is larger than the step size, try to span the 
					// tail of the path so that we cover the path entirely, while 
					// keeping the length of the last segment at 
					// pathSegmentLengthInWorldCoords.
					if( pathSegmentLengthInWorldCoords - stepSizeInWorldCoords > -m_Epsilon )
					{
						// Traverse in the forward direction starting from the start
						// of the path.
						TraverseDistanceInWorldCoords(this->StartOfInput(),
																					image,
																					pathSegmentLengthInWorldCoords,
																					true,
																					pathSegmentEndPoint,
																					traversedDistance);
						
						Pointer pathSegment = Self::New();
						SubPath(this->StartOfInput(), pathSegmentEndPoint, pathSegment);
						pathSegmentArray.push_back(pathSegment);
					}
					
					break;
				}
				
				// Create the path segment and add it to the list.
				Pointer pathSegment = Self::New();
				SubPath(pathSegmentStartPoint, i, pathSegment);
				pathSegmentArray.push_back(pathSegment);
				
				// Get the end point for the next segment by traversing step 
				// size along the path in the backward direction.
				if( !TraverseDistanceInWorldCoords(i,
																					 image,
																					 stepSizeInWorldCoords,
																					 false,
																					 i,
																					 traversedDistance) )
				{
					break;
				}
			}
			// Reverse the array to get the correct ordering of segments.
			std::reverse(pathSegmentArray.begin(), pathSegmentArray.end());
			
			
			// Determine the start point for starting from the middle and going in the 
			// forward direction. The first segment will NOT be centralized at the middle 
			// of the path, since such a segment has already been extracted before.
			i = 0.5*(this->EndOfInput() - this->StartOfInput());
			TraverseDistanceInWorldCoords(i,
																		image,
																		0.5 * pathSegmentLengthInWorldCoords,
																		false,
																		i,
																		traversedDistance);
			
			if( TraverseDistanceInWorldCoords(i,
																				image,
																				stepSizeInWorldCoords,
																				true,
																				i,
																				traversedDistance) )
			{
				// Now traverse in the forward direction.
				while((this->EndOfInput() - i)  > m_Epsilon)
				{
					// Get the end point of the current path segment.
					if( !TraverseDistanceInWorldCoords(i,
																						 image,
																						 pathSegmentLengthInWorldCoords,
																						 true,
																						 pathSegmentEndPoint,
																						 traversedDistance) )
					{
						// At this point, we reach end of the path!
						
						// If segment length is larger than the step size, try to span the 
						// tail of the path so that we cover the path entirely, while 
						// keeping the length of the last segment at 
						// pathSegmentLengthInWorldCoords.
						if( pathSegmentLengthInWorldCoords - stepSizeInWorldCoords > -m_Epsilon )
						{
							// Traverse in the backward direction starting from the end
							// of the path.
							TraverseDistanceInWorldCoords(this->EndOfInput(),
																						image,
																						pathSegmentLengthInWorldCoords,
																						false,
																						pathSegmentStartPoint,
																						traversedDistance);
							
							Pointer pathSegment = Self::New();
							SubPath(pathSegmentStartPoint, this->EndOfInput(), pathSegment);
							pathSegmentArray.push_back(pathSegment);
						}
						
						break;
					}
					
					// Create the path segment and add it to the list.
					Pointer pathSegment = Self::New();
					SubPath(i, pathSegmentEndPoint, pathSegment);
					pathSegmentArray.push_back(pathSegment);
					
					// Get the start point for the next segment by traversing step 
					// size along the path in the forward direction.
					if( !TraverseDistanceInWorldCoords(i,
																						 image,
																						 stepSizeInWorldCoords,
																						 true,
																						 i,
																						 traversedDistance) )
					{
						break;
					}
				}
			}
			
		}
		
		
		
		
		
		// Returns true if the specified distance has been sucessfully traversed
		// without prematurely hitting the end point of the path.
		template<class TImage>
		bool TraverseDistanceInWorldCoords(InputType startPoint,
																			 const TImage* image,
																			 double distInWorldCoords,
																			 bool traverseForward,
																			 InputType& endPoint,
																			 double& traversedDistance) const
		{
			if( distInWorldCoords <= 0 )
			{
				itkExceptionMacro("Distance to traverse can not be smaller than or equal to zero.");
			}
			
			traversedDistance = 0.0;
			InputType increment;
			if( traverseForward )
			{
				increment = 1.0;
				
				if( startPoint < this->StartOfInput() )
				{
					startPoint = this->StartOfInput();
				}
				
				if( startPoint >= this->EndOfInput() )
				{
					endPoint = this->EndOfInput();
					return false;
				}
			}
			else 
			{
				increment = -1.0;
				
				if( startPoint > this->EndOfInput() )
				{
					startPoint = this->EndOfInput();
				}
				
				if( startPoint <= this->StartOfInput() )
				{
					endPoint = this->StartOfInput();
					return false;
				}
			}

			InputType currentPoint = startPoint;
			InputType nextPoint;
			if( traverseForward )
			{
				nextPoint = (InputType)(((unsigned long)startPoint) + 1);
			}
			else
			{
				nextPoint = (InputType)(((unsigned long)startPoint));
			}
			while(vcl_fabs(distInWorldCoords - traversedDistance) > m_Epsilon)
			{
				double currentDistBtwSuccPoints = GetDistanceInWorldCoords(currentPoint, 
																																	 nextPoint, 
																																	 image);
				
				// Did we pass the point to stop.
				if(traversedDistance + currentDistBtwSuccPoints > distInWorldCoords)
				{
					// Find the portion of the line to traverse.
					endPoint = currentPoint + (nextPoint - currentPoint) * 
					((distInWorldCoords - traversedDistance)/currentDistBtwSuccPoints);
					
					traversedDistance = distInWorldCoords;
					
					return true;
				}
				else
				{
					traversedDistance += currentDistBtwSuccPoints;
				}
				
				endPoint = nextPoint;
				
				// Did we reach end of the path.
				if( traverseForward )
				{
					if( vcl_fabs(nextPoint - this->EndOfInput()) < m_Epsilon )
					{
						return false;
					}
				}
				else
				{
					if( vcl_fabs(nextPoint - this->StartOfInput()) < m_Epsilon )
					{
						return false;
					}
				}
				
				
				currentPoint = nextPoint;
				nextPoint += increment;
			}
			
			return true;
		}
		
		virtual OffsetType IncrementInput(InputType & input) const;
		
		virtual inline InputType EndOfInput() const
    {
			return m_VertexList->Size() - 1;
    }
		
		virtual void Initialize(void)
		{
			m_VertexList->Initialize();
			m_RadiusList.clear();
		}
	
		itkGetConstObjectMacro( VertexList, VertexListType );
		
		virtual unsigned long GetNumberOfVertices()
		{
			return static_cast<unsigned long>(m_VertexList->Size());
		}
		
		
		// TODO Add an EraseVertex() function.
		
		virtual bool InsertVertex(typename VertexListType::ElementIdentifier vertexIndex, 
															const VertexType& vertex,
															RadiusType radius)
		{
			if( (vertexIndex < 0) && (vertexIndex > m_VertexList->Size()) )
			{
				itkExceptionMacro(<<"Insertion index is out of bounds.");
			}
			
			// Check if the given vertex is close to its neighbors in the given 
			// position of the list and if so, do not insert it.
			if( vertexIndex > 0 )
			{
				const VertexType& prevVertex = m_VertexList->ElementAt(vertexIndex - 1);
				
				if( vertex.SquaredEuclideanDistanceTo(prevVertex) <  m_Epsilon )
				{
					return false;
				}
			}
			if( vertexIndex < m_VertexList->Size() )
			{
				const VertexType& nextVertex = m_VertexList->ElementAt(vertexIndex);
				
				if( vertex.SquaredEuclideanDistanceTo(nextVertex) <  m_Epsilon )
				{
					return  false;
				}
			}
			
			// TODO The below call to CastToSTLContainer() requires dynamic 
			// cast, which is costly. You can avoid this by replacing the 
			// VectorContainer type for the vertex list by std::vector<VertexType>.
			typename std::vector<VertexType>::size_type insertLoc = vertexIndex;
			std::vector<VertexType>& vertexList = m_VertexList->CastToSTLContainer();
			vertexList.insert(vertexList.begin() + insertLoc, vertex);
			m_RadiusList.insert(m_RadiusList.begin() + insertLoc, radius);
			
			return true;
		}
		
		virtual VertexType GetVertex(typename VertexListType::ElementIdentifier vertexIndex) const
		{
			VertexIndexBoundCheck(vertexIndex);
			return m_VertexList->ElementAt(vertexIndex);
		}
		
		virtual bool SetVertex(typename VertexListType::ElementIdentifier vertexIndex, 
													 const VertexType& vertex)
		{
			VertexIndexBoundCheck(vertexIndex);
		
			// Check if the given vertex is close to its neighbors in the given 
			// position of the list and if so, do not set it.
			if( vertexIndex > 0 )
			{
				const VertexType& prevVertex = m_VertexList->ElementAt(vertexIndex - 1);
				
				if( vertex.SquaredEuclideanDistanceTo(prevVertex) <  m_Epsilon )
				{
					return false;
				}
			}
			if( (vertexIndex + 1) < m_VertexList->Size() )
			{
				const VertexType& nextVertex = m_VertexList->ElementAt(vertexIndex + 1);
				
				if( vertex.SquaredEuclideanDistanceTo(nextVertex) <  m_Epsilon )
				{
					return false;
				}
			}
			
			m_VertexList->ElementAt(vertexIndex) = vertex;
			
			return true;
		}
		
		virtual bool SetVertex(typename VertexListType::ElementIdentifier vertexIndex, 
														const VertexType& vertex, 
														RadiusType radius)
		{
			if( !SetVertex(vertexIndex, vertex) )
			{
				return false;
			}
			
			SetVertexRadius(vertexIndex, radius);
			return true;
		}
		
		virtual const RadiusListType& GetRadiusList() const
		{
			return m_RadiusList;
		}
		
		virtual RadiusType GetVertexRadius(typename RadiusListType::size_type vertexIndex) const
		{
			VertexIndexBoundCheck(vertexIndex);
			return m_RadiusList[vertexIndex];
		}
		
		virtual void SetVertexRadius(typename RadiusListType::size_type vertexIndex, 
																 RadiusType radius)
		{
			VertexIndexBoundCheck(vertexIndex);
			m_RadiusList[vertexIndex] = radius;
		}
		
		virtual unsigned int GetEffectiveDimension() const
		{
			return m_EffectiveDimension;
		}
		
		virtual double GetEpsilon() const
		{
			return m_Epsilon;
		}
		virtual void SetEpsilon(double eps)
		{
			m_Epsilon = eps;
		}
		
		
		// TODO HACK Store the full feature vector for this path here! 
		itk::Array<float> m_FeatureVector;
		
	protected:
		PolyLineParametricPathExtended();
		virtual ~PolyLineParametricPathExtended(){}
		virtual void PrintSelf(std::ostream& os, Indent indent) const;
		
		virtual void VertexIndexBoundCheck(typename VertexListType::ElementIdentifier vertexIndex) const
		{
			if( (vertexIndex < 0) && (vertexIndex >= m_VertexList->Size()) )
			{
				itkExceptionMacro(<<"Vertex index is out of bounds.");
			}
		}
		
	private:
		PolyLineParametricPathExtended(const Self&); //purposely not implemented
		void operator=(const Self&);									//purposely not implemented
		

		template <class TPoint>
		double 
		GetEuclDistToLineSegment(const TPoint& lineVertex1, 
														 const TPoint& lineVertex2, 
														 const TPoint& vertex, 
														 InputType& minDistPoint) const;	// takes values between 0 (start vertex) and 1 (end vertex).
		
		
		// Check if a given point is in the tail region. 
		// A tail region for an end point of the path is 
		// represented as a hemisphere in world coordinate system.
		// No radius check is performed in this method.
		template <class TImage>
		bool
		IsPointInTailRegion(const PointType& point, 
												const InputType& closestPathInput,
												const TImage* image) const;
		
		/** Smooths the vertex locations and radius values. */
		template<class TImage>
		void Smooth(double avgWindowRadiusInWorldCoords,
								const TImage* image,
								bool smoothLocations,
								bool smoothRadii);
		
		
		template<class TImage>
		void OverlayFromNeighIter(typename TImage::PixelType overlayVal,
															ShapedNeighborhoodIterator<TImage>& iter,
															const Image<bool, TImage::ImageDimension>* binaryMap = 0) const;
		
		void ComputeBoundingImageRegionFromVertexList(const VertexListType* vertexList,
																									ImageRegionType& region) const;
			
		
		friend class boost::serialization::access;
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		
		template<class Archive>
		void save(Archive & ar, const unsigned int version) const
		{
			unsigned int dimension = Dimension;
			ar & dimension;
			
			const typename VertexListType::VectorContainerSizeType count(m_VertexList->Size());
			ar & BOOST_SERIALIZATION_NVP(count);
			for(unsigned int i = 0; i < count; i++)
			{
				const VertexType& cIndex= m_VertexList->ElementAt(i);
				RadiusType radi = m_RadiusList[i];
				
				for(unsigned int j = 0; j < Dimension; j++)
				{
					ar & cIndex[j];
				}
				
				ar & radi;
			}
		}
		
		template<class Archive>
		void load(Archive & ar, const unsigned int version)
		{
			this->Initialize();
			
			ar & m_EffectiveDimension;
			
			typename VertexListType::VectorContainerSizeType count;
			ar & BOOST_SERIALIZATION_NVP(count);
			
			for(unsigned int i = 0; i < count; i++)
			{
				RadiusType radi;
				VertexType cIndex; 
				cIndex.Fill(0);	// if Dimension > m_EffectiveDimension, set the rest of the elements to zero.
				for(unsigned int j = 0; j < std::min(VDimension, m_EffectiveDimension); j++)
				{
					ar & cIndex[j];
				}
				if( Dimension < m_EffectiveDimension ) // Do dummy reads  if needed.
				{
					typename VertexType::ValueType dummy;
					for(unsigned int i = Dimension; i < m_EffectiveDimension; i++)
					{
						ar & dummy;
					}
				}
				
				ar & radi;
				
				this->AddVertex( cIndex, radi );
			}
		}
		
		// Path point list in image coordinate system.
		VertexListPointer		m_VertexList;		
		
		// Path radius list in world coordinate syste.
		RadiusListType			m_RadiusList;
		
		
		unsigned int				m_EffectiveDimension; // Effective dimension read from the file.
		
		double							m_Epsilon;		// TODO INclude this in the serialization!
	};
	
}

#if ITK_TEMPLATE_TXX
# include "itkPolyLineParametricPathExtended.txx"
#endif

#endif
