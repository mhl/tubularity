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
				const ContinuousIndexType& endVertex = m_VertexList->ElementAt(m_VertexList->Size() - 1);
				
				if( vertex.EuclideanDistanceTo(endVertex) <  NumericTraits<double>::epsilon() )
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
		
		/** Extracts a subpath of this path */
		virtual void SubPath(InputType startPoint, InputType endPoint, Self* path)  const;
		
		/** Rounds the continuous indices this path traverses. */
		virtual void RoundIndices(  bool removeRepeatedVertices = true );
		
		/** Computes bounding image region of the path. */
		virtual ImageRegionType ComputeBoundingImageRegion() const;
		
		/** Compute minimum/maximum/mean radius along the path. */
		virtual RadiusType ComputeMinRadius() const;
		virtual RadiusType ComputeMaxRadius() const;
		virtual RadiusType ComputeMeanRadius() const;
		
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
		 */
		virtual unsigned long GetNoOfIntersectionPoints(const Self* path) const;
		
		/** Compute the sum of the pixels along the path. It is assumed that 
		 * the + operator for the type TSum is defined.
		 */
		template<class TSum, class TImage> 
		TSum PixelSum(const TImage* image ) const
		{
			typedef  PathConstIterator<TImage, Self>				IterType;
			typedef  TSum																		SumType;
			
			SumType pixelSum = NumericTraits<SumType>::ZeroValue();
			IterType iter(image, this);
			for(iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
			{
				pixelSum = pixelSum + static_cast<SumType>(iter.Get());
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
		
		/** Computes the distance between two points on the path */
		template<class TImage>
		double GetDistanceInWorldCoords(InputType startPoint, 
																						InputType endPoint, 
																						const TImage* image ) const
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
			
			if( vcl_fabs(endPoint - startPoint) < NumericTraits<InputType>::epsilon() )
			{
				return 0;
			}
			
			double dist = 0;
			ContinuousIndexType currentContIndx = Evaluate(startPoint);
			PointType currentPointLoc;
			image->TransformContinuousIndexToPhysicalPoint(currentContIndx, currentPointLoc);
			unsigned long nextVertexPoint = ((unsigned long)(startPoint - this->StartOfInput())) + 1;
			ContinuousIndexType nextContIndx = m_VertexList->ElementAt( nextVertexPoint );
			PointType nextPointLoc;
			while((endPoint - ((InputType)nextVertexPoint) - this->StartOfInput()) > 
						NumericTraits<InputType>::epsilon() )
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
		
		/** 
		 * Computes the length of the path in world coordinates along the pixels 
		 * of the image grid and not the lines of the path.
		 */
		template<class TImage>
		double GetLengthOnImageGridInWorldCoords(const TImage* image ) const
		{
			typedef itk::PathConstIterator<TImage, Self>			PathIterType;
			
			double distInWorldCoord = 0.0;
			PathIterType pathIter(image, this);
			
			bool isAtBegin = true;
			typename TImage::PointType prevWorldPoint;
			typename TImage::PointType currentWorldPoint;			
			
			for(pathIter.GoToBegin(); !pathIter.IsAtEnd(); ++pathIter)
			{
				image->TransformIndexToPhysicalPoint(pathIter.GetIndex(), currentWorldPoint);
				
				if( !isAtBegin )
				{
					// Compute the distance between the previous point and 
					// the current one in world coordinates.
					distInWorldCoord += currentWorldPoint.EuclideanDistanceTo(prevWorldPoint);
				}
				
				prevWorldPoint = currentWorldPoint;
				isAtBegin = false;
			}
			
			return distInWorldCoord;
		}
		
		template<class TImage>
		VectorType EvaluateDirectionVector(const InputType & input, 
																							 const TImage* image ) const
		{
			VectorType diffVector;
			if( m_VertexList->Size() < 2 )
			{
				diffVector.Fill(0.0);
				return diffVector;
			}
			
			bool nextExists = (input + 1.0) <= this->EndOfInput();
			bool prevExists = (input - 1.0) >= this->StartOfInput();	
			VertexType currentIndex = this->Evaluate( input );
			VertexType prevIndex;
			VertexType nextIndex;
			if( prevExists )
			{
				prevIndex = this->Evaluate( input - 1 );
			}
			else
			{
				prevIndex = currentIndex;
			}
			if( nextExists )
			{
				nextIndex = this->Evaluate( input + 1 );
			}
			else
			{
				nextIndex = currentIndex;
			}
			if( !prevExists && !nextExists  )
			{
				prevIndex = this->Evaluate( Math::Floor<InputType>(input) );
				nextIndex = this->Evaluate( Math::Ceil<InputType>(input) );
			}
			
			// Check if two points are the same.
			if(prevIndex.EuclideanDistanceTo(nextIndex) < NumericTraits<double>::epsilon())
			{
				// If all points are at the same location, then 
				// return zero vector.
				if(prevIndex.EuclideanDistanceTo(currentIndex) < NumericTraits<double>::epsilon())
				{
					diffVector.Fill(0.0);
					return diffVector;
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
		
		
		template<class TImage>
		void Overlay(TImage* image, 
								 typename TImage::PixelType overlayVal,
								 double radiusMultipFactor,
								 bool clipTails = false) const
		{
			typedef TImage																						ImageType;
			typedef ShapedNeighborhoodIterator<ImageType>							NeighborhoodIteratorType;		
			
			typedef Image<bool, ImageType::ImageDimension>						TailMapImageType;
			typedef ShapedNeighborhoodIterator<TailMapImageType>			TailMapNeighborhoodIteratorType;
			
			typename ImageType::RegionType region = image->GetLargestPossibleRegion();
			NeighborhoodIteratorType edgeNeighIter;
			
			// Creating a circular shaped neighborhood iterator creator.
			CircularShapedNeighborhoodIterCreator<NeighborhoodIteratorType> iterCreator;
			iterCreator.SetRadiusInImageCoordinates( false );		// radius values of the path are stored in world coords.
			iterCreator.SetRegion( region );
			iterCreator.SetImage( image );
			
			
			typename TailMapImageType::Pointer tailMapImage = 0;
			if( clipTails )
			{
				if( m_VertexList->Size() <= 1 )
				{
					itkWarningMacro(<<"Clip Tails flag is turned on but the path has less than two vertices. "
													<<"For overlaying with clipping, at least two vertices are required "
													<<"to be in the path so that path direction at the start and end of the "
													<<"path can be determined. Skipping the tail clipping step ...");
				}
				else
				{
					// Create a boolean map image for the tail regions.
					tailMapImage = TailMapImageType::New();
					tailMapImage->SetRegions( region );
					tailMapImage->SetSpacing( image->GetSpacing() );
					tailMapImage->SetOrigin( image->GetOrigin() );
					tailMapImage->SetDirection(image->GetDirection() );
					tailMapImage->SetNumberOfComponentsPerPixel(
																											image->GetNumberOfComponentsPerPixel() );
					tailMapImage->Allocate();
					tailMapImage->FillBuffer( true );

					
					// Creating a circular shaped neighborhood iterator creator only for the tail regions.
					CircularShapedNeighborhoodIterCreator<TailMapNeighborhoodIteratorType> tailMapIterCreator;
					tailMapIterCreator.SetRadiusInImageCoordinates( false );		// radius values of the path are stored in world coords.
					tailMapIterCreator.SetRegion( region );
					tailMapIterCreator.SetImage( tailMapImage );
					tailMapIterCreator.SetGenerateHalfCircularNeigh( true );
					
					// Create a half-circular neighborhood iterator for the start pixel 
					// of the path.
					TailMapNeighborhoodIteratorType tailMapEdgeNeighIter;
					VertexType firstVertex = m_VertexList->ElementAt(1);
					VertexType secondVertex = m_VertexList->ElementAt(0);
					
					tailMapIterCreator.SetRadius(this->EvaluateRadius(this->StartOfInput()) * 
																										 radiusMultipFactor);
					tailMapIterCreator.SetHalfCircleDirection( secondVertex - firstVertex );
					tailMapIterCreator.GenerateIterator(tailMapEdgeNeighIter);
					
					typename ImageType::IndexType startIndex;
					startIndex.CopyWithRound( secondVertex );
					if( !(region.IsInside( startIndex )) )
					{
						itkExceptionMacro(<<"Start index of the path " << startIndex 
															<<" is found to be outside the image largest possible region "
															<<region);
					}
					tailMapEdgeNeighIter.SetLocation( startIndex );
					
					// Now, mark the half-circular region for the start pixel 
					// of the path.
					this->OverlayFromNeighIter(region, 
																		 false, 
																		 tailMapEdgeNeighIter);
					
					
					// Create a half-circular neighborhood iterator for the end pixel 
					// of the path.
					firstVertex = m_VertexList->ElementAt(m_VertexList->Size() - 2);
					secondVertex = m_VertexList->ElementAt(m_VertexList->Size() - 1);
					
					tailMapIterCreator.SetRadius(this->EvaluateRadius(this->EndOfInput()) * 
																			 radiusMultipFactor);
					tailMapIterCreator.SetHalfCircleDirection( secondVertex - firstVertex );
					tailMapIterCreator.GenerateIterator(tailMapEdgeNeighIter);
					
					typename ImageType::IndexType endIndex;
					endIndex.CopyWithRound( secondVertex );
					if( !(region.IsInside( endIndex )) )
					{
						itkExceptionMacro(<<"End index of the path " << endIndex 
															<<" is found to be outside the image largest possible region "
															<<region);
					}
					tailMapEdgeNeighIter.SetLocation( endIndex );
					
					// Now, mark the half-circular region for the end pixel 
					// of the path.
					this->OverlayFromNeighIter(region, 
																		 false, 
																		 tailMapEdgeNeighIter);
					
					// Finally, be sure that the start and end indices 
					// of the path centerline are excluded from the tail map.
					tailMapImage->SetPixel( startIndex, true );
					tailMapImage->SetPixel( endIndex, true );
				}				
			}

			
			// Traverse each pixel along the path and for 
			// each such pixel traverse its circular neighborhood.
			PathConstIterator<ImageType, Self> iter(image, this );
			for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
			{				
				// Check if the index is inside the image region.
				if( !(region.IsInside( iter.GetIndex() )) )
				{
					itkExceptionMacro(<<"A path point index " << iter.GetIndex() 
														<<" is found to be outside the image largest possible region "
														<<region);
				}
				
				// Create a new neighborhood iterator with the radius attribute taken 
				// from the path point of the tubularity graph edge.
				iterCreator.SetRadius(this->EvaluateRadius(iter.GetPathPosition()) * 
															radiusMultipFactor);
				iterCreator.GenerateIterator(edgeNeighIter);
				
				edgeNeighIter.SetLocation( iter.GetIndex() );
				
				this->OverlayFromNeighIter(region, 
																	 overlayVal, 
																	 edgeNeighIter, 
																	 tailMapImage.GetPointer());
			}
		}
		
		
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
		
		
		// TODO Create a filter (Image to Image filter) out of this (Of course without writing to the file). Used only for debugging.
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
		void  ExtractPathSegments(std::vector<Pointer>& pathSegmentArray,
															double pathSegmentLengthInWorldCoords,
															double stepSizeInWorldCoords,
															const TImage* image) const
		{
			pathSegmentArray.clear();
			
			InputType pathSegmentEndPoint;
			for(InputType i = this->StartOfInput(); 
					this->EndOfInput() - i  < NumericTraits<InputType>::epsilon();)
			{
				// Get the end point of the current path segment.
				if( !TraverseDistanceInWorldCoords(i,
																					 image,
																					 pathSegmentLengthInWorldCoords,
																					 true,
																					 pathSegmentEndPoint) )
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
						// tail of the path so that we cover the path entirely.
						if( pathSegmentLengthInWorldCoords >= stepSizeInWorldCoords )
						{
							// Traverse in the backward direction starting from the end
							// of the path.
							TraverseDistanceInWorldCoords(this->EndOfInput(),
																						image,
																						pathSegmentLengthInWorldCoords,
																						false,
																						pathSegmentEndPoint);
							
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
																					 i) )
				{
					break;
				}
			}
			
		}
		
		// Returns true if the specified distance has been sucessfully traversed
		// without prematurely hitting the end point of the path.
		template<class TImage>
		bool TraverseDistanceInWorldCoords(InputType startPoint,
																			 const TImage* image,
																			 double distInWorldCoords,
																			 bool traverseForward = true,
																			 InputType& endPoint = 0,
																			 double& traversedDistance = 0) const
		{
			if( distInWorldCoords <= 0 )
			{
				itkExceptionMacro("Distance to traverse can not smaller than or equal to zero.");
			}
			
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
					traversedDistance = 0;
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
					traversedDistance = 0;
					endPoint = this->StartOfInput();
					return false;
				}
			}

			
			traversedDistance = 0;
			InputType currentPoint = startPoint;
			InputType nextPoint;
			if( traverseForward )
			{
				nextPoint = (InputType)(((unsigned long)startPoint) + 1);
			}
			else
			{
				if( IsIntegerInput(startPoint) )
				{
					nextPoint = (InputType)(((unsigned long)startPoint) - 1);
				}
				else
				{
					nextPoint = (InputType)(((unsigned long)startPoint));
				}
			}
			while(vcl_fabs(distInWorldCoords - traversedDistance) > NumericTraits<double>::epsilon())
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
					if( vcl_fabs(nextPoint - this->EndOfInput()) < NumericTraits<InputType>::epsilon() )
					{
						return false;
					}
				}
				else
				{
					if( vcl_fabs(nextPoint - this->StartOfInput()) < NumericTraits<InputType>::epsilon() )
					{
						return false;
					}
				}
				
				
				currentPoint = nextPoint;
				nextPoint += increment;
			}
			
			return true;
		}

		
		// TODO Call this function in all functions of this class. (search for the occurance of epsilon)
		virtual bool IsIntegerInput(InputType input) const
		{
			return vcl_fabs(input - Math::Round<InputType>(input)) < NumericTraits<InputType>::epsilon();
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
		
		virtual VertexType GetVertex(typename VertexListType::ElementIdentifier vertexIndex) const
		{
			if( (vertexIndex < 0) && (vertexIndex >= m_VertexList->Size()) )
			{
				itkExceptionMacro(<<"Vertex index is out of bounds.");
			}
			return m_VertexList->ElementAt(vertexIndex);
		}
		
		virtual void SetVertex(typename VertexListType::ElementIdentifier vertexIndex, 
													 VertexType vertex)
		{
			if( (vertexIndex < 0) && (vertexIndex >= m_VertexList->Size()) )
			{
				itkExceptionMacro(<<"Vertex index is out of bounds.");
			}
			m_VertexList->ElementAt(vertexIndex) = vertex;
		}
		
		virtual const RadiusListType& GetRadiusList() const
		{
			return m_RadiusList;
		}
		
		virtual RadiusType GetVertexRadius(typename RadiusListType::size_type vertexIndex) const
		{
			if( (vertexIndex < 0) && (vertexIndex >= m_RadiusList.size()) )
			{
				itkExceptionMacro(<<"Vertex index is out of bounds.");
			}
			return m_RadiusList[vertexIndex];
		}
		
		virtual void SetVertexRadius(typename RadiusListType::size_type vertexIndex, 
																 RadiusType radius)
		{
			if( (vertexIndex < 0) && (vertexIndex >= m_RadiusList.size()) )
			{
				itkExceptionMacro(<<"Vertex index is out of bounds.");
			}
			m_RadiusList[vertexIndex] = radius;
		}
		
		unsigned int GetEffectiveDimension()
		{
			return m_EffectiveDimension;
		}
		
		// TODO HACK Store the full feature vector for this path here! 
		itk::Array<float> m_FeatureVector;
		
	protected:
		PolyLineParametricPathExtended();
		virtual ~PolyLineParametricPathExtended(){}
		virtual void PrintSelf(std::ostream& os, Indent indent) const;
		
	private:
		PolyLineParametricPathExtended(const Self&); //purposely not implemented
		void operator=(const Self&);									//purposely not implemented
		

		template <class TPoint>
		double 
		GetEuclDistToLineSegment(const TPoint& lineVertex1, 
														 const TPoint& lineVertex2, 
														 const TPoint& vertex, 
														 InputType& minDistPoint) const;	// takes values between 0 (start vertex) and 1 (end vertex).
		
		
		template<class TImage>
		void OverlayFromNeighIter(typename TImage::RegionType& regionForBoundCheck, 
															typename TImage::PixelType overlayVal,
															ShapedNeighborhoodIterator<TImage>& iter,
															const Image<bool, TImage::ImageDimension>* binaryMap = 0) const;
		
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
	};
	
}

#if ITK_TEMPLATE_TXX
# include "itkPolyLineParametricPathExtended.txx"
#endif

#endif
