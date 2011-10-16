//////////////////////////////////////////////////////////////////////////////////
//																																							//
// Copyright (C) 2011 Engin Turetken																						//
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


#include "itkPolyLineParametricPathExtended.h"
#include "itkMath.h"
#include "itkBoundingBox.h"
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <set>

namespace itk
{
	
	template<unsigned int VDimension>
	typename PolyLineParametricPathExtended<VDimension>::OutputType
	PolyLineParametricPathExtended<VDimension>
	::Evaluate( const InputType & input ) const
	{
		OutputType output;
		VertexType vertex0;
		VertexType vertex1;
		double     fractionOfLineSegment;
		
		if(  input  <=  this->StartOfInput()  )
    {
			return m_VertexList->ElementAt(0); // the first vertex
		}
		else if(  input  >=  this->EndOfInput()  )
    {
			return m_VertexList->ElementAt(m_VertexList->Size() - 1); // the last vertex
    }
		
		vertex0 = m_VertexList->ElementAt( int(input) );
//		vertex1 = m_VertexList->ElementAt( int(input+1.0) );			// BUG 2  in the original PolyLineParametricPath!
		vertex1 = m_VertexList->ElementAt( ((int)input) + 1 );
		
		fractionOfLineSegment = input - int(input);
		
		// Optimized processing for 2D and 3D.
		if( VDimension == 2 )
		{
			output[0] = vertex0[0] + (vertex1[0]-vertex0[0]) * fractionOfLineSegment;
			output[1] = vertex0[1] + (vertex1[1]-vertex0[1]) * fractionOfLineSegment;
		}
		else if( VDimension == 3 )
		{
			output[0] = vertex0[0] + (vertex1[0]-vertex0[0]) * fractionOfLineSegment;
			output[1] = vertex0[1] + (vertex1[1]-vertex0[1]) * fractionOfLineSegment;
			output[2] = vertex0[2] + (vertex1[2]-vertex0[2]) * fractionOfLineSegment;
		}
		else
		{
			PointType outputPoint = vertex0 + (vertex1-vertex0)*fractionOfLineSegment;
			// For some stupid reason, there is no easy way to cast from a point to a
			// continuous index.
			for( unsigned int i=0; i<VDimension; i++ ) { output[i] = outputPoint[i]; }
		}
		
		return output;
	}
	
	// Performs linear interpolation of the radius values for a given point 
	// between two successive vertices of the path.
	template<unsigned int VDimension>
	typename PolyLineParametricPathExtended<VDimension>::RadiusType
	PolyLineParametricPathExtended<VDimension>
	::EvaluateRadius( const InputType & input ) const
	{
		RadiusType radius0;
		RadiusType radius1;
		double     fractionOfLineSegment;
		
		if(  input  <=  this->StartOfInput()  )
    {
			return m_RadiusList.front(); // radius for the first vertex
    }
		else if(  input  >=  this->EndOfInput()  )
    {
			return m_RadiusList.back(); // radius for the last vertex
    }
		
		radius0 = m_RadiusList[(int(input))];
		radius1 = m_RadiusList[((int)input) + 1];
		
		fractionOfLineSegment = input - int(input);
		
		return radius0 + (radius1-radius0)*fractionOfLineSegment;
	}
	
	template<unsigned int VDimension>
	typename PolyLineParametricPathExtended<VDimension>::OffsetType
	PolyLineParametricPathExtended<VDimension>
	::IncrementInput(InputType & input) const
	{
		int         iterationCount;
		bool        tooSmall;
		bool        tooBig;
		InputType   inputStepSize;
		InputType   finalInputValue;
		OffsetType  offset;
		IndexType   currentImageIndex;
		IndexType   nextImageIndex;
//		IndexType   finalImageIndex;
		InputType   nextInputValue;
		InputType		multipFactor;
		InputType		divFactor;
		unsigned long		mainIterationCount;
		
		iterationCount    = 0;
		mainIterationCount = 0;
		inputStepSize     = Superclass::m_DefaultInputStepSize;
		
		// Try to find a point, which is outside the current pixel and
		// smaller than or equal to the final point.
		finalInputValue   = this->EndOfInput();
		currentImageIndex = this->EvaluateToIndex( input );
		nextInputValue = input;
		offset = this->GetZeroOffset();
		while( (offset == this->GetZeroOffset()) && (nextInputValue < finalInputValue))
		{
			
			nextInputValue = floor(nextInputValue) + 1.0;
			if( nextInputValue > finalInputValue )
			{
				nextInputValue = finalInputValue;
			}
			nextImageIndex   = this->EvaluateToIndex( nextInputValue );
			offset            = nextImageIndex - currentImageIndex;		
		}
		
		// If we reach the end of the path, then the whole path is within a single pixel
		// so give up and return a zero offset.
		if( (offset == this->GetZeroOffset()) && !(nextInputValue < finalInputValue) )
		{
			return offset;
		}
		
		// At this point, we are sure that the offset is not equal to zero. 
		// Check if we have already found a solution.
		tooBig = false;
		for( unsigned int i=0; i<VDimension && !tooBig; i++ )
		{
			tooBig = ( offset[i] >= 2 || offset[i] <= -2 );
		}
		if( !tooBig )
		{
			input = nextInputValue;
			return offset;
		}
		
		// If we are not that lucky, then start the numeric search for a pixel that 
		// is in the neighborhood of the current one.
		
//		currentImageIndex = this->EvaluateToIndex( input );
//		finalImageIndex   = this->EvaluateToIndex( finalInputValue );
//		offset            = finalImageIndex - currentImageIndex;
//		if(  ( offset == this->GetZeroOffset() && input != this->StartOfInput() )  ||
//       ( input >=nextInputValue )  )
//    {
//			return this->GetZeroOffset();
//    }
		
		multipFactor = 2;
		divFactor = 1.5;
		do
    {
			if( iterationCount++ > 1000 )		// BUG 1  in the original PolyLineParametricPath!
			{
				multipFactor = 1.0 + (multipFactor - 1.0) / 1.5;
				divFactor = 1.0 + (divFactor - 1.0) / 1.5;
				mainIterationCount++;
				iterationCount = 0;
				
				if( mainIterationCount > 100 )
				{
					// TODO Delete the next line
					std::cout << input << "       " << nextInputValue << "     " 
										<< finalInputValue << "	" 
										<< this->EvaluateToIndex( input ) << "       " << this->EvaluateToIndex( nextInputValue ) << "     " 
										<< this->EvaluateToIndex( finalInputValue ) << "	" 
										<< this->EvaluateToIndex( input + inputStepSize ) << " " 
										<< this->Evaluate( input ) << "       " << this->Evaluate( nextInputValue ) << "     " 
										<< this->Evaluate( finalInputValue ) << "	" 
										<< this->Evaluate( input + inputStepSize ) << " " 					
										<< inputStepSize << "	"
										<< offset << "	"					
										<< multipFactor << " " << divFactor << std::endl;
					
					
					itkExceptionMacro(<<"Too many iterations");
				}
			}
			
			nextImageIndex    = this->EvaluateToIndex( input + inputStepSize );
			offset            = nextImageIndex - currentImageIndex;
			
			tooBig = false;
			tooSmall = ( offset == this->GetZeroOffset() );
			if( tooSmall )
      {
				// double the input step size, but don't go past the end of the input
				inputStepSize *= multipFactor;
				if(  (input + inputStepSize) >= nextInputValue  )
        {
					inputStepSize = nextInputValue - input;
        }
      }
			else
      {
				// Search for an offset dimension that is too big
				for( unsigned int i=0; i<VDimension && !tooBig; i++ )
        {
					tooBig = ( offset[i] >= 2 || offset[i] <= -2 );
        }
				
				if( tooBig )
        {
					inputStepSize /= divFactor;
        }
      }
    }
		while( tooSmall || tooBig );
		
		input += inputStepSize;
		return offset;
	}
	
	
	// TODO You can also create an itk filter (e.g., itk::PathDuplicator) 
	// that extends itk::PathToPathFilter out of this function.
	template <unsigned int VDimension>
	typename PolyLineParametricPathExtended<VDimension>::Pointer
	PolyLineParametricPathExtended<VDimension>::
	Clone() const
	{
		Pointer clone = Self::New();
		
		for(typename VertexListType::ElementIdentifier i = 0; 
				i < m_VertexList->Size(); 
				i++)
		{
			clone->AddVertex( m_VertexList->ElementAt(i), m_RadiusList[i] );
		}
		
		return clone;
	}
	
	/** Copies the given path input to this path. */
	template <unsigned int VDimension>
	void
	PolyLineParametricPathExtended<VDimension>::
	Copy(const Self* path)
	{
		if( path != this)
		{
			this->Initialize();
			
			for(typename VertexListType::ElementIdentifier i = 0; 
					i < path->GetVertexList()->Size(); 
					i++)
			{
				m_VertexList->InsertElement( m_VertexList->Size(), path->Evaluate(i) );
			}
			
			m_RadiusList = path->GetRadiusList();
			
			this->Modified();
		}
	}
	
	
	// TODO You can also create an itk filter 
	// that extends itk::PathToPathFilter out of this function.
	template <unsigned int VDimension>
	void
	PolyLineParametricPathExtended<VDimension>::
	AppendEnd(const Self* path)
	{
		const VertexListType* extension = path->GetVertexList();
		const RadiusListType& extensionRadi = path->GetRadiusList();
		
		// Add the first vertex with the costly AddVertex() function
		// which checks for duplication of the last element. 
		if( extension->Size() > 0 )
		{
			this->AddVertex( extension->ElementAt(0), extensionRadi[0] );
		}
		
		for(typename VertexListType::ElementIdentifier i = 1; 
				i < extension->Size(); 
				i++)
		{
			m_VertexList->InsertElement( m_VertexList->Size(), extension->ElementAt(i) );
			m_RadiusList.push_back( extensionRadi[i] );
		}
		
		this->Modified();
	}
	
	// TODO You can also create an itk filter 
	// that extends itk::PathToPathFilter out of this function.
	template <unsigned int VDimension>
	void
	PolyLineParametricPathExtended<VDimension>::
	Split(InputType splitPoint, Self* path1, Self* path2) const
	{
		if( (path1 == 0) && (path2 == 0) )
		{
			itkExceptionMacro(<<"Both paths are null.");
		}
		
		if( path1 != 0 )
		{
			path1->Initialize();
		}
		
		if( path2 != 0 )
		{
			path2->Initialize();
		}
		
		if( (splitPoint - this->StartOfInput()) < m_Epsilon )
		{
			// Split point is smaller than or equal to StartOfInput()
			// The first path will include only the first vertex on the path.											
			splitPoint = this->StartOfInput();
		}
		if( (this->EndOfInput() - splitPoint) < m_Epsilon )
		{
			// Split point is greater than or equal to EndOfInput()
			// The first path will include only the last vertex on the path.
			splitPoint = this->EndOfInput();
		}
	
		unsigned long splitIndex = static_cast<unsigned long>(splitPoint - this->StartOfInput());
		
		if( path1 != 0 )
		{
			for(unsigned long j = 0; j <= splitIndex; j++)
			{
				path1->AddVertex( m_VertexList->ElementAt(j), m_RadiusList[j] );
			}
			
			if(vcl_fabs(((InputType)splitIndex) - (splitPoint - this->StartOfInput())) >
				 m_Epsilon )
			{
				path1->AddVertex( this->Evaluate(splitPoint), this->EvaluateRadius(splitPoint) );
			}
		}

		if( path2 != 0 )
		{
			// Path1 and path 2 share the middle point.
			path2->AddVertex( this->Evaluate(splitPoint), this->EvaluateRadius(splitPoint) );
			for(unsigned long j = splitIndex+1; j < m_VertexList->Size(); j++)
			{
				path2->AddVertex( m_VertexList->ElementAt(j), m_RadiusList[j] );
			}
		}
	}
	
	
	
	/** Reverses the direction of the path. */
	template <unsigned int VDimension>
	void
	PolyLineParametricPathExtended<VDimension>::
	Reverse()
	{
		if( m_VertexList->Size() != 0 )
		{
			VertexListPointer oldVertexList = m_VertexList;
			
			m_VertexList = VertexListType::New();
			
			typename VertexListType::ElementIdentifier length;
			length = oldVertexList->Size();
			m_VertexList->Reserve( length );
			
			for(typename VertexListType::ElementIdentifier i = 0;
					i < length;
					i++)
			{
				m_VertexList->SetElement( i, oldVertexList->ElementAt(length - i - 1) );
			}
			
			std::reverse(m_RadiusList.begin(), m_RadiusList.end());
			
			this->Modified();
		}
	}
	
	
	template <unsigned int VDimension>
	void
	PolyLineParametricPathExtended<VDimension>::
	SubPath(InputType startPoint, InputType endPoint, Self* path) const
	{
		path->Initialize();
		
		if( startPoint > endPoint )
		{
			// TODO Delete this check and even the whole if statement.
			itkWarningMacro(<<"start point is smaller than the ene point: "<< startPoint<< " " << endPoint);
			
			std::swap(startPoint, endPoint);
		}
		
		
		if( startPoint < this->StartOfInput() )
		{
			startPoint = this->StartOfInput();
		}		
		if( startPoint > this->EndOfInput() )
		{
			startPoint = this->EndOfInput();
		}
		
		if( endPoint < this->StartOfInput() )
		{
			endPoint = this->StartOfInput();
		}		
		if( endPoint > this->EndOfInput() )
		{
			endPoint = this->EndOfInput();
		}
		
		// TODO Delete this check whenever it becomes annoying since it fills the console.
		if( vcl_fabs(startPoint - endPoint) < m_Epsilon )
		{
			itkWarningMacro(<<"Start and end points are the same: " 
											<< startPoint
											<<" The extracted path will include only one vertex.");
		}
		
		unsigned long startIndex = static_cast<unsigned long>(startPoint - this->StartOfInput()) + 1;
		unsigned long endIndex = static_cast<unsigned long>(endPoint - this->StartOfInput());
		if( (static_cast<InputType>(startIndex) - startPoint) < m_Epsilon )
		{
			startIndex++;
		}
		if( (endPoint - static_cast<InputType>(endIndex)) < m_Epsilon )
		{
			if( endIndex != 0 )
			{
				endIndex--;
			}
		}
		
		path->AddVertex( this->Evaluate(startPoint), this->EvaluateRadius(startPoint) );
		for(unsigned long j = startIndex; j <= endIndex; j++)
		{
			path->AddVertex( m_VertexList->ElementAt(j), m_RadiusList[j] );
		}
		path->AddVertex( this->Evaluate(endPoint), this->EvaluateRadius(endPoint) );
	}
	
	
	template <unsigned int VDimension>
	void
	PolyLineParametricPathExtended<VDimension>::
	RoundIndices( bool removeRepeatedVertices )
	{
		if( m_VertexList->empty() )
		{
			return;
		}
		
		for(typename VertexListType::ElementIdentifier i = 0;
					i < m_VertexList->Size();
					i++)
		{
			VertexType vertex = m_VertexList->ElementAt(i);
		
			for(unsigned int j = 0; j < VDimension; j++)
			{
				vertex[j] = Math::Round<typename VertexType::ValueType>(vertex[j]);
			}
			
			m_VertexList->SetElement(i, vertex);
		}
		
		// Traverse the list and remove the repeated vertices.
		if( removeRepeatedVertices )
		{
			VertexListPointer newVertexList = VertexListType::New();
			newVertexList->Initialize();
			VertexType prevVertex = m_VertexList->ElementAt(0);
			newVertexList->InsertElement(newVertexList->Size(), prevVertex);
			for(typename VertexListType::ElementIdentifier i = 1;
						i < m_VertexList->Size();
						i++)
			{
				VertexType vertex = m_VertexList->ElementAt(i);

				// Check if the two vertices are the same.
				bool isEqual = 
				(vertex.SquaredEuclideanDistanceTo( prevVertex ) < 
				 m_Epsilon);
				
				if( !isEqual )
				{
					newVertexList->InsertElement(newVertexList->Size(), vertex);
				}
				
				prevVertex = vertex;
			}
			
			m_VertexList = newVertexList;
		}
	}
	
	
	template <unsigned int VDimension>
	template<class TImage>
	void 
	PolyLineParametricPathExtended<VDimension>::
	Resample(double stepInWorldCoords, const TImage* image)
	{
		if( m_VertexList->Size() < 2 )
		{
			return;
		}
	
		Pointer clone = this->Clone();
		this->Initialize(); // clears the vertex and radius lists.
		
		// Add the first vertex.
		this->AddVertex(clone->GetVertexList()->ElementAt(0), 
										clone->GetRadiusList()[0]);
		
		// Resample the rest.
		bool endReached = false;
		InputType currentPos = clone->StartOfInput();
		InputType nextPos;
		double traversedDistance; // dummy
		while( !endReached )
		{
			endReached = 
			!(clone->TraverseDistanceInWorldCoords(currentPos,
																						 image,
																						 stepInWorldCoords,
																						 true,
																						 nextPos,
																						 traversedDistance));
			if( !endReached )
			{
				this->AddVertex(clone->Evaluate(nextPos),
												clone->EvaluateRadius(nextPos));
			}
			
			currentPos = nextPos;
		}
		
		// Add the last vertex.
		this->AddVertex(clone->GetVertexList()->ElementAt(clone->GetVertexList()->Size()-1), 
										clone->GetRadiusList().back());
	}
	
	/** Smooths the vertex locations and the radius values. */
	template <unsigned int VDimension>
	template<class TImage>
	void 
	PolyLineParametricPathExtended<VDimension>::	
	SmoothVertexLocationsAndRadii(double avgWindowRadiusInWorldCoords,
																const TImage* image)
	{
		Smooth(avgWindowRadiusInWorldCoords, image, true, true);
	}
	
	/** Smooths the vertex locations. */
	template <unsigned int VDimension>
	template<class TImage>
	void 
	PolyLineParametricPathExtended<VDimension>::	
	SmoothVertexLocations(double avgWindowRadiusInWorldCoords,
												const TImage* image)
	{
		Smooth(avgWindowRadiusInWorldCoords, image, true, false);
	}
	
	/** Smooths the radius values along the path. */
	template <unsigned int VDimension>
	template<class TImage>
	void 
	PolyLineParametricPathExtended<VDimension>::	
	SmoothRadiusValues(double avgWindowRadiusInWorldCoords,
										 const TImage* image)
	{
		Smooth(avgWindowRadiusInWorldCoords, image, false, true);
	}
	
	/** Smooths the vertex locations. */
	template <unsigned int VDimension>
	template<class TImage>
	void 
	PolyLineParametricPathExtended<VDimension>::	
	Smooth(double avgWindowRadiusInWorldCoords,
				 const TImage* image,
				 bool smoothLocations,
				 bool smoothRadii)
	{
		if( (!smoothLocations) && (!smoothRadii) )	
		{
			return;
		}
		
		Pointer clone = this->Clone();
		
		// Visit and smooth all the vertices except the first and the last one.
		for(typename VertexListType::ElementIdentifier i = 1; 
				(i+1) < m_VertexList->Size(); 
				i++)
		{
			double backDistance;
			InputType backPoint;
			double forwardDistance;
			InputType forwardPoint;
			
			clone->TraverseDistanceInWorldCoords(i,
																					 image,
																					 avgWindowRadiusInWorldCoords,
																					 false,
																					 backPoint,
																					 backDistance);
			
			clone->TraverseDistanceInWorldCoords(i,
																					 image,
																					 avgWindowRadiusInWorldCoords,
																					 true,
																					 forwardPoint,
																					 forwardDistance);
			
			if(vcl_fabs(backDistance - forwardDistance) > m_Epsilon)
			{
				if(backDistance < forwardDistance)
				{
					clone->TraverseDistanceInWorldCoords(i,
																							 image,
																							 backDistance,
																							 true,
																							 forwardPoint,
																							 forwardDistance);
				}
				else
				{
					clone->TraverseDistanceInWorldCoords(i,
																							 image,
																							 forwardDistance,
																							 false,
																							 backPoint,
																							 backDistance);
				}
			}
			
			// Extract the segment for the moving average window.
			Pointer segment = Self::New();
			clone->SubPath(backPoint, forwardPoint, segment);
			
			if( smoothLocations )
			{
				// Compute the centroid location of the segment and set it as 
				// the new path vertex.
				m_VertexList->SetElement(i, segment->ComputeCentroidVertex(image));
			}
			
			if( smoothRadii )
			{
				// Compute the centroid radius of the segment and set it as 
				// the new radius.
				m_RadiusList[i] = segment->ComputeCentroidRadius(image);
			}
		}
	}


	template <unsigned int VDimension>
	typename PolyLineParametricPathExtended<VDimension>::ImageRegionType 
	PolyLineParametricPathExtended<VDimension>::	
	ComputeBoundingImageRegionWithoutRadius() const
	{
		ImageRegionType imageRegion;
		
		ComputeBoundingImageRegionFromVertexList(m_VertexList.GetPointer(),
																						 imageRegion);
																						 
		return imageRegion;
	}

	template <unsigned int VDimension>
	typename PolyLineParametricPathExtended<VDimension>::ImageRegionType 
	PolyLineParametricPathExtended<VDimension>::	
	ComputeBoundingImageRegionWithRadius(const SpacingType& spacing) const
	{
		ImageRegionType imageRegion;
			
		VertexListPointer vertexListWithRadiusOffsets = VertexListType::New();
																							 
		for(typename VertexListType::ElementIdentifier i = 0; 
				i < m_VertexList->Size(); 
				i++)
		{
			RadiusType radius = m_RadiusList[i];
		
			// Compute the offset vector for the given radius value.
			RadiusType offset;
			VertexType vertex = m_VertexList->ElementAt(i);
			VertexType vertexTopLeft;
			VertexType vertexBottomRight;
			for(unsigned int j = 0; j < Dimension; j++)
			{
				offset = radius / spacing[j];
				
				vertexTopLeft[j] = vertex[j] + offset;
				vertexBottomRight[j] = vertex[j] - offset;
			}
		
			// Insert two points in the list one with the offset 
			// added and the other subtracted.
			vertexListWithRadiusOffsets->InsertElement(
			vertexListWithRadiusOffsets->Size(), vertexTopLeft);
			vertexListWithRadiusOffsets->InsertElement(
			vertexListWithRadiusOffsets->Size(), vertexBottomRight);
		}
		
		ComputeBoundingImageRegionFromVertexList(vertexListWithRadiusOffsets.GetPointer(),
																						 imageRegion);
																						 
		return imageRegion;
	}
	
	template <unsigned int VDimension>
	typename PolyLineParametricPathExtended<VDimension>::VertexType 
	PolyLineParametricPathExtended<VDimension>::	
	ComputeMeanVertex() const
	{
		VertexType meanVertex;
		meanVertex.Fill( NumericTraits<typename VertexType::ValueType>::ZeroValue() );
		if( m_VertexList->Size() == 0 )
		{
			return meanVertex;
		}
		
		for(typename VertexListType::ElementIdentifier i = 0; 
				i < m_VertexList->Size(); 
				i++)
		{
			VertexType vertex = m_VertexList->ElementAt(i);
			for(unsigned int j = 0; j < Dimension; j++)
			{
				meanVertex[j] += vertex[j];
			}
		}
		
		for(unsigned int j = 0; j < Dimension; j++)
		{
			meanVertex[j] /= static_cast<typename VertexType::ValueType>(m_VertexList->Size());
		}
		
		return meanVertex;
	}
	
	template <unsigned int VDimension>
	template<class TImage>
	typename PolyLineParametricPathExtended<VDimension>::VertexType 
	PolyLineParametricPathExtended<VDimension>::
	ComputeCentroidVertex(TImage* image) const
	{
		VertexType meanVertex;
		meanVertex.Fill( NumericTraits<typename VertexType::ValueType>::ZeroValue() );
		if( m_VertexList->Size() == 0 )
		{
			return meanVertex;
		}
	
		if( m_VertexList->Size() == 1 )
		{
			return m_VertexList->ElementAt(0);
		}
		
		double segmentLength;
		double pathLength = 0.0;
		for(typename VertexListType::ElementIdentifier i = 1; 
				i < m_VertexList->Size(); 
				i++)
		{
			VertexType midVertex = Evaluate(static_cast<InputType>(i) - 0.5);
			segmentLength = GetDistanceInWorldCoords(i-1, i, image);
			pathLength += segmentLength;
			for(unsigned int j = 0; j < Dimension; j++)
			{
				meanVertex[j] += (midVertex[j] * segmentLength);
			}
		}
		
		for(unsigned int j = 0; j < Dimension; j++)
		{
			meanVertex[j] /= static_cast<typename VertexType::ValueType>(pathLength);
		}
		
		return meanVertex;
	}
	
	
	
	template <unsigned int VDimension>
	typename PolyLineParametricPathExtended<VDimension>::RadiusType 
	PolyLineParametricPathExtended<VDimension>::	
	ComputeMinRadius() const
	{
		if( m_RadiusList.empty() )
		{
			return NumericTraits<RadiusType>::ZeroValue();
		}
	
		RadiusType minRadius = m_RadiusList[0];
		for(typename RadiusListType::size_type i = 1; 
				i < m_RadiusList.size(); 
				i++)
		{
			if( minRadius > m_RadiusList[i] )
			{
				minRadius = m_RadiusList[i];
			}
		}
		
		return minRadius;
	}
	
	template <unsigned int VDimension>
	typename PolyLineParametricPathExtended<VDimension>::RadiusType 
	PolyLineParametricPathExtended<VDimension>::	
	ComputeMaxRadius() const
	{
		if( m_RadiusList.empty() )
		{
			return NumericTraits<RadiusType>::ZeroValue();
		}
	
		RadiusType maxRadius = m_RadiusList[0];
		for(typename RadiusListType::size_type i = 1; 
				i < m_RadiusList.size(); 
				i++)
		{
			if( maxRadius < m_RadiusList[i] )
			{
				maxRadius = m_RadiusList[i];
			}
		}
		
		return maxRadius;
	}

	template <unsigned int VDimension>
	typename PolyLineParametricPathExtended<VDimension>::RadiusType 
	PolyLineParametricPathExtended<VDimension>::	
	ComputeMeanRadius() const
	{
		RadiusType meanRadius = NumericTraits<RadiusType>::ZeroValue();
		if( m_RadiusList.empty() )
		{
			return meanRadius;
		}
	
		for(typename RadiusListType::size_type i = 0; 
				i < m_RadiusList.size(); 
				i++)
		{
			meanRadius += m_RadiusList[i];
		}
		
		meanRadius /= static_cast<RadiusType>(m_RadiusList.size());
		
		return meanRadius;
	}
	
	template <unsigned int VDimension>
	template<class TImage>
	typename PolyLineParametricPathExtended<VDimension>::RadiusType 
	PolyLineParametricPathExtended<VDimension>::
	ComputeCentroidRadius(TImage* image) const
	{
		RadiusType meanRadius = NumericTraits<RadiusType>::ZeroValue();
		if( m_VertexList->Size() == 0 )
		{
			return meanRadius;
		}
		
		if( m_VertexList->Size() == 1 )
		{
			return m_RadiusList[0];
		}
		
		double segmentLength;
		double pathLength = 0.0;
		for(typename VertexListType::ElementIdentifier i = 1; 
				i < m_VertexList->Size(); 
				i++)
		{
			RadiusType midRadius = EvaluateRadius(static_cast<InputType>(i) - 0.5);
			segmentLength = GetDistanceInWorldCoords(i-1, i, image);
			pathLength += segmentLength;
			meanRadius += (midRadius * segmentLength);
		}
		
		meanRadius /= static_cast<RadiusType>(pathLength);
		
		return meanRadius;
	}
	
	
	template <unsigned int VDimension>
	void 
	PolyLineParametricPathExtended<VDimension>::	
	ScaleAndShiftRadii(double scaleFactor,
										 double shift)
	{
		for(typename RadiusListType::iterator it = m_RadiusList.begin(); 
				it != m_RadiusList.end(); 
				it++)
		{
			(*it) = static_cast<RadiusType>(
			(static_cast<double>(*it) * scaleFactor) + shift);
		}
	}
		
	template <unsigned int VDimension>
	void 
	PolyLineParametricPathExtended<VDimension>::
	ScaleAndShiftVertices(double scaleFactor,
												double shift)
	{
		for(typename VertexListType::ElementIdentifier i = 0; 
				i < m_VertexList->Size(); 
				i++)
		{
			VertexType& vertex = m_VertexList->ElementAt(i);
			for(unsigned int j = 0; j < Dimension; j++)
			{
				vertex[j] = (vertex[j] * scaleFactor) + shift;
			}
		}
	}
	
	
	template <unsigned int VDimension>
	template<class TImage>
	void 
	PolyLineParametricPathExtended<VDimension>::	
	GetUniqueImageIndices(const TImage* image, 
												std::vector<IndexType>& uniqueIndices) const
	{
		typedef typename OffsetType::OffsetValueType	OffsetValueType;
	
		// Get the list of image pixel offsets iterated by the path iterator.
		std::vector<OffsetValueType> offsets;
		PathConstIterator<TImage, Self> iter(image, this);
		for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
		{
			offsets.push_back( image->ComputeOffset( iter.GetIndex() ) );
		}
		
		// Remove repeated elements from the list (complexity nlog(n)).
		typename std::vector< OffsetValueType >::iterator arrayIter1;
		typename std::vector< OffsetValueType >::iterator arrayIter2;
		std::set< OffsetValueType > tmpIndexSet;
		for(arrayIter1 = offsets.begin() , arrayIter2 = offsets.begin() ; 
				arrayIter1 != offsets.end() ; ++arrayIter1 )
		{
			if( tmpIndexSet.insert( *arrayIter1 ).second )
			{
				*arrayIter2++ = *arrayIter1;
			}
		}
		offsets.erase( arrayIter2 , offsets.end() );
		
		// Convert the unique offset list to an index list.
		uniqueIndices.resize( offsets.size() );
		for(typename std::vector< OffsetValueType >::size_type i = 0; 
				i < offsets.size();
				i++ )
		{
			uniqueIndices[i] = image->ComputeIndex( offsets[i] );
		}
	}
	
	template <unsigned int VDimension>	
	template<class TImage>
	typename PolyLineParametricPathExtended<VDimension>::Pointer 
	PolyLineParametricPathExtended<VDimension>::	
	CreatePathFromUniqueImageIndices(const TImage* image) const
	{
		Pointer path = Self::New();
		
		// Get the unique index list.
		std::vector<IndexType> uniqueIndices;
		this->GetUniqueImageIndices(image, uniqueIndices);
		
		// Check if the last unique index is the end of this path.
		// If not, find the location of the end index in the unique index list, 
		// and remove the rest of the unique indices. This essentially removes 
		// the loop that includes the end vertex.
		IndexType endIndex;
		endIndex.CopyWithRound( m_VertexList->ElementAt(m_VertexList->Size() - 1) );
		if( uniqueIndices.back() != endIndex )
		{
			typename std::vector<IndexType>::reverse_iterator revIndexIter;
			for(revIndexIter = uniqueIndices.rbegin(); 
					revIndexIter < uniqueIndices.rend();
					revIndexIter++ )
			{
				if( (*revIndexIter) == endIndex )
				{
					uniqueIndices.erase(uniqueIndices.end() - (revIndexIter - uniqueIndices.rbegin()), 
					uniqueIndices.end());
					break;
				}
			}
		}
		
		// Find the elements of the vertex list and the radius list
		// that correspond to the unique indices.
		typename std::vector<IndexType>::iterator indexIter;
		for(indexIter = uniqueIndices.begin(); 
				indexIter != uniqueIndices.end(); 
				indexIter++ )
		{
			// We don't evaluate the continuous vertex index
			// and add it (by path->AddVertex( this->Evaluate(pos) ...), which 
			// may be different than the index. That is why, for the radius values 
			// we estimate the path position that is closest to the current 
			// discrete index and then evaluate radius at that point of the path. 
			// We get the path position corresponding to the closest point 
			// by callng the GetEuclDistToPointInImageCoords() method.
			InputType pos = this->GetEuclDistToPointInImageCoords( *indexIter).second;
			path->AddVertex( *indexIter, this->EvaluateRadius(pos) );
		}
		
		return path;
	}
	
	template <unsigned int VDimension>
	std::pair<double,  typename PolyLineParametricPathExtended<VDimension>::InputType>
	PolyLineParametricPathExtended<VDimension>::	
	GetEuclDistToPointInImageCoords(const VertexType& vertex) const
	{
		InputType minDistPoint = 0.0;
		
		const VertexListType* vertexList = this->GetVertexList();

		double minDist = NumericTraits<double>::max();
		
		// Check if there is a only a single element in the list.
		// Note that if there is no element in the list, then the function 
		// will return max value for double.
		if( this->GetVertexList()->Size() == 1 )
		{
			minDist = vertexList->ElementAt(0).EuclideanDistanceTo(vertex);
		}
		else if( this->GetVertexList()->Size() > 1 )
		{
			// Traverse the list of edges and find the shorthest distance for each of them.
			VertexType prevVertex = vertexList->ElementAt(0);
			for(typename VertexListType::ElementIdentifier i = 1;
				i < vertexList->Size(); 
				i++)
			{
				VertexType currentVertex = vertexList->ElementAt(i);
				InputType minDistPointForLine;
				double minDistForLine = 
				this->GetEuclDistToLineSegment(prevVertex, currentVertex, 
																				vertex, minDistPointForLine);
			
				if( minDist > minDistForLine )
				{
					minDist = minDistForLine;
					minDistPoint = static_cast<InputType>(i-1) + minDistPointForLine;
				}
				
				prevVertex = currentVertex;
			}
		}
		
		return std::make_pair(minDist, minDistPoint);
	}
	
	template <unsigned int VDimension>
	template <class TImage>
	std::pair<double,  typename PolyLineParametricPathExtended<VDimension>::InputType>
	PolyLineParametricPathExtended<VDimension>::	
	GetEuclDistToPointInWorldCoords(const VertexType& vertex, 
																	const TImage* image) const
	{
		PointType point;
		PointType prevPoint;
		PointType currentPoint;
		InputType minDistPoint = 0.0;
		
		image->TransformContinuousIndexToPhysicalPoint(vertex, point);
		
		double minDist = NumericTraits<double>::max();
		
		// Check if there is a only a single element in the list.
		// Note that if there is no element in the list, then the function 
		// will return max value for double.
		if( m_VertexList->Size() == 1 )
		{
			image->TransformContinuousIndexToPhysicalPoint(m_VertexList->ElementAt(0), prevPoint);
			minDist = prevPoint.EuclideanDistanceTo(point);
		}
		else if( m_VertexList->Size() > 1 )
		{
			// Traverse the list of edges and find the shorthest distance for each of them.
			image->TransformContinuousIndexToPhysicalPoint(m_VertexList->ElementAt(0), prevPoint);
			for(typename VertexListType::ElementIdentifier i = 1;
					i < m_VertexList->Size(); 
					i++)
			{
				image->TransformContinuousIndexToPhysicalPoint(m_VertexList->ElementAt(i), currentPoint);
				
				InputType minDistPointForLine;
				double minDistForLine = 
				this->GetEuclDistToLineSegment(prevPoint, currentPoint, 
																			 point, minDistPointForLine);
				
				if( minDist > minDistForLine )
				{
					minDist = minDistForLine;
					minDistPoint = static_cast<InputType>(i-1) + minDistPointForLine;
				}
				
				prevPoint = currentPoint;
			}
		}
		
		return std::make_pair(minDist, minDistPoint);
	}
	
	template <unsigned int VDimension>
	template <class TImage>
	bool
	PolyLineParametricPathExtended<VDimension>::	
	IsPointInTailRegion(const PointType& point, 
											const InputType& closestPathInput,
											const TImage* image) const
	{
		// Check if this is an end point and if so, compute the angle between
		// the path direction vector and the vector that is formed by the path point 
		// and the given one.
		bool isInTailRegion = false;
		bool isEndVertex = false;
		VertexType pathEndVertex;
		if(vcl_fabs(closestPathInput - this->StartOfInput()) < m_Epsilon)
		{
			isEndVertex = true;
			
			pathEndVertex = m_VertexList->ElementAt(0); // start vertex		
		}
		else if(vcl_fabs(closestPathInput - this->EndOfInput()) < m_Epsilon)
		{
			isEndVertex = true;
			
			pathEndVertex = 
			m_VertexList->ElementAt(m_VertexList->Size() - 1);	// end vertex
		}
		
		if( isEndVertex )
		{
			// Check if the given point is on the path's crossection 
			// plane that contains the end point.
			PointType endPoint;
			image->TransformContinuousIndexToPhysicalPoint(pathEndVertex, endPoint);
			
			VectorType pathEndPointToGivenPointVector;
			pathEndPointToGivenPointVector = point - endPoint;
			
			VectorType pathDirectionVector = 
			this->EvaluateDirectionInWorldCoords(closestPathInput, image);
			
			// Check if the vectors are orthogonal to each other, if not
			// then the given vertex is in the tail region of the path
			// and not in the path.
			if(vcl_fabs(pathDirectionVector * pathEndPointToGivenPointVector) > 
				 m_Epsilon)
			{
				isInTailRegion = true;
			}
		}
		
		return isInTailRegion;
	}
	
	template <unsigned int VDimension>
	template <class TImage>
	void
	PolyLineParametricPathExtended<VDimension>::
	ExtractClosestPathPointsThatContainsPointInWorldCoords(const VertexType& vertex, 
																												 const TImage* image,
																												 std::pair< double, std::vector<InputType> >& outputPairs,
																												 bool excludeTails) const
	{
		PointType point;
		PointType prevPoint;
		PointType currentPoint;
		double minDist = NumericTraits<double>::max();
		std::vector<InputType> closestPoints;
		outputPairs = std::make_pair(minDist, closestPoints);
		
		// Check if there is a only a single element in the list.
		// Note that if there is no element in the list, then the function 
		// will return max value of double as the distance and an empty set
		// as the set of closest points.
		image->TransformContinuousIndexToPhysicalPoint(vertex, point);
		if( m_VertexList->Size() == 1 )
		{
			image->TransformContinuousIndexToPhysicalPoint(m_VertexList->ElementAt(0), 
																										 prevPoint);
			minDist = prevPoint.EuclideanDistanceTo(point);
			
			if( (minDist - static_cast<double>(m_RadiusList[0])) < m_Epsilon )
			{
				outputPairs.first = minDist;
				outputPairs.second.push_back( 0.0 );
			}
		}
		else if( m_VertexList->Size() > 1 )
		{
			InputType minDistPointForLine;
			double minDistForLine;
			InputType minDistPoint = NumericTraits<InputType>::max();
			
			// Traverse the list of edges and find the shorthest distance for each of them.
			image->TransformContinuousIndexToPhysicalPoint(m_VertexList->ElementAt(0), prevPoint);
			for(typename VertexListType::ElementIdentifier i = 1;
					i < m_VertexList->Size(); 
					i++)
			{
				image->TransformContinuousIndexToPhysicalPoint(m_VertexList->ElementAt(i), currentPoint);
				
				minDistForLine = this->GetEuclDistToLineSegment(prevPoint, currentPoint, 
																												point, minDistPointForLine);
				// Check if this is the same path point as the one returned from the previous
				// line segment.
				if(vcl_fabs(minDistPoint - (static_cast<InputType>(i-1) + minDistPointForLine)) > 
					 m_Epsilon)
				{
					minDistPoint = static_cast<InputType>(i-1) + minDistPointForLine;
					
					// Check if the path point contains the given one.
					if((minDistForLine - static_cast<double>(this->EvaluateRadius(minDistPoint))) < 
						 m_Epsilon )
					{
						// Check if the given point is in the tail region.
						bool isTailRegion = false;
						if( excludeTails )
						{
							isTailRegion = this->IsPointInTailRegion(point, minDistPoint, image);
						}
						
						if( !isTailRegion )
						{					
							// Check if this is one of the closest points so far.
							if( vcl_fabs(minDist - minDistForLine) < m_Epsilon )
							{
								outputPairs.second.push_back( minDistPoint );
							}
							else if( minDist > minDistForLine ) // check if this is the closest point.
							{
								minDist = minDistForLine;
								
								outputPairs.first = minDist;
								outputPairs.second.clear();
								outputPairs.second.push_back( minDistPoint );
							}
						}
					}
				}
				
				prevPoint = currentPoint;
			}
		}
	}

	
	template <unsigned int VDimension>
	template <class TImage>
	bool  
	PolyLineParametricPathExtended<VDimension>::
	DoesContainPoint(const VertexType& vertex,
									 const TImage* image,
									 bool excludeTails) const
	{
		PointType point;
		PointType prevPoint;
		PointType currentPoint;
	
		// Check if there is a only a single element in the list.
		// Note that if there is no element in the list, then the function 
		// will return false.
		bool doesContain = false;
		image->TransformContinuousIndexToPhysicalPoint(vertex, point);
		if( m_VertexList->Size() == 1 )
		{
			image->TransformContinuousIndexToPhysicalPoint(m_VertexList->ElementAt(0), prevPoint);
			doesContain = (static_cast<double>(prevPoint.EuclideanDistanceTo(point)) - 
			static_cast<double>(this->GetRadiusList()[0])) < m_Epsilon;
		}
		else if( m_VertexList->Size() > 1 )
		{
			// Traverse the list of edges and find the shorthest distance for each of them.
			image->TransformContinuousIndexToPhysicalPoint(m_VertexList->ElementAt(0), prevPoint);
			for(typename VertexListType::ElementIdentifier i = 1;
					(i < m_VertexList->Size()) && (!doesContain); 
					i++)
			{
				image->TransformContinuousIndexToPhysicalPoint(m_VertexList->ElementAt(i), currentPoint);
				
				InputType minDistPoint;
				double minDistForLine = 
				this->GetEuclDistToLineSegment(prevPoint, currentPoint, 
																			 point, minDistPoint);
				minDistPoint = static_cast<InputType>(i-1) + minDistPoint;
				
				// Check if the given point is inside the line segment.
				if((minDistForLine - static_cast<double>(this->EvaluateRadius(minDistPoint))) < 
					 m_Epsilon )
				{
					doesContain = true;
					
					if( excludeTails )
					{
						doesContain = !(this->IsPointInTailRegion(point, minDistPoint, image));
					}
				}
				
				prevPoint = currentPoint;
			}
		}
	 
		return doesContain;
	}
	
	template <unsigned int VDimension>
	template <class TImage>
	bool  
	PolyLineParametricPathExtended<VDimension>::
	DoesContainPointsOfPath(const Self* path,
													const TImage* image,
													bool excludeTails) const
	{
		bool containVertices = true;
	
		// Check if this path contains all the vertices of 
		// the given path. 
		const VertexListType* vertexList = path->GetVertexList();
		typename VertexListType::ElementIdentifier noOfVertices = 
		vertexList->Size();
		for(typename VertexListType::ElementIdentifier i = 0;
				(i < noOfVertices) && containVertices; 
				i++)
		{
			containVertices = DoesContainPoint(vertexList->ElementAt(i), 
																				 image, excludeTails);
		}
		
		return containVertices;
	}
	
	
	template <unsigned int VDimension>
	template <class TImage>
	double  
	PolyLineParametricPathExtended<VDimension>::
	ComputeMaxPathSegmLengthOutsidePath(const Self* path,
																			const TImage* image,
																			bool excludeTails) const
	{
		// Compute the minimum image spacing.
		const typename TImage::SpacingType& spacing = image->GetSpacing();
		double minSpacing = spacing.GetVnlVector().min_value();
		
		double length = 0.0;
		double max_length = 0.0;
		const VertexListType* vertexList = this->GetVertexList();
		typename VertexListType::ElementIdentifier noOfVertices = 
		vertexList->Size();
		PointType currPoint;
		PointType prevPoint;
		bool prevFlag = false;
		bool currFlag = false;
		image->TransformContinuousIndexToPhysicalPoint
		(vertexList->ElementAt(0), prevPoint);
		for(typename VertexListType::ElementIdentifier i = 1;
				i < noOfVertices; 
				i++)
		{
			// Determine the step size to take.
			image->TransformContinuousIndexToPhysicalPoint
			(vertexList->ElementAt(i), currPoint);
			
			double dist = currPoint.EuclideanDistanceTo(prevPoint);
			unsigned long numberOfSteps = 
			itk::Math::Ceil<unsigned long>(dist / minSpacing);
			double worldStep = dist / static_cast<double>(numberOfSteps);
			InputType inputStep = 1.0 / static_cast<InputType>(numberOfSteps);
			
			InputType currInput = 
			this->StartOfInput() + static_cast<InputType>(i-1) + 
			inputStep;	// discard the first point.
			for(unsigned long j = 0;
					j < numberOfSteps; 
					j++)
			{
				currFlag = !(path->DoesContainPoint
										 (this->Evaluate(currInput), image, excludeTails));
				
				if( currFlag )
				{
					length += worldStep;
				}
				else if( (!currFlag) && prevFlag )
				{
					if( max_length < length )
					{
						max_length = length;
					}
					length = 0.0;		// reset the accumulated length.
				}
				
				currInput += inputStep;
				prevFlag = currFlag;
			}
			if( max_length < length )
			{
				max_length = length;
			}			
			
			prevPoint = currPoint;
		}
		
		return max_length;
	}
	
	template <unsigned int VDimension>
	template <class TImage>
	double  
	PolyLineParametricPathExtended<VDimension>::
	ComputePathLengthInsidePath(const Self* path,
														  const TImage* image,
															bool excludeTails) const
	{
		// Compute the minimum image spacing.
		const typename TImage::SpacingType& spacing = image->GetSpacing();
		double minSpacing = spacing.GetVnlVector().min_value();
	
		double length = 0.0;
		const VertexListType* vertexList = this->GetVertexList();
		typename VertexListType::ElementIdentifier noOfVertices = 
		vertexList->Size();
		PointType currPoint;
		PointType prevPoint;
		image->TransformContinuousIndexToPhysicalPoint
		(vertexList->ElementAt(0), prevPoint);
		for(typename VertexListType::ElementIdentifier i = 1;
				i < noOfVertices; 
				i++)
		{
			// Determine the step size to take.
			image->TransformContinuousIndexToPhysicalPoint
			(vertexList->ElementAt(i), currPoint);
			
			double dist = currPoint.EuclideanDistanceTo(prevPoint);
			unsigned long numberOfSteps = 
			itk::Math::Ceil<unsigned long>(dist / minSpacing);
			double worldStep = dist / static_cast<double>(numberOfSteps);
			InputType inputStep = 1.0 / static_cast<InputType>(numberOfSteps);
			
			InputType currInput = 
			this->StartOfInput() + static_cast<InputType>(i-1) + 
			inputStep;	// discard the first point.
			for(unsigned long j = 0;
					j < numberOfSteps; 
					j++)
			{
				if(path->DoesContainPoint
					 (this->Evaluate(currInput), image, excludeTails))
				{
					length += worldStep;
				}
				currInput += inputStep;
			}
			
			prevPoint = currPoint;
		}
		
		return length;
	}
	
	template <unsigned int VDimension>
	template <class TImage>
	double  
	PolyLineParametricPathExtended<VDimension>::
	ComputePathLengthOutsidePath(const Self* path,
															 const TImage* image,
															 bool excludeTails) const
	{
		double length = 
		this->GetLengthInWorldCoords(image) - 
		this->ComputePathLengthInsidePath(path, image, excludeTails);
		
		if( length < 0.0 )
		{
			return 0.0;
		}
		
		return length;
	}
	
	template <unsigned int VDimension>
	template <class TImage>
	double  
	PolyLineParametricPathExtended<VDimension>::
	ComputePathLengthInsideMask(const TImage* maskImage) const
	{
		// Compute the minimum image spacing.
		const typename TImage::SpacingType& spacing = maskImage->GetSpacing();
		double minSpacing = spacing.GetVnlVector().min_value();
	
		double length = 0.0;
		const VertexListType* vertexList = this->GetVertexList();
		typename VertexListType::ElementIdentifier noOfVertices = 
		vertexList->Size();
		PointType currPoint;
		PointType prevPoint;
		maskImage->TransformContinuousIndexToPhysicalPoint
		(vertexList->ElementAt(0), prevPoint);
		for(typename VertexListType::ElementIdentifier i = 1;
				i < noOfVertices; 
				i++)
		{
			// Determine the step size to take.
			maskImage->TransformContinuousIndexToPhysicalPoint
			(vertexList->ElementAt(i), currPoint);
		
			double dist = currPoint.EuclideanDistanceTo(prevPoint);
			unsigned long numberOfSteps = 
			itk::Math::Ceil<unsigned long>(dist / minSpacing);
			double worldStep = dist / static_cast<double>(numberOfSteps);
			InputType inputStep = 1.0 / static_cast<InputType>(numberOfSteps);
			
			InputType currInput = 
			this->StartOfInput() + static_cast<InputType>(i-1) + 
			inputStep;	// discard the first point.
			IndexType currIndex;
			for(unsigned long j = 0;
					j < numberOfSteps; 
					j++)
			{
				currIndex.CopyWithRound( this->Evaluate(currInput) );
				if( static_cast<bool>(maskImage->GetPixel(currIndex)) )
				{
					length += worldStep;
				}
				currInput += inputStep;
			}
			
			prevPoint = currPoint;
		}
		
		return length;
	}
	
	
	template <unsigned int VDimension>
	template <class TImage>
	double  
	PolyLineParametricPathExtended<VDimension>::
	ComputePathLengthOutsideMask(const TImage* maskImage) const
	{
		double length =  
		this->GetLengthInWorldCoords(maskImage) - 
		this->ComputePathLengthInsideMask(maskImage);
		
		if( length < 0.0 )
		{
			return 0.0;
		}
	}
	
	template <unsigned int VDimension>
	template <class TPoint>
	double  
	PolyLineParametricPathExtended<VDimension>::	
	GetEuclDistToLineSegment(const TPoint& lineVertex1, 
													 const TPoint& lineVertex2, 
													 const TPoint& vertex, 
													 InputType& minDistPoint)	const// takes values between 0 (start vertex) and 1 (end vertex).
	{
		double distAB;
		double distBC;
		double distBX;		
		double dot1;
		double dot2;
		double finalDist;
		
		// Check if the line points are very close to each other.
		distAB = lineVertex1.EuclideanDistanceTo(lineVertex2);
		if( distAB < m_Epsilon )
		{
			minDistPoint = 0.0;
			return lineVertex1.EuclideanDistanceTo(vertex);
		}
		
		// Supose that lineVertex1 := A, lineVertex2 := B and vertex := C.
		// Find the dot product of BA and AC.
		dot1 = (lineVertex1 - lineVertex2) * (vertex - lineVertex1) ;
		if( dot1 > -m_Epsilon )
		{
			minDistPoint = 0.0;
			return lineVertex1.EuclideanDistanceTo(vertex);
		}
		
		// Find the dot product of AB and BC.
		dot2 = (lineVertex2 - lineVertex1) * (vertex - lineVertex2) ;
		if( dot2 > -m_Epsilon )
		{
			minDistPoint = 1.0;
			return lineVertex2.EuclideanDistanceTo(vertex);
		}
		
		// At this point, we are sure that the minimum distance 
		// point is inside the line segment.
		distBC = lineVertex2.EuclideanDistanceTo(vertex);
		distBX = - dot2 / distAB;
		minDistPoint = (distAB - distBX) / distAB;
		if( minDistPoint > 1.0 )
		{
			minDistPoint = 1.0;
		}
		if( minDistPoint < 0.0 )
		{
			minDistPoint = 0.0;
		}
		finalDist = distBC*distBC - distBX*distBX;
		if( finalDist < m_Epsilon )
		{
			return 0;
		}
		else
		{
			return vcl_sqrt(finalDist);
		}
	}
	
	
	
	template <unsigned int VDimension>
	template<class TImage>
	void 
	PolyLineParametricPathExtended<VDimension>::
	OverlayFromNeighIter(typename TImage::PixelType overlayVal,
											 ShapedNeighborhoodIterator<TImage>& iter,
											 const Image<bool, TImage::ImageDimension>* binaryMap) const
	{
		typedef TImage																						ImageType;
		typedef ShapedNeighborhoodIterator<ImageType>							NeighborhoodIteratorType;

		ImageRegionType region = iter.GetRegion();
		
		if( binaryMap == 0 )
		{
			for(typename NeighborhoodIteratorType::Iterator neighIter = iter.Begin(); 
					!neighIter.IsAtEnd(); 
					++neighIter)
			{		
				typename ImageType::IndexType imageIndex = 
				iter.GetIndex(neighIter.GetNeighborhoodIndex());
				
				if( region.IsInside( imageIndex ) )
				{
					iter.SetPixel(neighIter.GetNeighborhoodIndex(), overlayVal);
				}
			}
		}
		else
		{
			// The image of the input iterator can be different than the 
			// binaryMap image. That is why, we need to use their overlapping 
			// region.
			bool isInside = region.Crop(binaryMap->GetBufferedRegion());
			
			// Check if the region is outside the buffered 
			// region of the image.
			if( !isInside )
			{
				itkExceptionMacro(<<"The iterator region "
													<<region
													<<" is outside the  "
													<<"BufferedRegion "
													<< binaryMap->GetBufferedRegion()
													<<" of the input binary image.");
			}
			
			for(typename NeighborhoodIteratorType::Iterator neighIter = iter.Begin(); 
					!neighIter.IsAtEnd(); 
					++neighIter)
			{		
				typename ImageType::IndexType imageIndex = 
				iter.GetIndex(neighIter.GetNeighborhoodIndex());
				
				if(region.IsInside( imageIndex ) )
				{
					if( binaryMap->GetPixel( imageIndex ) )
					{
						iter.SetPixel(neighIter.GetNeighborhoodIndex(), overlayVal);
					}
				}
			}			
		}
	}
	
	template <unsigned int VDimension>
	void
	PolyLineParametricPathExtended<VDimension>::	
	ComputeBoundingImageRegionFromVertexList(const VertexListType* vertexList,
																					 ImageRegionType& region) const
	{
		typedef BoundingBox<typename VertexListType::ElementIdentifier,
		Dimension, typename VertexType::CoordRepType, VertexListType> BoundingBoxType;
		
		IndexType				regionStartIndex;
		IndexType				regionEndIndex;
		SizeType				regionSize;
		
		if( vertexList->Size() == 0 )
		{
			itkWarningMacro(<<"ComputeBoundingImageRegionFromVertexList method "
											<<"is called with an empty vertex list. Returning an "
											<<"empty image region!");
			
			regionStartIndex.Fill( 0 );
			regionSize.Fill( 0 );
		}
		else
		{		
			// Compute the bounding box of the path.
			typename BoundingBoxType::Pointer boundingBox = BoundingBoxType::New();
			boundingBox->SetPoints( vertexList );
			boundingBox->ComputeBoundingBox();
			typename BoundingBoxType::PointType bbMin = boundingBox->GetMinimum();
			typename BoundingBoxType::PointType bbMax = boundingBox->GetMaximum();			
			for(unsigned int dim = 0; dim < Dimension; dim++)
			{
				// Note that a pixel at index indx takes the image 
				// region up to (indx + 0.5). That is why we use round here.
				regionStartIndex[dim] = 
				Math::Round<typename IndexType::IndexValueType>(bbMin[dim]);
				regionEndIndex[dim] = 
				Math::Round<typename IndexType::IndexValueType>(bbMax[dim]);
				
				regionSize[dim] = 
				static_cast<typename SizeType::SizeValueType>
				(regionEndIndex[dim] - regionStartIndex[dim] + 1);
			}
		}		
		
		region.SetIndex( regionStartIndex );
		region.SetSize( regionSize );
	}
	
	
	template <unsigned int VDimension>
	typename PolyLineParametricPathExtended<VDimension>::NMinus1DPathType::Pointer 
	PolyLineParametricPathExtended<VDimension>::	
	ConvertToNMinus1DPath(double radiusOrigin, double radiusSpacing) const
	{
		typename NMinus1DPathType::Pointer nMinus1DPath = NMinus1DPathType::New();
		nMinus1DPath->Initialize();
		
		const VertexListType* nDVertexList = this->GetVertexList();
		
		for(typename VertexListType::ElementIdentifier i = 0; 
				i < nDVertexList->Size(); 
				i++ )
		{
			typename NMinus1DPathType::VertexType vertex;
			VertexType nDVertex = nDVertexList->ElementAt(i);
			for(unsigned int j = 0; 
					j < NMinus1DPathType::Dimension; 
					j++ )
			{
				vertex[j] = nDVertex[j];
			}
			
			// Compute the point radius in world coordinates. 
			typename NMinus1DPathType::RadiusType radius = 
			nDVertex[Dimension-1] * radiusSpacing + radiusOrigin;
			
			// Add the (N-1)-D vertex with the Nth element treated 
			// as the index along the radius dimension.
			nMinus1DPath->AddVertex( vertex, radius );
		}
		
		return nMinus1DPath;
	}
	
	
	template <unsigned int VDimension>
	typename PolyLineParametricPathExtended<VDimension>::NPlus1DPathType::Pointer 
	PolyLineParametricPathExtended<VDimension>::	
	ConvertToNPlus1DPath(double radiusOrigin, double radiusSpacing) const
	{
		typename NPlus1DPathType::Pointer nPlus1DPath = NPlus1DPathType::New();
		nPlus1DPath->Initialize();
		
		const VertexListType* nDVertexList = this->GetVertexList();
		const RadiusListType& nDRadiusList = this->GetRadiusList();
		
		for(typename VertexListType::ElementIdentifier i = 0; 
				i < nDVertexList->Size(); 
				i++ )
		{
			typename NPlus1DPathType::VertexType vertex;
			VertexType nDVertex = nDVertexList->ElementAt(i);
			for(unsigned int j = 0; 
					j < Dimension; 
					j++ )
			{
				vertex[j] = nDVertex[j];
			}
			
			// Compute the point radius index in image coordinates. 
			vertex[NPlus1DPathType::Dimension-1] = (nDRadiusList[i] - radiusOrigin) / radiusSpacing;
			
			nPlus1DPath->AddVertex( vertex, nDRadiusList[i] );	// keep the original radius value.
		}
		
		return nPlus1DPath;
	}
	
	
	template <unsigned int VDimension>
	unsigned long 
	PolyLineParametricPathExtended<VDimension>::
	GetNoOfIntersectionPoints(const Self* path) const
	{
		Pointer pathClone1 = this->Clone();
		Pointer pathClone2 = path->Clone();
		
		pathClone1->RoundIndices(true);
		pathClone2->RoundIndices(true);
		
		typename VertexListType::ConstPointer vertexList1 =  pathClone1->GetVertexList();
		typename VertexListType::ConstPointer vertexList2 =  pathClone2->GetVertexList();
		
		unsigned long intersectionCount = 0;
		for(typename VertexListType::ElementIdentifier i = 0; 
				i < vertexList1->Size(); 
				i++)
		{
			for(typename VertexListType::ElementIdentifier j = 0; 
					j < vertexList2->Size(); 
					j++)
			{
				if(vertexList1->ElementAt(i).SquaredEuclideanDistanceTo(vertexList2->ElementAt(j)) <
					 m_Epsilon)
				{
					intersectionCount++;
					break;
				}
			}
		}
		
		return intersectionCount;
	}

	template <unsigned int VDimension>
	template<class TImage>
	void
	PolyLineParametricPathExtended<VDimension>::
	ComputeIntersectionAndUnionInImageCoords(const Self* path,
																					 const TImage* image,
																					 unsigned long& intersectionVolume,
																					 unsigned long& unionVolume,
																					 bool excludeTails) const
	{
		typedef TImage																ImageType;
		typedef ImageRegionConstIterator<ImageType>		RegionIterType;
		
		intersectionVolume = 0;
		unionVolume = 0;
		
		// Compute the union of bounding image regions.
		// Its damn stupid not to have an itk method for this 
		// inside the ImageRegion class.
		ImageRegionType unionOfBoundingRegions;
		ImageRegionType boundingRegion1 = 
		this->ComputeBoundingImageRegionWithRadius(image->GetSpacing());
		ImageRegionType boundingRegion2 = 
		path->ComputeBoundingImageRegionWithRadius(image->GetSpacing());
		boundingRegion1.Crop(image->GetBufferedRegion());
		boundingRegion2.Crop(image->GetBufferedRegion());
		typename ImageRegionType::SizeValueType noOfBBPixels1 = 
		boundingRegion1.GetNumberOfPixels();
		typename ImageRegionType::SizeValueType noOfBBPixels2 = 
		boundingRegion2.GetNumberOfPixels();
		typename ImageRegionType::SizeValueType zeroValue = 
		NumericTraits<typename ImageRegionType::SizeValueType>::ZeroValue();
		if((noOfBBPixels1 != zeroValue) && (noOfBBPixels2 != zeroValue))
		{
			IndexType topLeftIndex1 = boundingRegion1.GetIndex();
			IndexType topLeftIndex2 = boundingRegion2.GetIndex();
			IndexType bottomRightIndex1 = 
			boundingRegion1.GetIndex() + boundingRegion1.GetSize();
			IndexType bottomRightIndex2 = 
			boundingRegion2.GetIndex() + boundingRegion2.GetSize();
			for(unsigned int i = 0; i < Dimension; i++)
			{
				bottomRightIndex1[i]--;
				bottomRightIndex2[i]--;
			}
			if( !boundingRegion1.IsInside( topLeftIndex2 ) )
			{
				for(unsigned int i = 0; i < Dimension; i++)
				{
					if( topLeftIndex1[i] > topLeftIndex2[i] )
					{
						topLeftIndex1[i] = topLeftIndex2[i];
					}
				}
			}
			SizeType unionSize;
			for(unsigned int i = 0; i < Dimension; i++)
			{
				if( bottomRightIndex1[i] < bottomRightIndex2[i] )
				{
					unionSize[i] = bottomRightIndex2[i] - topLeftIndex1[i] + 1;
				}
				else
				{
					unionSize[i] = bottomRightIndex1[i] - topLeftIndex1[i] + 1;
				}
			}
			unionOfBoundingRegions.SetIndex( topLeftIndex1 );
			unionOfBoundingRegions.SetSize( unionSize );
		}
		else if((noOfBBPixels1 == zeroValue) && (noOfBBPixels2 != zeroValue))
		{
			unionOfBoundingRegions = boundingRegion2;
		}
		else if((noOfBBPixels1 != zeroValue) && (noOfBBPixels2 == zeroValue))
		{
			unionOfBoundingRegions = boundingRegion1;
		}
		else
		{
			return;
		}
		
		
		// Traverse the union of bounding regions and compute the 
		// number of pixels in the intersection and union regions.
		RegionIterType iter(image, unionOfBoundingRegions);
		for(iter.GoToBegin(); 
				!iter.IsAtEnd(); 
				++iter)
		{	
			// Find the closest path point to this pixel. 
			IndexType currentIndex = iter.GetIndex();
			OutputType currentContIndex( currentIndex );
			
			bool isInsidePath1 = 
			this->DoesContainPoint(currentContIndex, image, excludeTails);
			
			bool isInsidePath2 = 
			path->DoesContainPoint(currentContIndex, image, excludeTails);
			
			if(isInsidePath1 && isInsidePath2)
			{
				intersectionVolume++;
			}
			if(isInsidePath1 || isInsidePath2)
			{
				unionVolume++;
			}
		}
	}
																					 
	
	template <unsigned int VDimension>
	template<class TImage>
	void
	PolyLineParametricPathExtended<VDimension>::
	ComputeIntersectionAndUnionInWorldCoords(const Self* path,
																					 const TImage* image,
																					 double& intersectionVolume,
																					 double& unionVolume,
																					 bool excludeTails) const
	{
		// Compute first in image coordinates.
		unsigned long intersectPixelCounter;
		unsigned long unionPixelCounter;
		this->ComputeIntersectionAndUnionInImageCoords(path, image,
																									 intersectPixelCounter,
																									 unionPixelCounter,
																									 excludeTails);

		// Multiply the number of pixel with the volume of a pixel to 
		// get the actual volume of the regions.
		SpacingType spacing = image->GetSpacing();
		double pixelVolume = 1.0;
		for(unsigned int i = 0; i < Dimension; i++)
		{
			pixelVolume *= static_cast<double>(spacing[i]);
		}
		intersectionVolume = 
		static_cast<double>(intersectPixelCounter) * pixelVolume;
		unionVolume = 
		static_cast<double>(unionPixelCounter) * pixelVolume;
	}

	
	
	template <unsigned int VDimension>
	PolyLineParametricPathExtended<VDimension>
	::PolyLineParametricPathExtended()
	{
		this->SetDefaultInputStepSize( 0.1 );
		this->SetEpsilon( NumericTraits<float>::epsilon() );
		m_VertexList = VertexListType::New();
		m_EffectiveDimension = Dimension;
	}
	
	template <unsigned int VDimension>
	void
	PolyLineParametricPathExtended<VDimension>
	::PrintSelf( std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf( os, indent );
		os << indent << "Verticies:  " << m_VertexList << std::endl;
		os << indent << "Radius List:  " << &m_RadiusList << std::endl;
		
		os << indent << "Effective Dimension:  " << m_EffectiveDimension << std::endl;
		
		os << indent << "Epsilon:  " << m_Epsilon << std::endl;
	}
	
} // end namespace itk