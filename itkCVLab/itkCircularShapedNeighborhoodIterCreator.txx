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


#include "itkCircularShapedNeighborhoodIterCreator.h"
#include "itkMath.h"
#include "itkNumericTraits.h"

#include <vcl_algorithm.h>

namespace itk
{
	template<class TShapedNeighborhoodIter>
	CircularShapedNeighborhoodIterCreator<TShapedNeighborhoodIter>
	::CircularShapedNeighborhoodIterCreator()
	{
		m_Radius = 1.0;
		m_RadiusInImageCoordinates = false;
		m_GenerateHalfCircularNeigh = false;
		m_HalfCircleDirection.Fill( 0.0 );
		m_HalfCircleDirection[0] = 1.0;
	}
	
	template<class TShapedNeighborhoodIter>
	void
	CircularShapedNeighborhoodIterCreator<TShapedNeighborhoodIter>
	::GenerateIterator( IteratorType& neighIter )
	{
		// Check if the input image is null.
		if( m_Image.IsNull() )
		{
			itkGenericExceptionMacro(<<"Input image is null!");
		}
		
		if( m_GenerateHalfCircularNeigh )
		{
			if( m_HalfCircleDirection.GetNorm() < itk::NumericTraits
					<typename DirectionVectorType::RealValueType>::epsilon() )
			{
				itkGenericExceptionMacro(<<"The input half-circle direction vector is all zero.");
			}
		}
		
		// Check if the region is empty.
		if(m_Region.GetNumberOfPixels() == itk::NumericTraits
			 <typename IterRegionType::SizeValueType>::ZeroValue() )
		{
			itkGenericOutputMacro(<<"The region is empty, possibly not set "
														<<"by the caller. Setting it to the "
														<<"BufferedRegion of the input image.");
			
			m_Region = m_Image->GetBufferedRegion();
		}
		else
		{
			// Crop the region by the buffered possible region of the image.
			bool isInside = m_Region.Crop(m_Image->GetBufferedRegion());
			
			// Check if the region is outside the buffered 
			// region of the image.
			if( !isInside )
			{
				itkGenericOutputMacro(<<"The region is outside the  "
															<<"BufferedRegion of the input image. "
															<<"Setting it to the "
															<<"BufferedRegion of the image.");
				
				m_Region = m_Image->GetBufferedRegion();
			}
		}
		
		// Get the image spacing.
		typename ImageType::SpacingType spacing = m_Image->GetSpacing();
		
		IterSizeType size;
		typename ImageType::SpacingType spacingSqr;
		if( m_RadiusInImageCoordinates )
		{	
			size.Fill( Math::Ceil<IterSizeValueType>(m_Radius) );
			
			// Do not take into account spacing if the radius is 
			// already in image coordinates.
			spacingSqr.Fill( itk::NumericTraits
											<typename ImageType::SpacingValueType>::OneValue() );
		}		
		else
		{	
			// Find the minimum spacing and the corresponding maximum number of 
			// pixels that this radius corresponds to in the image coordinates.
			// This is required to correctly set the radius of the neighborhood 
			// iterator.
			IterSizeValueType maxPossibleIntegerRadius = 
			NumericTraits<IterSizeValueType>::ZeroValue();
			for(unsigned int dim = 0; dim < Dimension; dim++)
			{
				IterSizeValueType dimRadius = 
				vcl_min(Math::Ceil<IterSizeValueType>(m_Radius / spacing[dim]),
								static_cast<IterSizeValueType>(m_Region.GetSize()[dim]));
				
				if( maxPossibleIntegerRadius < dimRadius )
				{
					maxPossibleIntegerRadius = dimRadius;
				}
				
				spacingSqr[dim] = spacing[dim] * spacing[dim];
			}
			
			size.Fill( maxPossibleIntegerRadius );
		}	
		IteratorType dummyNeighIter(size, m_Image, m_Region);
		neighIter = dummyNeighIter;
		neighIter.ClearActiveList();
		
		
		CircleRadiusType epsilon = NumericTraits<CircleRadiusType>::epsilon();
		CircleRadiusType radiusSqrThrsh = m_Radius * m_Radius + epsilon;
		unsigned int centerPixelIndex = neighIter.GetCenterNeighborhoodIndex();
		unsigned int totalNumberOfPixels = 2 * centerPixelIndex + 1;
		for(unsigned int indx = 0; indx < totalNumberOfPixels; indx++)
		{
			IterOffsetType offset = neighIter.GetOffset( indx );
			
			// Compute the Euclidean distance of this point from the center one 
			// in image coordindates.
			CircleRadiusType distSqr = 0;
			
			if( !m_GenerateHalfCircularNeigh )
			{
				for( unsigned int dim = 0; dim < Dimension; dim++ )
				{
					distSqr += spacingSqr[dim] * offset[dim] * offset[dim];
				}
				if( !(distSqr > radiusSqrThrsh ) )
				{
					neighIter.ActivateOffset( offset );
				}
			}
			else
			{
				typename DirectionVectorType::RealValueType dotProd = 0;
				for( unsigned int dim = 0; dim < Dimension; dim++ )
				{
					distSqr += spacingSqr[dim] * offset[dim] * offset[dim];
					dotProd += offset[dim] * m_HalfCircleDirection[dim] * spacingSqr[dim]; 
					// Note that, in the above line, spacingSqr essentially converts 
					// both vectors to world coordinate system. 
				}
				if( (!(distSqr > radiusSqrThrsh )) && (dotProd > -epsilon) )
				{
					neighIter.ActivateOffset( offset );
				}			
			}
		}
	}
	
} // end namespace itk