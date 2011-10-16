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

#ifndef __itkCircularShapedNeighborhoodIterCreator_h
#define __itkCircularShapedNeighborhoodIterCreator_h


#include "itkShapedNeighborhoodIterator.h"
#include "itkImage.h"
#include "itkVector.h"

namespace itk
{
	/** \class CircularShapedNeighborhoodIterCreator
	 * \brief  Creates a circular shaped neighborhood iterator from 
	 * a radius value given either in image or world coordinates.
	 * 
	 * The template argument TShapedNeighborhoodIter can be either 
	 * ConstShapedNeighborhoodIterator or ShapedNeighborhoodIterator.
	 * 
	 * \author Engin Turetken
	 *
	 */
	template <class TShapedNeighborhoodIter>
	class CircularShapedNeighborhoodIterCreator
	{
	public:
		/** Standard class typedefs. */
		typedef CircularShapedNeighborhoodIterCreator   Self;	
		
		/** Basic data-structure types used */
		typedef TShapedNeighborhoodIter								IteratorType;
		typedef typename IteratorType::ImageType			ImageType;
		typedef typename IteratorType::OffsetType			IterOffsetType;
		typedef typename IteratorType::SizeType				IterSizeType;
		typedef typename IterSizeType::SizeValueType	IterSizeValueType;
		typedef typename IteratorType::RegionType			IterRegionType;
		
		typedef double																CircleRadiusType;
		
		typedef itk::Vector<double, 
		ImageType::ImageDimension>										DirectionVectorType;
		
		/** Dimension constant. */
		itkStaticConstMacro(Dimension, unsigned int, ImageType::ImageDimension);
		
		/** Constructor and Destructor. */		
		CircularShapedNeighborhoodIterCreator();
		virtual ~CircularShapedNeighborhoodIterCreator(){}
			
		/** Setter/Getter for the radius value. */
		virtual void SetRadius( CircleRadiusType radius )
		{
			m_Radius = radius;
		}
		virtual CircleRadiusType GetRadius() const
		{
			return m_Radius;
		}
		
		/** Setter/Getter for the coordinate system of the radius value. */
		virtual void SetRadiusInImageCoordinates( bool flag )
		{
			m_RadiusInImageCoordinates = flag;
		}
		virtual bool GetRadiusInImageCoordinates() const
		{
			return m_RadiusInImageCoordinates;
		}
		
		/** Setter/Getter for the input image. */
		virtual void SetImage( const ImageType* image )
		{
			m_Image = image;
		}
		virtual const ImageType* GetImage() const
		{
			return m_Image.GetPointer();
		}
		
		/** Setter/Getter for the image region of the iterator. */
		virtual void SetRegion( const IterRegionType& region )
		{
			m_Region = region;
		}
		virtual const IterRegionType& GetRegion() const
		{
			return m_Region;
		}
		
		/** 
		 * Setter/Getter for the flag that determines if a 
		 * half-circular neighborhood will be generated 
		 * instead of a full one. 
		 */
		virtual void SetGenerateHalfCircularNeigh( bool flag )
		{
			m_GenerateHalfCircularNeigh = flag;
		}
		virtual bool GetGenerateHalfCircularNeigh() const
		{
			return m_GenerateHalfCircularNeigh;
		}

		/** 
		 * Setter/Getter for the direction of the half circle 
		 * in image coordinates. 
		 */
		virtual void SetHalfCircleDirection( const DirectionVectorType& direction )
		{			
			m_HalfCircleDirection = direction;
		}
		virtual const DirectionVectorType& GetHalfCircleDirection() const
		{
			return m_HalfCircleDirection;
		}
		
		
		virtual void GenerateIterator( IteratorType& neighIter );
		
	private:
		CircularShapedNeighborhoodIterCreator(const Self&); //purposely not implemented
		void operator=(const Self&);												//purposely not implemented		
		
		// Circular neighborhood radius. By default it is 1.0.
		CircleRadiusType									m_Radius;
		
		// Is the radius value given in image or world coordinate system.
		// By default, it is false.
		bool															m_RadiusInImageCoordinates;
		
		// Image region for the iterator to be generated. By default its 
		// index and size values are all zero.
		IterRegionType										m_Region;
		
		// Iterator image.
		typename ImageType::ConstPointer	m_Image;
		
		// Flag that determines if a half-circular neighborhood 
		// will be generated instead of a full one. By default, 
		// it is false, i.e., full circular neighborhood is generated.
		bool															m_GenerateHalfCircularNeigh;
		
		// Direction of the half circle in image coordinates. This member 
		// is used only if the m_GenerateHalfCircularNeigh flag is turned on.
		// By default, it is in the direction of the first image axis.
		DirectionVectorType								m_HalfCircleDirection;
	};
	
}

#if ITK_TEMPLATE_TXX
# include "itkCircularShapedNeighborhoodIterCreator.txx"
#endif

#endif
