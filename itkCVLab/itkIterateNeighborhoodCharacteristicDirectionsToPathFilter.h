//**********************************************************
//Copyright 2011 Fethallah Benmansour
//
//Licensed under the Apache License, Version 2.0 (the "License"); 
//you may not use this file except in compliance with the License. 
//You may obtain a copy of the License at
//
//http://www.apache.org/licenses/LICENSE-2.0 
//
//Unless required by applicable law or agreed to in writing, software 
//distributed under the License is distributed on an "AS IS" BASIS, 
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
//See the License for the specific language governing permissions and 
//limitations under the License.
//**********************************************************

#ifndef __itkIterateNeighborhoodCharacteristicDirectionsToPathFilter_h
#define __itkIterateNeighborhoodCharacteristicDirectionsToPathFilter_h

#include "itkNumericTraits.h"
#include "itkExceptionObject.h"
#include "itkContinuousIndex.h"
#include "itkPolyLineParametricPath.h"
#include "itkCommand.h"
#include "itkImageToPathFilter.h"
#include "itkSingleValuedNonLinearOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkNeighborhoodIterator.h"



namespace itk
{
	
	
	/** \class IterateNeighborhoodCharacteristicDirectionsToPathFilter
	 * \brief Extracts a path from Characteristic Directions estimated using
	 * the Eikonal solver.
	 *
	 * This filter extracts the geodesic (minimal) path between the given
	 * end-point and a start-point (which is implicitly embedded in the
	 * given arrival function). The path is extracted by back-propagating
	 * following the steepest direction given by the Eikonal solver
	 * from the end-point to the global minimum of the arrival function
	 * (ie. the start-point). 
	 * This is done just by iterativelly looking for the neighboor with
	 *  the smallest value
	 *
	 * The user must provide the following:
	 *    1. the Geodesic distance from the source
	 *    2. At least one path end point
	 *			 (AddEndPoint() should be called at least once).
	 *
	 * This optimizer is designed to operate on a monotonic cost function
	 * WITHOUT using gradient information (derivatives). The user must set
	 * the Neighborhood size, and optionally the connectivity.
	 *
	 * \author Fethallah Benmansour, CVLAB EPFL, fethallah[at]gmail.com
	 *
	 * \ingroup ImageToPathFilters
	 */
	
	template <class TInputImage,
	class TOutputPath = PolyLineParametricPath<TInputImage::ImageDimension> >
	class ITK_EXPORT IterateNeighborhoodCharacteristicDirectionsToPathFilter :
	public ImageToPathFilter< TInputImage, TOutputPath >
	{
	public:
		/** Standard class typedefs. */
		typedef IterateNeighborhoodCharacteristicDirectionsToPathFilter                Self;
		typedef ImageToPathFilter<TInputImage,TOutputPath>					Superclass;
		typedef SmartPointer<Self>																	Pointer;
		typedef SmartPointer<const Self>														ConstPointer;
		
		/** Run-time type information (and related methods). */
		itkTypeMacro( IterateNeighborhoodCharacteristicDirectionsToPathFilter, ImageToPathFilter );
		
		/** Method for creation through the object factory. */
		itkNewMacro(Self);
		
		/** ImageDimension constants */
		itkStaticConstMacro(SetDimension, unsigned int,
												TInputImage::ImageDimension);
		
		/** Some image typedefs. */
		typedef TInputImage																					InputImageType;
		typedef typename InputImageType::Pointer										InputImagePointer;
		typedef typename InputImageType::ConstPointer								InputImageConstPointer;
		typedef typename InputImageType::RegionType									InputImageRegionType; 
		typedef typename InputImageType::PixelType									InputImagePixelType;
		typedef typename InputImageType::IndexType									InputImageIndexType;
		typedef typename InputImageType::SpacingType								SpacingType;
		typedef InputImagePixelType																	VectorType;
		
		/** Some path typedefs. */
		typedef TOutputPath OutputPathType;
		typedef typename OutputPathType::Pointer										OutputPathPointer;
		typedef typename OutputPathType::ConstPointer								OutputPathConstPointer;
		
		/** Some convenient typedefs. */
		typedef Index< SetDimension >																IndexType;
		typedef ContinuousIndex< double, SetDimension >							ContinuousIndexType;
		typedef Point< double, SetDimension >												PointType;
		
		
		
		void SetStartPoint(IndexType startPointIndex)
		{
			PointType startPoint;
			this->GetInput()->TransformIndexToPhysicalPoint( startPointIndex, startPoint );
			
			m_StartPoint = startPoint;
		}
		
		/** Clears the list of end points and adds the given point to the list. */
		virtual void SetPathEndPoint( const IndexType & point )
		{
			this->ClearPathEndPoints();
			this->AddPathEndPoint( point );
		}
		
		/** Adds the given point to the list. */
		virtual void AddPathEndPoint( const IndexType & index )
		{
			PointType point;
			this->GetInput()->TransformIndexToPhysicalPoint( index, point );
			m_EndPointList.push_back( point );
			this->Modified();
		};
		
		/** Clear the list of end points. */
		virtual void ClearPathEndPoints()
		{
			if (m_EndPointList.size() > 0)
      {
				m_EndPointList.clear();
				this->Modified();
      }
		};
		
		itkSetMacro(NbMaxIter, unsigned int);
		
	protected:
		IterateNeighborhoodCharacteristicDirectionsToPathFilter();
		~IterateNeighborhoodCharacteristicDirectionsToPathFilter();
		virtual void PrintSelf(std::ostream& os, Indent indent) const;
		
		/** Override since the filter needs all the data for the algorithm */
		void GenerateInputRequestedRegion();
		
		/** Implemention of algorithm */
		void GenerateData(void);
		
		/** Get the arrival function from which to extract the path. */
		virtual unsigned int GetNumberOfPathsToExtract( ) const;
		
		/** Get the next end point from which to back propagate. */
		virtual const PointType & GetNextEndPoint( );
	private:
		IterateNeighborhoodCharacteristicDirectionsToPathFilter(const Self&); //purposely not implemented
		void operator=(const Self&); //purposely not implemented
		
		unsigned int																								m_NbMaxIter;
		std::vector<PointType>																			m_EndPointList;
		PointType																										m_StartPoint;
		unsigned int																								m_CurrentOutput;
		
		InputImageIndexType																					m_StartIndex;
		InputImageIndexType																					m_LastIndex;
		
	};
	
}

#include "itkIterateNeighborhoodCharacteristicDirectionsToPathFilter.txx"

#endif