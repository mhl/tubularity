//**********************************************************
//Copyright 2010 Engin Turetken
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


#ifndef __itkHessianMainPrincipleCurvatureObjectnessImageFilter_h
#define __itkHessianMainPrincipleCurvatureObjectnessImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "itkSymmetricEigenAnalysis.h"
#include "itkArray.h"

namespace itk
{
	
	/** This functor computes the principle curvature for every input symmetric
	 * matrix/tensor pixel from its eigen values. The principle curvature measure 
	 * is a variant of the Lindeberg's tubularity measure (Lindeberg: "Edge 
	 * detection and ridge detection with automatic scale selection", 
	 * International Journal of Computer Vision, vol 30, number 2, pp. 
	 * 117–154, 1998.) In essence, for bright(dark) structures, it computes 
	 * the normalized absolute value of the minimum(maximum) eigenvalue 
	 * of the hessian matrix. 
	 *
	 * The input pixel type must provide the API for the [][] operator and 
	 * also for the definition of 'ScalarRealType' type.
	 * 
	 * \author Engin Turetken
	 *
	 */
	namespace Functor 
	{  	
		template< typename TInput, typename TOutput >
		class HessianMainPrincipleCurvatureObjectnessFunction
		{
		public:
			typedef typename NumericTraits<TInput>::ScalarRealType	ScalarRealType;
			
			HessianMainPrincipleCurvatureObjectnessFunction()
			{
				// Set the default values of the attributes.
				m_BrightObject = true;
				m_InputTensorSize = 2;
				m_Sigma = 1.0;
				m_Gamma = 1.5;
			}
			~HessianMainPrincipleCurvatureObjectnessFunction() {}
			
			bool operator!=( const HessianMainPrincipleCurvatureObjectnessFunction & ) const
			{
				return false;
			}
			bool operator==( const HessianMainPrincipleCurvatureObjectnessFunction & other ) const
			{
				return !(*this != other);
			}
			
			inline TOutput operator()( const TInput & input ) const
			{
				typedef Array<ScalarRealType>								EigenValArrayType;
				typedef SymmetricEigenAnalysis< TInput, EigenValArrayType>	CalculatorType;
				
				// Create the eigen value calculator, the eigenvalue array and the output buffer.
				CalculatorType calculator(m_InputTensorSize);
				EigenValArrayType eigenValues(m_InputTensorSize);
				TOutput output;
				
				// Order the eigen values in ascending order.
				// This is, lambda_1 < lambda_2 < ....
				calculator.SetOrderEigenValues( true );
				
				//Compute the objectness measure.
				calculator.ComputeEigenValues( input, eigenValues );
				if( m_BrightObject )
				{
					output = static_cast<TOutput>(0.5 * 
																				vcl_pow(m_Sigma, m_Gamma) * 
																				vcl_fabs(eigenValues[0]));
				}
				else
				{
					output = static_cast<TOutput>(0.5 * 
																				vcl_pow(m_Sigma, m_Gamma) * 
																				vcl_fabs(eigenValues[m_InputTensorSize-1]));
				}
				
				return output;
			}
			
			/** Set the size of the input pixel tensor. By default, it is 2. */
			void SetInputTensorSize(unsigned int InputTensorSize)
			{
				m_InputTensorSize = InputTensorSize;
			}
			
			/** Get the size of the input pixel tensor. */
			unsigned int GetInputTensorSize() const
			{
				return m_InputTensorSize;
			}
			
			
			/** Method to set the bright object flag.  */
			void SetBrightObject( bool bIsBrightObject )
			{
				m_BrightObject = bIsBrightObject;
			}
			
			/** Method to get the bright object flag.  */
			bool GetBrightObject() const
			{
				return m_BrightObject;
			}
			
			/** Method to set the sigma value of the Gaussian 
			 * function that is used to compute the input Hessian 
			 * matrix pixels. By default, it is one and hence 
			 * no normalization is carried out for computing 
			 * the objectness measure. */
			void SetSigma( ScalarRealType Sigma )
			{
				m_Sigma = Sigma;
			}
			
			/** Method to get the sigma value of the Gaussian 
			 * function that is used to compute the input Hessian 
			 * matrix pixels. By default, it is one. */
			ScalarRealType GetSigma() const
			{
				return m_Sigma;
			}
			
			/** Method to set the gamma value of the normalization 
			 * constant. By default, it is 1.5 . */
			void SetGamma( ScalarRealType Gamma )
			{
				m_Gamma = Gamma;
			}
			
			/** Method to get the gamma value of the normalization 
			 * constant. By default, it is 1.5 . */
			ScalarRealType GetGamma() const
			{
				return m_Gamma;
			}
			
		private:
			unsigned int m_InputTensorSize;
			bool m_BrightObject;
			ScalarRealType m_Sigma;
			ScalarRealType m_Gamma;
		}; 
		
	}  // end namespace functor
	
	
	/** \class HessianMainPrincipleCurvatureObjectnessImageFilter
	 * \brief Computes the principle curvature for every input symmetric matrix pixel.
	 *
	 * The principle curvature measure is a variant of the Lindeberg's 
	 * tubularity measure (Lindeberg: "Edge detection and ridge detection with 
	 * automatic scale selection", International Journal of Computer Vision, 
	 * vol 30, number 2, pp. 117–154, 1998.) In essence, for bright(dark) 
	 * structures, it computes the normalized absolute value of the 
	 * minimum(maximum) eigenvalue of the hessian matrix. 
	 * 
	 * 
	 * \author Engin Turetken 
	 *
	 * \ingroup  Multithreaded  TensorObjects
	 */
	template <typename  TInputImage, 
	typename  TOutputImage = itk::Image<typename NumericTraits<
	typename TInputImage::PixelType::ValueType>::RealType, 
	TInputImage::ImageDimension > >
	class ITK_EXPORT HessianMainPrincipleCurvatureObjectnessImageFilter :
	public
	UnaryFunctorImageFilter<TInputImage,TOutputImage, 
	Functor::HessianMainPrincipleCurvatureObjectnessFunction< 
	typename TInputImage::PixelType,
	typename TOutputImage::PixelType> >
	{
	public:
		/** Standard class typedefs. */
		typedef HessianMainPrincipleCurvatureObjectnessImageFilter  Self;
		typedef UnaryFunctorImageFilter<
		TInputImage,TOutputImage, 
		Functor::HessianMainPrincipleCurvatureObjectnessFunction< 
		typename TInputImage::PixelType,
		typename TOutputImage::PixelType> >				Superclass;
		typedef SmartPointer<Self>											Pointer;
		typedef SmartPointer<const Self>								ConstPointer;
		
		typedef TOutputImage														OutputImageType;
		typedef TInputImage															InputImageType;		
		typedef typename OutputImageType::PixelType     OutputPixelType;
		typedef typename InputImageType::PixelType      InputPixelType;
		typedef typename InputPixelType::ValueType      InputValueType;
		typedef typename Superclass::FunctorType        FunctorType; 
		
		typedef typename NumericTraits<InputPixelType>::ScalarRealType  ScalarRealType;
		
		/** Run-time type information (and related methods).   */
		itkTypeMacro( HessianMainPrincipleCurvatureObjectnessImageFilter, UnaryFunctorImageFilter );
		
		/** Method for creation through the object factory. */
		itkNewMacro(Self);
		
		/** Print internal ivars */
		void PrintSelf(std::ostream& os, Indent indent) const
		{ 
			this->Superclass::PrintSelf( os, indent ); 
			// TODO
		}
		
		/** Set the size of the input pixel tensor. By default, it is set to 
		 * the input image dimension. */
		void SetInputTensorSize(unsigned int InputTensorSize)
		{
			if( this->GetFunctor().GetInputTensorSize() != InputTensorSize )
			{
				this->GetFunctor().SetInputTensorSize(InputTensorSize);
				this->Modified();
			}
		}
		
		/** Get the size of the input pixel tensor. */
		unsigned int GetInputTensorSize() const
		{
			return this->GetFunctor().GetInputTensorSize();
		}
		
		/** Method to set the sigma value of the Gaussian 
		 * function that is used to compute the input Hessian 
		 * matrix pixels. By default, it is one (i.e., 
		 * no smoothing) and hence no normalization is carried 
		 * out for computing the objectness measure. */
		void SetSigma( ScalarRealType Sigma )
		{
			if( this->GetFunctor().GetSigma() != Sigma )
			{
				this->GetFunctor().SetSigma(Sigma);
				this->Modified();
			}
		}
		
		/** Method to get the sigma value of the Gaussian 
		 * function that is used to compute the input Hessian 
		 * matrix pixels. By default, it is one. */
		ScalarRealType GetSigma() const
		{
			return this->GetFunctor().GetSigma();
		}
		
		
		/** Method to set the gamma value of the normalization 
		 * constant. By default, it is 1.5 . */
		void SetGamma( ScalarRealType Gamma )
		{
			if( this->GetFunctor().GetGamma() != Gamma )
			{
				this->GetFunctor().SetGamma(Gamma);
				this->Modified();
			}
		}
		
		/** Method to get the gamma value of the normalization 
		 * constant. By default, it is 1.5 . */
		ScalarRealType GetGamma() const
		{
			return this->GetFunctor().GetGamma();
		}
		
		/** Set the bright object flag. */		
		void SetBrightObject( bool bIsBrightObject )
		{
			if( this->GetFunctor().GetBrightObject() != bIsBrightObject )
			{
				this->GetFunctor().SetBrightObject(bIsBrightObject);
				this->Modified();
			}
		}
		
		/** Get the bright object flag. */				
		bool GetBrightObject() const
		{
			return this->GetFunctor().GetBrightObject();
		}
		
#ifdef ITK_USE_CONCEPT_CHECKING
		/** Begin concept checking */
		itkConceptMacro(InputHasNumericTraitsCheck,
										(Concept::HasNumericTraits<InputValueType>));
		/** End concept checking */
#endif
		
	protected:
		HessianMainPrincipleCurvatureObjectnessImageFilter()
		{
			this->GetFunctor().SetInputTensorSize( TInputImage::ImageDimension );
			this->GetFunctor().SetSigma( 1.0 );
			this->GetFunctor().SetGamma( 1.5 );
		}
		virtual ~HessianMainPrincipleCurvatureObjectnessImageFilter() {};
		
	private:
		HessianMainPrincipleCurvatureObjectnessImageFilter(const Self&); //purposely not implemented
		void operator=(const Self&);									 //purposely not implemented
		
	};
} // end namespace itk

#endif
