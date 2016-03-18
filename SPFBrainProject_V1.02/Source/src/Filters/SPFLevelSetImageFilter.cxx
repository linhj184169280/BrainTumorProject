#ifndef SPFLEVELSETIMAGEFILTER_HXX_INCLUDED
#define SPFLEVELSETIMAGEFILTER_HXX_INCLUDED


#include "SPFLevelSetImageFilter.h"

#include "itkNumericTraits.h"
//#include "itkProgressReporter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkNeighborhoodAlgorithm.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"


namespace MySpace
{

SPFLevelSetImageFilter::
SPFLevelSetImageFilter()
{
    //Initial the filter..
    m_CropFilter = CropFilterType::New();
    m_IsotropicResampler = IsotropicResamplerType::New();
    m_ResampleThickSliceData = true;
    m_AnisotropyThreshold = 1.0;

    this->m_MaximumNumberOfIterations = 20;

    this->m_Sigma = 0.3;
    this->m_NeiborSize = 3;
    this->m_delta = 1.0;
    this->m_Propagation= 10.0;

}


void SPFLevelSetImageFilter
::GenerateInputRequestedRegion() throw(InvalidRequestedRegionError)
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  if ( !this->GetInput() )
    {
      InputImageType::Pointer inputPtr  =
      const_cast< ImageType *>( this->GetInput() );

    // Request the entire input image
    inputPtr->SetRequestedRegion(inputPtr->GetLargestPossibleRegion());
    }
}



void SPFLevelSetImageFilter
::GenerateOutputInformation()
{
  // get pointers to the input and output
    Superclass::OutputImagePointer      outputPtr = this->GetOutput();
    Superclass::InputImageConstPointer  inputPtr  = this->GetInput();
  if ( !outputPtr || !inputPtr)
    {
    return;
    }

  // Minipipeline is :
  //   Input -> Crop -> Resample_if_too_anisotropic -> Segment

  m_CropFilter->SetInput(inputPtr);
  m_CropFilter->SetRegionOfInterest(m_RegionOfInterest);

  // Compute the spacing after isotropic resampling.
  double minSpacing = NumericTraits< double >::max();
  for (int i = 0; i < ImageDimension; i++)
    {
    minSpacing = (minSpacing > inputPtr->GetSpacing()[i] ?
                  inputPtr->GetSpacing()[i] : minSpacing);
    }

  // Try and reduce the anisotropy.
  SpacingType outputSpacing = inputPtr->GetSpacing();
  for (int i = 0; i < ImageDimension; i++)
    {
    if (outputSpacing[i]/minSpacing > m_AnisotropyThreshold && m_ResampleThickSliceData)
      {
      outputSpacing[i] = minSpacing * m_AnisotropyThreshold;
      }
    }

  if (m_ResampleThickSliceData)
    {
    m_IsotropicResampler->SetInput( m_CropFilter->GetOutput() );
    m_IsotropicResampler->SetOutputSpacing( outputSpacing );
    m_IsotropicResampler->GenerateOutputInformation();
    outputPtr->CopyInformation( m_IsotropicResampler->GetOutput() );
    }
  else
    {
    outputPtr->CopyInformation( m_CropFilter->GetOutput() );
    }
}




void SPFLevelSetImageFilter
::GenerateData()
{


  // Allocate the output
  this->GetOutput()->SetBufferedRegion( this->GetOutput()->GetRequestedRegion() );
  this->GetOutput()->Allocate();

  // Get the input image
    InputImageType::ConstPointer  input  = this->GetInput();

  // Crop and perform thin slice resampling (done only if necessary)
  m_CropFilter->Update();

    InputImageType::Pointer inputImage = NULL;
  if (m_ResampleThickSliceData)
    {
    m_IsotropicResampler->Update();
    inputImage = this->m_IsotropicResampler->GetOutput();
    }
  else
    {
    inputImage = m_CropFilter->GetOutput();
    }

  // Convert the output of resampling (or cropping based on
  // m_ResampleThickSliceData) to a spatial object that can be fed into
  // the lesion segmentation method
    inputImage->DisconnectPipeline();

    OutputImageType *featureImage = inputImage;

    typedef DiscreteGaussianImageFilter<
		OutputImageType, OutputImageType >		  GaussianFilterType;

    //Get Seed Point...
      SeedSpatialObjectType::Pointer seedSpatialObject =
        SeedSpatialObjectType::New();
    seedSpatialObject->SetPoints(m_Seeds);
	const unsigned int numberOfPoints = seedSpatialObject->GetNumberOfPoints();
	typedef   InputSpatialObjectType::PointListType		PointListType;
	const PointListType &points = seedSpatialObject->GetPoints();

    typedef   OutputImageType::IndexType					IndexType;
	IndexType  index;

    //Output
    outputImage = OutputImageType::New();
	outputImage->SetRegions( featureImage->GetRequestedRegion() );
	outputImage->CopyInformation( featureImage );
	outputImage->Allocate();

	//Gradient Map
      InternalImageType::Pointer  gradientImage = InternalImageType::New();
	gradientImage->SetRegions( featureImage->GetRequestedRegion() );
	gradientImage->CopyInformation( featureImage );
	gradientImage->Allocate();

	//SPF Map
      InternalImageType::Pointer  SPFImage = InternalImageType::New();
	SPFImage->SetRegions( featureImage->GetRequestedRegion() );
	SPFImage->CopyInformation( SPFImage );
	SPFImage->Allocate();

	//Declare of Feature Image Iterator
	ImageRegionIterator< FeatureImageType >
		featureIt( featureImage, featureImage->GetRequestedRegion() );

	//Declare of SPF Image Iterator
	ImageRegionIterator< InternalImageType >
		SPFIt( SPFImage, SPFImage->GetRequestedRegion() );

	//Declare of Output Image Iterator
	ImageRegionIterator< OutputImageType >
		outputIt( outputImage, outputImage->GetRequestedRegion() );

    //Initial the contour within Neighbor distance 2
    for(outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt )
	{
		outputIt.Set(-1.0);
	}
	for( unsigned int i=0; i< numberOfPoints; i++ )
	{
		  NeighborhoodIterator< OutputImageType >::RadiusType  radius;
		radius.Fill(2);				//Neighbor with in 2
		NeighborhoodIterator< OutputImageType >
			outputNeiborIt( radius, outputImage, outputImage->GetRequestedRegion() );

		outputImage->TransformPhysicalPointToIndex( points[i].GetPosition(), index );

		outputNeiborIt.SetLocation( index );

		for(int i=0; i<outputNeiborIt.Size();i++)
		{
			outputNeiborIt.SetPixel(i,1.0);
		}
	}

    unsigned int iteratorNum, changedNumber;
	iteratorNum=changedNumber=0;		//Initial the iterator Number

	//Start algorithm of SPF
    while( iteratorNum<=m_MaximumNumberOfIterations )
	{
	    //Initial the Gradient Filter..
        typedef GradientMagnitudeImageFilter< OutputImageType, OutputImageType > GradientFilterType;
		  GradientFilterType::Pointer m_gradientFilter = GradientFilterType::New();

        //Initial the Gaussian Filter..
        typedef DiscreteGaussianImageFilter< InternalImageType, OutputImageType > GaussianFilterType;
		  GaussianFilterType::Pointer m_gaussianFilter = GaussianFilterType::New();
		m_gaussianFilter->SetUseImageSpacingOn();
	//	m_gaussianFilter->SetMaximumKernelWidth( m_NeiborSize );
		m_gaussianFilter->SetVariance( m_Sigma );

		//Calculate the gradient of level set map..
        m_gradientFilter->SetInput(outputImage);
		m_gradientFilter->Update();
		gradientImage = m_gradientFilter->GetOutput();

		//Initial the Iterator of Gradient map
        ImageRegionIterator< InternalImageType >
			gradientIt( gradientImage, gradientImage->GetRequestedRegion() );

        float c1,c2;		//c1,c2 store the average intensity value of points inside or outside the contour
		c1=c2=0;

        float c1_OnNum, c1_DownNum, c2_OnNum, c2_DownNum;
		c1_OnNum=c1_DownNum=c2_OnNum=c2_DownNum=0;	//temp parameters.

        //Calculate the average value of the inside and outside
        ImageRegionIterator< OutputImageType >
			outputIt( outputImage, outputImage->GetRequestedRegion() );
		for( outputIt.GoToBegin(), featureIt.GoToBegin() ; !outputIt.IsAtEnd(); ++outputIt, ++featureIt )
		{
			if( outputIt.Value()<0 )
			{
			//	c1_OnNum+= outputIt.Value()*featureIt.Value();
			//	c1_DownNum+= outputIt.Value();
				c1_OnNum+= featureIt.Value();
				c1_DownNum++;
			}
			else
			{
			//	c2_OnNum+= outputIt.Value()*featureIt.Value();
			//	c2_DownNum+= outputIt.Value();
				c2_OnNum+= featureIt.Value();
				c2_DownNum++;
			}
		}
        c1= c1_OnNum/c1_DownNum;
		c2= c2_OnNum/c2_DownNum;

        //Calculate the SPF Function value & update the level set map(OutputImage)..
        float max_SPF = 0;
		for( SPFIt.GoToBegin(), featureIt.GoToBegin() ; !SPFIt.IsAtEnd(); ++SPFIt, ++featureIt )
		{
			SPFIt.Set( featureIt.Value()-(c1+c2)/2 );
			max_SPF = max_SPF>abs(SPFIt.Value())?max_SPF:abs(SPFIt.Value());
		}
        for( gradientIt.GoToBegin(), SPFIt.GoToBegin(), outputIt.GoToBegin(); !gradientIt.IsAtEnd(); ++SPFIt, ++gradientIt, ++outputIt )
		{
			SPFIt.Set( SPFIt.Value()/max_SPF );
			float temp = m_delta*( m_Propagation*SPFIt.Value()*gradientIt.Value() );
			outputIt.Set( outputIt.Value() + temp );
			if ( outputIt.Value()<0.0 )
			{
				outputIt.Set( -1.0 );
			}
			else
			{
				outputIt.Set( 1.0 );
			}
		}

		//Gaussian Regular the level set map (OutputImage)..
        m_gaussianFilter->SetInput(outputImage);
		m_gaussianFilter->Update();
		outputImage = m_gaussianFilter->GetOutput();
		iteratorNum++;

	}

	
	OutputImageType::Pointer tempImage = OutputImageType::New();
	tempImage->SetRegions( this->GetInput()->GetLargestPossibleRegion() );
	tempImage->CopyInformation( this->GetInput() );
	tempImage->Allocate();
	
	RegionType    t_ROIregion(m_RegionOfInterest);


{//copy result from origin Image to tempImage
	//Declare of Feature Image Iterator
	ImageRegionConstIterator<OutputImageType >
		inputIt( this->GetInput(), this->GetInput()->GetLargestPossibleRegion() );
	//Declare of Feature Image Iterator
	ImageRegionIterator<OutputImageType >
		tempIt( tempImage, tempImage->GetLargestPossibleRegion() );

	for( inputIt.GoToBegin(), tempIt.GoToBegin() ; !inputIt.IsAtEnd(); ++inputIt, ++tempIt )
	{
		tempIt.Set( inputIt.Value() );
	}
	//此处加入各项异性调整
//	tempImage->DisconnectPipeline();
	if (m_ResampleThickSliceData)
	{
		CropFilterType::Pointer                 t_ROIFilter;
		t_ROIFilter = CropFilterType::New();

		t_ROIFilter->SetInput( tempImage );
		t_ROIFilter->SetRegionOfInterest( tempImage->GetLargestPossibleRegion() );
		



		// Compute the spacing after isotropic resampling.
		double minSpacing = NumericTraits< double >::max();
		for (int i = 0; i < ImageDimension; i++)
		{
			minSpacing = (minSpacing > tempImage->GetSpacing()[i] ?
				tempImage->GetSpacing()[i] : minSpacing);
		}

		// Try and reduce the anisotropy.
		SpacingType outputSpacing = tempImage->GetSpacing();
		for (int i = 0; i < ImageDimension; i++)
		{
			if (outputSpacing[i]/minSpacing > m_AnisotropyThreshold && m_ResampleThickSliceData)
			{
				outputSpacing[i] = minSpacing * m_AnisotropyThreshold;
			}
		}

		m_IsotropicResampler->SetInput( t_ROIFilter->GetOutput() );
		m_IsotropicResampler->SetOutputSpacing( outputSpacing );
		m_IsotropicResampler->GenerateOutputInformation();
		m_IsotropicResampler->Update();
		tempImage = m_IsotropicResampler->GetOutput();
		tempImage->SetRegions(tempImage->GetLargestPossibleRegion());
		ImageType::PointType t_point;
		this->GetInput()->TransformIndexToPhysicalPoint(m_RegionOfInterest.GetIndex(), t_point);
		ImageType::IndexType t_index;
		tempImage->TransformPhysicalPointToIndex(t_point, t_index);
		t_ROIregion.SetIndex(t_index);
		t_ROIregion.SetSize(outputImage->GetRequestedRegion().GetSize());
		 
 	}

}

{//Copy result from outputImage to tempImage
	//Declare of Output Image Iterator
	ImageRegionIterator< OutputImageType >
		outputIt( outputImage, outputImage->GetRequestedRegion() );

	ImageRegionIterator<OutputImageType >
		tempIt( tempImage, t_ROIregion );


	

	for( outputIt.GoToBegin(), tempIt.GoToBegin() ; !outputIt.IsAtEnd(); ++outputIt, ++tempIt )
	{
		if(outputIt.Value()>0)
			tempIt.Set( 5000 );
	}
}

  // Seeds

  // Graft the output.
  //  OutputImageType::Pointer tempImage = outputImage;

  tempImage->DisconnectPipeline();
  this->GraftOutput(tempImage);

   // DEBUGGING CODE
/*  typedef ImageFileWriter< OutputImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();
  writer->SetFileName("e:/output.mha");
  writer->SetInput(tempImage);
  writer->UseCompressionOn();
  writer->Write();*/
}


void SPFLevelSetImageFilter
::SetAbortGenerateData( bool abort )
{
  this->Superclass::SetAbortGenerateData(abort);
  this->m_CropFilter->SetAbortGenerateData(abort);
  this->m_IsotropicResampler->SetAbortGenerateData(abort);
}


void SPFLevelSetImageFilter
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}

}//end MySpace





#endif // SPFLEVELSETIMAGEFILTER_HXX_INCLUDED
