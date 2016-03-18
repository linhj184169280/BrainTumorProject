#ifndef SPFLEVELSETIMAGEFILTER_H_INCLUDED
#define SPFLEVELSETIMAGEFILTER_H_INCLUDED


#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkFixedArray.h"
#include "itkCommand.h"
#include "itkImageSpatialObject.h"
#include "itkLandmarkSpatialObject.h"

#include "itkRegionOfInterestImageFilter.h"
#include "IsotropicResamplerImageFilter.h"



#include <string>
using namespace itk;

namespace MySpace
{
/** \class SPFLevelSetImageFilter
 * \ingroup ITKLesionSizingToolkit
 */
typedef float PixelType;
const static unsigned int ImageDimension = 3;
typedef itk::Image< PixelType, ImageDimension > ImageType;


class SPFLevelSetImageFilter  : public ImageToImageFilter<ImageType, ImageType>
{
public:
  /** Standard "Self" & Superclass typedef.  */
  typedef SPFLevelSetImageFilter                    Self;
  typedef ImageToImageFilter<ImageType, ImageType>     Superclass;

  /** Image typedef support   */
  typedef ImageType  InputImageType;
  typedef ImageType OutputImageType;
  typedef ImageType FeatureImageType;
  typedef ImageType InternalImageType;

  /** SmartPointer typedef support  */
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Define pixel types. */
  typedef ImageType::PixelType         InputImagePixelType;
  typedef ImageType::PixelType        OutputImagePixelType;
  typedef ImageType::IndexType         IndexType;
  typedef InputImageType::SpacingType    SpacingType;


  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Typedef to describe the output image region type. */
  typedef ImageType::RegionType OutputImageRegionType;
  typedef ImageType::RegionType RegionType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(SPFLevelSetImageFilter, ImageToImageFilter);


  /** ImageDimension constant    */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      ImageType::ImageDimension);

  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      ImageType::ImageDimension);

  typedef LandmarkSpatialObject< ImageDimension >     InputSpatialObjectType;

  virtual void GenerateInputRequestedRegion()
            throw(InvalidRequestedRegionError);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<InputImagePixelType>));
  itkConceptMacro(OutputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<OutputImagePixelType>));
  itkConceptMacro(SameDimensionCheck,
    (Concept::SameDimension<ImageDimension, OutputImageDimension>));
  itkConceptMacro(OutputIsFloatingPointCheck,
    (Concept::IsFloatingPoint<OutputImagePixelType>));
  /** End concept checking */
#endif

  /** Set the ROI */
  itkSetMacro( RegionOfInterest, RegionType );
  itkGetMacro( RegionOfInterest, RegionType );


  /** Turn On/Off isotropic resampling prior to running the segmentation */
  itkSetMacro( ResampleThickSliceData, bool );
  itkGetMacro( ResampleThickSliceData, bool );
  itkBooleanMacro( ResampleThickSliceData );

  /** If ResampleThickSliceData is ON, set the maximum anisotropy. Defauls
   * to twice the minimum spacing of the data. That will mean that
   * 0.7 x 0.7 x 1.25 mm CT data will not be resampled. However
   * 0.7 x 0.7 x 2.5 mm CT data will be resampled to 0.7 x 0.7 x 1.25 mm
   * thick slices. Values less than 1.0 result in supersampling of the
   * data. A value of 1 results in Isotropic resampling of the data.
   */
  itkSetMacro( AnisotropyThreshold, double );
  itkGetMacro( AnisotropyThreshold, double );

  /** Maximum number of iterations that the level set solve will run. */
  itkSetMacro( MaximumNumberOfIterations, unsigned int );
  itkGetMacro( MaximumNumberOfIterations, unsigned int );

  /** Propagation Speed. */
  itkSetMacro( Propagation, double );
  itkGetMacro( Propagation, double );

    /** Parameter that controls the updata : u = u + m_delta*(increase value) . */
  itkSetMacro( delta, double );
  itkGetMacro( delta, double );

    /** Gaussian sigma. */
  itkSetMacro( Sigma, double );
  itkGetMacro( Sigma, double );

  /** Gaussian NeiborSize */
  itkSetMacro( NeiborSize, unsigned int );
  itkGetMacro( NeiborSize, unsigned int );


  typedef itk::LandmarkSpatialObject< ImageDimension >    SeedSpatialObjectType;
  typedef SeedSpatialObjectType::PointListType   PointListType;

  void SetSeeds( PointListType p ) { this->m_Seeds = p; }
  PointListType GetSeeds() { return m_Seeds; }


  /** Override the superclass implementation so as to set the flag on all the
   * filters within our lesion segmentation pipeline */
  virtual void SetAbortGenerateData( const bool );

protected:
  SPFLevelSetImageFilter();
  SPFLevelSetImageFilter(const Self&) {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  virtual void GenerateOutputInformation();
  void GenerateData();


  //Definition of Filters
  typedef RegionOfInterestImageFilter< InputImageType, InputImageType > CropFilterType;
  typedef IsotropicResamplerImageFilter IsotropicResamplerType;


  //Output Image
  OutputImageType::Pointer outputImage;


private:
  virtual ~SPFLevelSetImageFilter(){};

  bool                                             m_ResampleThickSliceData;
  double                                           m_AnisotropyThreshold;
  RegionType                                       m_RegionOfInterest;
  SeedSpatialObjectType::PointListType    m_Seeds;

  CropFilterType::Pointer                 m_CropFilter;
  IsotropicResamplerType::Pointer         m_IsotropicResampler;

  double			m_Sigma;
  unsigned short    m_NeiborSize;
  double			m_delta;
  double			m_Propagation;
  unsigned int      m_MaximumNumberOfIterations;


};

}//end MySpace

#ifndef ITK_MANUAL_INSTANTIATION
//#include "SPFLevelSetImageFilter.hxx"
#endif


#endif // SPFLEVELSETIMAGEFILTER_H_INCLUDED
