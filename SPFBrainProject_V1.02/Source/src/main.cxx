#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMImageIOFactory.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkMetaImageIOFactory.h"
#include "itkImageSeriesReader.h"
#include "itkEventObject.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkOrientImageFilter.h"
#include "itkLandmarkSpatialObject.h"
#include "itkImageToVTKImageFilter.h"
#include "vtkImageData.h"
#include <fstream>
#include "time.h"

#include "SPFLevelSetImageFilter.h"



using namespace std;

// This needs to come after the other includes to prevent the global definitions
// of PixelType to be shadowed by other declarations.
//#include "SPFLesionSegmentationFilter.h"
#include "LesionSegmentationCLI.h"

#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
/*	LesionSegmentationCLI args( argc, argv );
	///////////////////////////////////////////////////////////////
	string OutputDirPath;
	string command;
	cout<<args.GetValueAsString("OutputDir")<<endl;
	if (!args.GetValueAsString("OutputDir").empty())
	{
		OutputDirPath = args.GetValueAsString("OutputDir");

	}
	//mkdir(OutputDirPath.c_str());
	command = "mkdir \""+OutputDirPath+"\"";
	system(command.c_str());
//	ofstream tcout( OutputDirPath+"/record.txt" );
//	ofstream tcout( "e:/record.txt" );
	///////////////////////////////////////////////////////////////
  //register DICOM and META IO factories
  itk::ObjectFactoryBase::RegisterFactory(itk::GDCMImageIOFactory::New());
  itk::ObjectFactoryBase::RegisterFactory(itk::MetaImageIOFactory::New());



  typedef LesionSegmentationCLI::ImageType ImageType;
  const unsigned int ImageDimension = LesionSegmentationCLI::ImageDimension;

  typedef itk::ImageFileReader< ImageType > InputReaderType;
  typedef itk::ImageFileWriter< ImageType > OutputWriterType;
  typedef MySpace::SPFLevelSetImageFilter SegmentationFilterType;
  */

  // Read the volume
/*  InputReaderType::Pointer reader = InputReaderType::New();
  ImageType::Pointer image;

  std::cout << "Reading " << args.GetValueAsString("InputImage") << ".." << std::endl;
  if (!args.GetValueAsString("InputDICOMDir").empty())
    {
    std::cout << "Reading from DICOM dir " << args.GetValueAsString("InputDICOMDir") << ".." << std::endl;
    image = GetImage(
      args.GetValueAsString("InputDICOMDir"),
      args.GetValueAsBool("IgnoreDirection"));

    if (!image)
      {
      std::cerr << "Failed to read the input image" << std::endl;
      return EXIT_FAILURE;
      }
    }

  if (!args.GetValueAsString("InputImage").empty())
    {
    reader->SetFileName(args.GetValueAsString("InputImage"));
    reader->Update();
    image = reader->GetOutput();
    }
	*/

  //To make sure the tumor polydata aligns with the image volume during
  //vtk rendering in ViewImageAndSegmentationSurface(),
  //reorient image so that the direction matrix is an identity matrix.
/*  typedef itk::ImageFileReader< ImageType > InputReaderType;
  itk::OrientImageFilter<ImageType,ImageType>::Pointer orienter =
  itk::OrientImageFilter<ImageType,ImageType>::New();
  orienter->UseImageDirectionOn();
  ImageType::DirectionType direction;
  direction.SetIdentity();
  orienter->SetDesiredCoordinateDirection (direction);
  orienter->SetInput(image);
  orienter->Update();
  image = orienter->GetOutput();

  // Set the image object on the args
  args.SetImage( image );*/

  // Compute the ROI region

/*  double *roi = args.GetROI();
  ImageType::RegionType region = image->GetLargestPossibleRegion();

  // convert bounds into region indices
  ImageType::PointType p1, p2;
  ImageType::IndexType pi1, pi2;
  ImageType::IndexType startIndex;
  for (unsigned int i = 0; i < ImageDimension; i++)
    {
    p1[i] = roi[2*i];
    p2[i] = roi[2*i+1];
    }

  image->TransformPhysicalPointToIndex(p1, pi1);
  image->TransformPhysicalPointToIndex(p2, pi2);

  ImageType::SizeType roiSize;
  for (unsigned int i = 0; i < ImageDimension; i++)
    {
    roiSize[i] = fabs((double)(pi2[i] - pi1[i]));
    startIndex[i] = (pi1[i]<pi2[i])?pi1[i]:pi2[i];
    }
  ImageType::RegionType roiRegion( startIndex, roiSize );
  std::cout << "ROI region is " << roiRegion << std::endl;
  if (!roiRegion.Crop(image->GetBufferedRegion()))
    {
    std::cerr << "ROI region has no overlap with the image region of"
              << image->GetBufferedRegion() << std::endl;
    return EXIT_FAILURE;
    }*/

/*  std::cout << "ROI region is " << roiRegion <<  " : covers voxels = "
    << roiRegion.GetNumberOfPixels() << " : " <<
    image->GetSpacing()[0]*image->GetSpacing()[1]*image->GetSpacing()[2]*roiRegion.GetNumberOfPixels()
   << " mm^3" << std::endl;*/
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
/*  // Write ROI if requested
  if (args.GetOptionWasSet("OutputROI"))
    {
    typedef itk::RegionOfInterestImageFilter<
      ImageType, ImageType > ROIFilterType;


    ROIFilterType::Pointer roiFilter = ROIFilterType::New();
    roiFilter->SetRegionOfInterest( roiRegion );

    typedef itk::ImageFileWriter< ImageType > ROIWriterType;
    ROIWriterType::Pointer roiWriter = ROIWriterType::New();
    roiWriter->SetFileName( OutputDirPath +"/"+ args.GetValueAsString("OutputROI") );
    roiFilter->SetInput( image );
    roiWriter->SetInput( roiFilter->GetOutput() );

    std::cout << "Writing the ROI image as: " <<
      args.GetValueAsString("OutputROI") << std::endl;
    try
      {
      roiWriter->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "ExceptionObject caught !" << err << std::endl;
      return EXIT_FAILURE;
      }
    }
*/

  // Progress reporting
//  typedef itk::LesionSegmentationCommandLineProgressReporter ProgressReporterType;
//  ProgressReporterType::Pointer progressCommand =
//    ProgressReporterType::New();


  //Parameters:
  int seedx = 200;
  int seedy = 304;
  int seedz = 11;
  double radius = 30; //physical radius
  string inputDirPath = "../Data/tumor";
  string outputPath = "../Data/tumor_SPF";
  //End of Parameters

  const static unsigned int ImageDimension = 3;
  typedef itk::Image< float, ImageDimension > ImageType;
  typedef itk::ImageFileWriter< ImageType > OutputWriterType;
  typedef MySpace::SPFLevelSetImageFilter SegmentationFilterType;

  ImageType::Pointer image;

  std::cout << "Reading from DICOM dir " <<inputDirPath << ".." << std::endl;
  image = GetImage(inputDirPath);


  typedef itk::LandmarkSpatialObject< 3 >    SeedSpatialObjectType;
  typedef SeedSpatialObjectType::PointListType   PointListType;

  typedef ImageType::IndexValueType  IndexValueType;
  typedef ImageType::IndexType		 IndexType;
  //seed point 
  IndexType t_index = {{
	 static_cast< IndexValueType >(seedx),
	 static_cast< IndexValueType >(seedy),
	 static_cast< IndexValueType >(seedz) }};
  ImageType::PointType point;
  image->TransformIndexToPhysicalPoint(t_index, point);
  double sx = point[0];
  double sy = point[1];
  double sz = point[2];
  PointListType seeds(1);
  seeds[0].SetPosition(sx,sy,sz);			//transform to physical seed

  //calculate ROI
  ImageType::PointType p1, p2;
  ImageType::IndexType pi1, pi2;
  ImageType::IndexType startIndex;
  for (int i = 0; i < ImageDimension; i++)
  {
	  p1[i] = seeds[0].GetPosition()[i] - radius;
	  p2[i] = seeds[0].GetPosition()[i] + radius;
  }
  image->TransformPhysicalPointToIndex(p1, pi1);
  image->TransformPhysicalPointToIndex(p2, pi2);
  ImageType::SizeType roiSize;
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
	  roiSize[i] = fabs((double)(pi2[i] - pi1[i]));
	  startIndex[i] = (pi1[i]<pi2[i])?pi1[i]:pi2[i];
  }
  ImageType::RegionType roiRegion( startIndex, roiSize );

  //Running segmentation
  std::cout << "\n Running the segmentation filter." << std::endl;
  SegmentationFilterType::Pointer seg = SegmentationFilterType::New();
  seg->SetInput(image);
  seg->SetSeeds(seeds);
//  seg->SetSeeds(args.GetSeeds());
  seg->SetRegionOfInterest(roiRegion);
  seg->SetResampleThickSliceData(false);		//Interpolation
//  seg->SetResampleThickSliceData(args.GetValueAsBool("ResampleThickSliceData"));

  clock_t start, finish;
  double  duration;
  start = clock();
try
{
	seg->Update();
}
catch (itk::ExceptionObject e)
{
	cout<<e<<endl;
}


  finish = clock();
  duration = (double)(finish - start);

  cout<<"Time = "<<duration/1000<<"sec, = "<<duration<<"ms"<<endl;

  //transform to vtkImage
  typedef itk::ImageToVTKImageFilter< ImageType > RealITKToVTKFilterType;
  RealITKToVTKFilterType::Pointer itk2vtko = RealITKToVTKFilterType::New();
  itk2vtko->SetInput( seg->GetOutput() );
  itk2vtko->Update();

  

  //output the result
  OutputWriterType::Pointer writer = OutputWriterType::New();
  writer->SetFileName(outputPath+"/result.mha");
  writer->SetInput(seg->GetOutput());
  writer->Update();

/*  if (!args.GetValueAsString("OutputImage").empty())
    {
    std::cout << "Writing the output segmented level set."
      << args.GetValueAsString("OutputImage") <<std::endl;
    OutputWriterType::Pointer writer = OutputWriterType::New();
    writer->SetFileName(OutputDirPath +"/" + args.GetValueAsString("OutputImage"));
    writer->SetInput(seg->GetOutput());
    writer->Update();
    }
*/


  return EXIT_SUCCESS;
}
