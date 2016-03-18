#ifndef __LesionSegmentationCLI_h
#define __LesionSegmentationCLI_h

#include "metaCommand.h"
#include "itkFixedArray.h"
#include "itkLandmarkSpatialObject.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include <fstream>
#include <sys/types.h>
#if !defined(_MSC_VER)
  #include <dirent.h> // exists only on POSIX type compilers
#else
  #include "dirent_win.h" // exists only on POSIX type compilers
#endif
#include <errno.h>
#include <vector>
#include <string>
#include <iostream>

class LesionSegmentationCLI : public MetaCommand
{
public:
  typedef float PixelType;
  const static unsigned int ImageDimension = 3;
  typedef itk::Image< PixelType, ImageDimension > ImageType;

  typedef itk::LandmarkSpatialObject< 3 >    SeedSpatialObjectType;
  typedef SeedSpatialObjectType::PointListType   PointListType;

  LesionSegmentationCLI( int argc, char *argv[] ) : MetaCommand()
  {
    m_Image = NULL;
    this->DisableDeprecatedWarnings();

    this->AddArgument("InputImage",false,"Input image to be segmented.");
    this->AddArgument("InputDICOMDir",false,"DICOM directory containing series of the Input image to be segmented.");
	this->AddArgument("OutputDir",true,"The File Directory to contain the results.");
    this->AddArgument("OutputImage", false, "Output segmented image name");
  //  this->AddArgument("OutputMesh", false, "Output segmented surface (STL filename expected)");
    this->AddArgument("OutputROI", false, "Write the ROI within which the segmentation will be confined to (for debugging purposes)");
 //   this->AddArgument("Visualize", false, "Visualize the input image and the segmented surface.", MetaCommand::BOOL, "0");
  //  this->AddArgument("Wireframe", false, "Visualize the input image and the segmented surface as a wireframe. Only valid if the Visualize flag is also enabled.", MetaCommand::BOOL, "0");
    this->AddArgument("IgnoreDirection", false, "Ignore the direction of the DICOM image", MetaCommand::BOOL, "0");
 //    this->AddArgument("PartSolid", false, "Default is to assume parameters for a solid lesion. Specify this if the lesion is part-solid.", MetaCommand::BOOL, "0");
    this->AddArgument("ROI", false,
                    "Bounds of the ROI if any, 6 parameters", MetaCommand::LIST);
  //  this->AddArgument("Sigma", false,
  //    "Manually specify sigma. This is an array with 3 values in physical units. This defaults to the maximum spacing in the dataset, if unspecified",
 //     MetaCommand::LIST);
    this->AddArgument("Seeds", false,
      "Manually specify seeds in physical coordinates. At least one seed must be specified using for a segmentation to be generated. Usage is of the form --Seeds 3 X1 Y1 Z1 (for 1 seed) or --Seeds 6 X1 Y1 Z1 X2 Y2 Z2 (for 2 seeds) etc..",
      MetaCommand::LIST);
    this->AddArgument("SeedUnitsInPixels", false,
      "Are the seeds specified in pixel coordinates ? (Note that pixel coords start at 0 index). If so, use this flag. By default seeds are assumed to be in physical coordinates.", MetaCommand::BOOL, "0");
    this->AddArgument("MaximumRadius", false, "Maximum radius of the lesion in mm. This can be used as alternate way of specifying the bounds. You specify a seed and a value of say 20mm, if you know the lesion is smaller than 20mm..", MetaCommand::FLOAT, "30");
  //  this->AddArgument("Screenshot",false,"Screenshot PNG file of the final segmented surface (requires \"Visualize\" to be ON.");
  //  this->AddArgument("ShowBoundingBox", false,
  //    "Show the ROI used for the segmentation as a bounding box on the visualization.", MetaCommand::BOOL, "0");
  //  this->AddArgument("GetZSpacingFromSliceNameRegex",false,
  //    "This option was added for the NIST Biochange challenge where the Z seed index was specified by providing the filename of the DICOM slice where the seed resides. Hence if this option is specified, the Z value of the seed is ignored.");
    this->AddArgument("ResampleThickSliceData", false,
                    "If we want to resample for isotropic, set it to 1, else set it to 0.(Default:1) ", MetaCommand::BOOL,"1");

    if(!this->Parse(argc, argv))
      {
      // We can't invoke errors from constructors..
      exit(-1);
      }
    }


  double* GetROI();
  PointListType GetSeeds();

  ImageType::Pointer GetImage( std::string dir, bool ignoreDirection );

  void SetImage( ImageType * image );


  typedef ImageType::IndexType IndexType;
  typedef IndexType::IndexValueType IndexValueType;

protected:

  void AddArgument( std::string name,
                    bool required,
                    std::string description,
                    TypeEnumType type = MetaCommand::STRING,
                    std::string defVal = "" );



  // get files in dir
  std::vector< std::string > GetFilesInDirectory(std::string dir);



 /* double GetZPositionFromFile( std::string file, std::string &sopInstanceUID )
    {
    typedef itk::ImageFileReader< ImageType > ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( file );
    reader->Update();

    itk::MetaDataDictionary dict = reader->GetMetaDataDictionary();
    itk::ExposeMetaData<std::string>( dict, "0008|0018", sopInstanceUID );

    return reader->GetOutput()->GetOrigin()[2];
    }*/



  double ROI[6];
  ImageType * m_Image;
};

LesionSegmentationCLI::ImageType::Pointer GetImage( std::string dir, bool ignoreDirection = false );

//#include "LesionSegmentationCLI.hxx"

#endif
