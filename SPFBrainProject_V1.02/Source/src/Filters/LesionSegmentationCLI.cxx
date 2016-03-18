#ifndef LESIONSEGMENTATIONCLI_HXX_INCLUDED
#define LESIONSEGMENTATIONCLI_HXX_INCLUDED

#include "LesionSegmentationCLI.h"
#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"

LesionSegmentationCLI::ImageType::Pointer GetImage( std::string dir, bool ignoreDirection )
{
  const unsigned int Dimension = LesionSegmentationCLI::ImageDimension;
  typedef itk::Image< LesionSegmentationCLI::PixelType, Dimension >         ImageType;

  typedef itk::ImageSeriesReader< ImageType >        ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

  typedef itk::GDCMImageIO       ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();

  reader->SetImageIO( dicomIO );

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails( true );
  nameGenerator->AddSeriesRestriction("0008|0021" );
  nameGenerator->SetDirectory( dir );

  try
    {
    std::cout << std::endl << "The directory: " << std::endl;
    std::cout << std::endl << dir << std::endl << std::endl;
    std::cout << "Contains the following DICOM Series: ";
    std::cout << std::endl << std::endl;

    typedef std::vector< std::string >    SeriesIdContainer;

    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();

    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
    while( seriesItr != seriesEnd )
      {
      std::cout << seriesItr->c_str() << std::endl;
      seriesItr++;
      }


    std::string seriesIdentifier;
    seriesIdentifier = seriesUID.begin()->c_str();


    std::cout << std::endl << std::endl;
    std::cout << "Now reading series: " << std::endl << std::endl;
    std::cout << seriesIdentifier << std::endl;
    std::cout << std::endl << std::endl;


    typedef std::vector< std::string >   FileNamesContainer;
    FileNamesContainer fileNames;

    fileNames = nameGenerator->GetFileNames( seriesIdentifier );

    FileNamesContainer::const_iterator  fitr = fileNames.begin();
    FileNamesContainer::const_iterator  fend = fileNames.end();

    while( fitr != fend )
      {
      std::cout << *fitr << std::endl;
      ++fitr;
      }


    reader->SetFileNames( fileNames );

    try
      {
      reader->Update();
      }
    catch (itk::ExceptionObject &ex)
      {
      std::cout << ex << std::endl;
      return NULL;
      }


    ImageType::Pointer image = reader->GetOutput();
    ImageType::DirectionType direction;
    direction.SetIdentity();
    image->DisconnectPipeline();
    std::cout << "Image Direction:" << image->GetDirection() << std::endl;


    if (ignoreDirection)
      {
      std::cout << "Ignoring the direction of the DICOM image and using identity." << std::endl;
      image->SetDirection(direction);
      }
    return image;
    }
  catch (itk::ExceptionObject &ex)
    {
    std::cout << ex << std::endl;
    return NULL;
    }

  return NULL;
}

double* LesionSegmentationCLI::GetROI()
{
    if (this->GetOptionWasSet("ROI"))
      {
     // Default to be physical units
     //TO DO: deal with ROI input in pixel units
      std::list< std::string > bounds = this->GetValueAsList("ROI");
      std::list< std::string >::const_iterator fit = bounds.begin();
      for (unsigned int i = 0; fit != bounds.end(); ++fit, ++i)
        {
        this->ROI[i] = (float)atof((*fit).c_str());
        }
      }
    else
      {
      PointListType seeds = this->GetSeeds();
      seeds[0];
      for (int i = 0; i < 3; i++)
        {
        this->ROI[2*i] = seeds[0].GetPosition()[i] - this->GetValueAsFloat("MaximumRadius");
        this->ROI[2*i+1] = seeds[0].GetPosition()[i] + this->GetValueAsFloat("MaximumRadius");
        }
      }
    return this->ROI;
}

LesionSegmentationCLI::PointListType LesionSegmentationCLI::GetSeeds()
{
    std::list< std::string > seedsString = this->GetValueAsList("Seeds");
    std::list< std::string >::const_iterator fit = seedsString.begin();
    const unsigned int nb_of_markers = seedsString.size() / 3;
    PointListType seeds(nb_of_markers);
    for (unsigned int i = 0; i < nb_of_markers; i++)
      {
      double sx = (double)atof((*fit).c_str());
      ++fit;
      double sy = (double)atof((*fit).c_str());
      ++fit;
      double sz = (double)atof((*fit).c_str());
      ++fit;

      if (this->GetOptionWasSet("SeedUnitsInPixels"))
        {

        // Convert seeds from pixel units to physical units
        IndexType index = {{
          static_cast< IndexValueType >(sx),
          static_cast< IndexValueType >(sy),
          static_cast< IndexValueType >(sz) }};

        ImageType::PointType point;
        m_Image->TransformIndexToPhysicalPoint(index, point);
        sx = point[0];
        sy = point[1];
        sz = point[2];

        // Get the z spacing from the slice name regex..
        if (this->GetOptionWasSet("GetZSpacingFromSliceNameRegex"))
          {
          std::string substring =
            this->GetValueAsString("GetZSpacingFromSliceNameRegex");
          std::vector< std::string > filesInDir =
            this->GetFilesInDirectory(this->GetValueAsString("InputDICOMDir"));
          bool found = false;
          for (std::vector< std::string >::iterator it = filesInDir.begin();
              it != filesInDir.end(); ++it)
            {
            if (it->find(substring) != std::string::npos &&
                it->find("vvi") == std::string::npos)
              {
              std::string file = this->GetValueAsString("InputDICOMDir");
              file += "/";
              file += (*it);
              }
            }
          if (!found)
            {
            // Loop now and check the SOP instance UID. Some datasets in the
            // biochange challenge rely on the filename, yet others rely on
            // the SOP instance UID present in the file.

            for (std::vector< std::string >::iterator it = filesInDir.begin();
                it != filesInDir.end(); ++it)
              {
              if (it->find("vvi") == std::string::npos)
                {
                std::string file = this->GetValueAsString("InputDICOMDir");
                file += "/";
                file += (*it);
                }
              }

            if (!found)
              {
              std::cerr << "Could not find a file with matching SOP "
                        << substring << std::endl;
              exit(-1);
              }
            }
          }
        }

      // Sanity check
      std::cout << "Seed position in physical units: (" << sx << ","
                << sy << "," << sz << ")" << std::endl;
      ImageType::PointType pointSeed;
      pointSeed[0] = sx;
      pointSeed[1] = sy;
      pointSeed[2] = sz;
      IndexType indexSeed;
      m_Image->TransformPhysicalPointToIndex(pointSeed, indexSeed);
      if (!this->m_Image->GetBufferedRegion().IsInside(indexSeed))
        {
        std::cerr << "Seed with pixel units of index: " <<
            indexSeed << " does not lie within the image. The images extents are"
          << this->m_Image->GetBufferedRegion() << std::endl;
        exit(-1);
        }

      seeds[i].SetPosition(sx,sy,sz);
      }
    return seeds;
}

void LesionSegmentationCLI::SetImage( ImageType * image )
{
    this->m_Image = image;
}


void LesionSegmentationCLI::AddArgument( std::string name,
                    bool required,
                    std::string description,
                    TypeEnumType type ,
                    std::string defVal )
{
    this->SetOption(name,name,required,description);
    this->SetOptionLongTag(name,name);
    this->AddOptionField(name,name,type,true, defVal);
}



  std::vector< std::string > LesionSegmentationCLI::GetFilesInDirectory(std::string dir)
    {
    std::vector< std::string > files;
    DIR *dp;
    struct dirent *dirp;
    if((dp = opendir(dir.c_str())) == NULL)
      {
      std::cerr << "Error(" << errno << ") opening " << dir << std::endl;
      return files;
      }

    while ((dirp = readdir(dp)) != NULL)
      {
      files.push_back(std::string(dirp->d_name));
      }

    closedir(dp);
    return files;
    }

#endif // LESIONSEGMENTATIONCLI_HXX_INCLUDED
