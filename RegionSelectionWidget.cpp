#include "RegionSelectionWidget.h"

RegionSelectWidget::RegionSelectWidget( ftk::Image::Pointer Input_channels, SelectedRegion *in_region,
						QWidget *parent, Qt::WindowFlags flags ): QMainWindow(parent,flags)
{
	//Set input pointer and get info
	region = in_region;
	input_im = Input_channels;
	info = input_im->GetImageInfo();
	originalHt = info->numRows;
	originalWd = info->numColumns;

	if( originalHt > originalWd ){
		rescaledHt = 400;
		rescaledWd = 400*(originalWd/originalHt);
		factor     = originalHt/400;
	} else {
		rescaledHt = 400*(originalHt/originalWd);
		rescaledWd = 400;
		factor     = originalWd/400;
	}

	generateBaseImages();

	connect(region, SIGNAL(change()), this, SLOT(regionChange()));

	setUpDisplayImage();

	imageLabel = new QLabel();
	imageLabel->setBackgroundRole(QPalette::Base);
	imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	imageLabel->setMouseTracking(true);
	imageLabel->setScaledContents(true);
	imageLabel->resize(0,0);

	scrollArea = new QScrollArea;
	scrollArea->setBackgroundRole(QPalette::Dark);
	scrollArea->setWidget(imageLabel);
	setCentralWidget(scrollArea);

//	this->createActions();
//	this->createMenus();

	setWindowTitle(tr("Slide Viewer"));
	this->resize(500,500);
}

RegionSelectWidget::~RegionSelectWidget(){
}

void RegionSelectWidget::setUpDisplayImage(void){
	input_2d = new QImage(rescaledWd, rescaledHt, QImage::Format_ARGB32_Premultiplied);
	input_2d->fill(qRgb(0,0,0));
	QPainter painter(input_2d);

	painter.setCompositionMode(QPainter::CompositionMode_Plus);

	for (int i=0; i < (*info).numChannels; i++){
		QImage gray(rescaledWd, rescaledHt, QImage::Format_ARGB32_Premultiplied);
		std::vector<unsigned char> color = (*info).channelColors[i];
		gray.fill(qRgb(color[0],color[1],color[2]));

		//Copy channel into a continuous mem array
		unsigned char *mem_p;
		mem_p = (unsigned char*) malloc( rescaledWd*rescaledHt*sizeof(unsigned char) ) ;
		typedef itk::ImageRegionConstIterator< IOType > ConstIteratorType;
		ConstIteratorType pix_buf( displayImageITK.at(i), displayImageITK.at(i)->GetRequestedRegion() );
		pix_buf.GoToBegin();
		size_t ind = 0;
		for ( pix_buf.GoToBegin(); !pix_buf.IsAtEnd(); ++pix_buf, ++ind )
			mem_p[ind] = (pix_buf.Get());
		QImage img(mem_p,rescaledWd,rescaledHt,QImage::Format_Indexed8);
		gray.setAlphaChannel(img);
		painter.drawImage(0,0,gray);
	}
	this->repaint();
}

void RegionSelectWidget::generateBaseImages(void){
//Adapted from code in ITK/Examples/Filtering/SubsampleVolume.cxx

	//GetItk pointer returns image with dim=3, it is really a 2-d image
	typedef itk::Image< float, 3 >   PrcType; //Processing type

	for (int i=0; i <((*info).numChannels); i++){
		IOType::Pointer inputImage =
				input_im->GetItkPtr<unsigned char>(0, i);

		typedef itk::CastImageFilter< IOType, PrcType > CastFilterType;
		CastFilterType::Pointer caster = CastFilterType::New();
		caster->SetInput( inputImage );

		typedef itk::RecursiveGaussianImageFilter< PrcType, PrcType
								 > GaussianFilterType;
		GaussianFilterType::Pointer smoother = GaussianFilterType::New();
		smoother->SetInput( caster->GetOutput() );
		smoother->SetNormalizeAcrossScale( false );
		const IOType::SpacingType& inputSpacing = inputImage->GetSpacing();
		const double sigma = inputSpacing[0] * factor;
		smoother->SetSigma( sigma );

		typedef itk::ResampleImageFilter< PrcType, IOType >
								ResampleFilterType;
		ResampleFilterType::Pointer resampler = ResampleFilterType::New();

		typedef itk::IdentityTransform< double, 3 >  TransformType;
		TransformType::Pointer transform = TransformType::New();
		transform->SetIdentity();
		resampler->SetTransform( transform );

		typedef itk::LinearInterpolateImageFunction< PrcType, double >
								InterpolatorType;
		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		resampler->SetInterpolator( interpolator );
		resampler->SetDefaultPixelValue( 0 );

		IOType::SpacingType spacing;
		spacing[0] = inputSpacing[0] * factor;
		spacing[1] = inputSpacing[1] * factor;
		spacing[2] = inputSpacing[2];

		resampler->SetOutputSpacing( spacing );
		resampler->SetOutputOrigin( inputImage->GetOrigin() );
		resampler->SetOutputDirection( inputImage->GetDirection() );

		IOType::SizeType inputSize =
				inputImage->GetLargestPossibleRegion().GetSize();
		typedef IOType::SizeType::SizeValueType SizeValueType;
		IOType::SizeType size;
		size[0] = static_cast< SizeValueType >( inputSize[0] / factor );
		size[1] = static_cast< SizeValueType >( inputSize[1] / factor );
		size[2] = static_cast< SizeValueType >( inputSize[2] );

		resampler->SetSize( size );
		resampler->SetInput( smoother->GetOutput() );
		IOType::Pointer displayImage = resampler->GetOutput();

		//Testing
		typedef itk::ImageFileWriter< IOType >  WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetInput( resampler->GetOutput() );
		std::stringstream ss; ss<<i;
		std::string out_file; out_file = "/home/gramak/downsample" + ss.str() + ".tif";
		writer->SetFileName( out_file );

		try{
			resampler->Update();
			writer->Update();  //Testing
		}
		catch( itk::ExceptionObject & excep ){
			std::cerr << "Exception caught!" << std::endl;
			std::cerr << excep << std::endl;
		}
		displayImageITK.push_back( displayImage );
	}
}

void RegionSelectWidget::paintEvent(QPaintEvent * event){
	QWidget::paintEvent(event);

	if(input_2d->height() <= 0 || input_2d->width() <= 0)
		return;

	displayImage = *input_2d;

	QPainter painter(&displayImage);
	painter.setRenderHint(QPainter::Antialiasing);
	painter.setPen(Qt::red);

	RegionInfo regInf = region->SelectedRegion();
	painter.drawRect ( region->StartX, region->StartY, region->Width, region->Height );


}

void RegionSelectWidget::regionChange(void){
	
}

SelectedRegion::SelectedRegion(){
	region.StartX=0; region.StartY=0;
	region.Width =0; region.Height=0;
	region.ZoomScale=0;
}

SelectedRegion::~SelectedRegion(){
}


void SelectedRegion::SetRegion( RegionInfo &input_region ){
	region = input_region;
	emit change();
t

