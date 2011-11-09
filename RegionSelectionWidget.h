#ifndef REGIONSELECTION_WIDGET
#define REGIONSELECTION_WIDGET

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

//QT Includes
#include <QtGui>
#include <QtGui/QMainWindow>
#include <QObject>
#include <QtGui/QMenuBar>
#include <QtGui/QMenu>
#include <QtGui/QLabel>
#include <QtCore/QSignalMapper>

//ITK Includes to subsample
#include "itkImage.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageRegionConstIterator.h"

//FTK
#include <ftkImage/ftkImage.h>

class SelectedRegion;

class RegionSelectWidget : public QMainWindow
{
	Q_OBJECT;

	typedef itk::Image< unsigned char, 3 >   IOType;
public:
	RegionSelectWidget( ftk::Image::Pointer Input_channels, SelectedRegion *in_region, QWidget *parent=0, Qt::WindowFlags flags=0 );
	~RegionSelectWidget();
private:
	size_t originalHt, originalWd, rescaledWd, rescaledHt;
	double factor;
	const ftk::Image::Info *info;
	void generateBaseImages(void);
	void setUpDisplayImage(void);

protected:
	QImage *input_2d; //Rescaled input image
	QImage displayImage;
	QLabel *imageLabel;
	QScrollArea *scrollArea;
	ftk::Image::Pointer input_im; //Input image pointer
	std::vector<IOType::Pointer> displayImageITK;
	SelectedRegion *region;
	void paintEvent(QPaintEvent *event);

protected slots:
	void regionChange(void);

};

class SelectedRegion : public QObject
{
	Q_OBJECT;
public:
	SelectedRegion();
	~SelectedRegion();
	typedef struct{
		size_t StartX, StartY, Width, Height;
		double ZoomScale;
	} RegionInfo;

	void       SetRegion( RegionInfo& );
	RegionInfo GetRegion(){ RegionInfo cpy = region; return cpy; }

signals:
	void change();

private:
	RegionInfo region;

};
#endif
