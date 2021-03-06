/**
 * \file
 * \author Karsten Rink
 * \date   2011-02-07
 * \brief  Implementation of the GeoTabWidget class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GEOTABWIDGET_H
#define GEOTABWIDGET_H

// ** INCLUDES **
#include "ui_GeoTabWidgetBase.h"

/**
 * \brief Widget containing GeoTreeView-objects.
 */
class GeoTabWidget : public QWidget, public Ui_GeoTabWidgetBase
{
	Q_OBJECT

public:
	GeoTabWidget(QWidget* parent = 0);

private:

private slots:
	void enableSaveButton(bool enable) { this->saveGeoPushButton->setEnabled(enable); };
	void enableRemoveButton(bool enable) { this->removeGeoPushButton->setEnabled(enable); };

};

#endif // GEOTABWIDGET_H
