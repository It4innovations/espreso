#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "../configuration/globalconfiguration.h"

#include <QStandardItemModel>

static QStandardItem* addConfiguration(const espreso::Configuration &configuration, QStandardItem *item)
{
	for (size_t i = 0; i < configuration.orderedParameters.size(); i++) {
		item->appendRow(new QStandardItem(configuration.orderedParameters[i]->name.c_str()));
	}

	for (size_t i = 0; i < configuration.orderedSubconfiguration.size(); i++) {
		item->appendRow(addConfiguration(
				*configuration.orderedSubconfiguration[i],
				new QStandardItem(configuration.orderedSubconfiguration[i]->name.c_str())
		));
	}
	return item;
}

MainWindow::MainWindow(int *argc, char ***argv)
: QMainWindow(0), ui(new Ui::MainWindow)
{
	ui->setupUi(this);

	espreso::GlobalConfiguration configuration(argc, argv);
	QStandardItemModel *model = new QStandardItemModel();

	model->appendRow(addConfiguration(configuration, new QStandardItem("ESPRESO")));
	ui->treeView->setModel(model);
}

MainWindow::~MainWindow()
{
	delete ui;
}

void MainWindow::on_pushButton_clicked()
{
	ui->treeView->move(QPoint(20, 20));
}

void MainWindow::mousePressEvent(QMouseEvent *ev)
{
	last = ev->pos();
}

void MainWindow::mouseMoveEvent(QMouseEvent *ev)
{
	QPoint diff = last - ev->pos();
	ui->treeView->move(ui->treeView->pos() - diff);
	last = ev->pos();
}
