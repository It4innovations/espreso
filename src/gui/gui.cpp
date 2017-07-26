#include "mainwindow.h"
#include <QApplication>

#include "../configuration/environment.h"

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	QApplication a(argc, argv);
	MainWindow w(&argc, &argv);
	w.show();

	a.exec();

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}
