
#include "gui/parallel/mpimanager.h"
#include "gui/workspace/workspacewindow.h"

#include "esinfo/mpiinfo.h"

#include <locale>
#include <QDebug>
#include <QApplication>

using namespace espreso;

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    QApplication a(argc, argv);

    MpiManager mpim(argc, argv);

    WorkspaceWindow ww(&mpim);

    setlocale(LC_ALL, "C");
    setlocale(LC_CTYPE, "C");

    if (info::mpi::rank == 0)
    {
        ww.init();
        ww.showMaximized();

        // CSS
//        QCoreApplication::setAttribute(Qt::AA_UseStyleSheetPropagationInWidgetStyles, true);
//        QFile css(":/stylesheets/stylesheets/general.css");
//        css.open(QFile::ReadOnly);
//        QTextStream stream(&css);
//        a.setStyleSheet(stream.readAll());
    }

    mpim.loop();

    if (info::mpi::rank == 0)
    {
        a.exec();
        mpim.masterExit();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
