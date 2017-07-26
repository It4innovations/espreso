#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QMouseEvent>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	explicit MainWindow(int *argc, char ***argv);
	~MainWindow();

private slots:
	void mouseMoveEvent(QMouseEvent *ev);
	void mousePressEvent(QMouseEvent *ev);

	void on_pushButton_clicked();

private:
	Ui::MainWindow *ui;
	QPoint last;
};

#endif // MAINWINDOW_H
