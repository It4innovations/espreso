#ifndef PLOT_H
#define PLOT_H

#include <QWidget>
#include <QGraphicsScene>

namespace Ui {
class Plot;
}

struct PlotLabel
{
    qreal val;
    QString text;
};

class Plot : public QWidget
{
    Q_OBJECT

public:
    explicit Plot(QWidget *parent = 0);
    ~Plot();

private:
    Ui::Plot *ui;

    QGraphicsScene* scene;

    QFont font;
    qreal fontSize;
    qreal fontSizeHalf;

    qreal fnXLeftBoundary;
    qreal fnXRightBoundary;
    qreal fnYTopBoundary;
    qreal fnYBottomBoundary;
    qreal fnXAxisLen;
    qreal fnYAxisLen;
    qreal rectFnXRatio;
    qreal rectFnYRatio;

    QGraphicsRectItem* rect;
    qreal rectX;
    qreal rectY;
    qreal rectWidth;
    qreal rectHeight;
    qreal plotXMiddle;
    qreal plotYMiddle;

    qreal computePrecision(qreal intervalLength);

    qreal fnXToRect(qreal x);
    qreal fnYToRect(qreal y);
    qreal rectXToFn(qreal x);
    qreal rectYToFn(qreal y);

    qreal fn(qreal x);

    QVector<PlotLabel> axisLabels(qreal leftBound, qreal rightBound, int limit);
    void generateAxisLabels(qreal leftBound, qreal rightBound, int limit,
                            QVector<PlotLabel>& labels, qreal offset);
    qreal maxLabelLength(QVector<PlotLabel>& labels);

    void drawPoint(QPointF p);
    void drawAxises();
    void drawFn();
    void drawXLabels(QVector<PlotLabel>& labels, int labelPointLength = 10);
    void drawYLabels(QVector<PlotLabel>& labels, int labelPointLength = 10);
};

#endif // PLOT_H
