
#include "tabletypewidget.h"

#include "gui/validators/validatordelegate.h"
#include "gui/validators/validatorfactory.h"

#include <QRegExp>
#include <QString>
#include <QStringList>
#include <QDebug>

using namespace espreso;

QStringList TableTypeWidget::headlines()
{
    QStringList result;
    result << QObject::tr("x") << QObject::tr("f(x)");

    return result;
}

TableTypeWidget::TableTypeWidget(QWidget* parent) :
    TableWidget(2, TableTypeWidget::headlines(), parent)
{
    this->mTable->setItemDelegateForColumn(0, new ValidatorDelegate(new DoubleValidatorFactory()));
    this->mTable->setItemDelegateForColumn(1, new ValidatorDelegate(new DoubleValidatorFactory()));
}

void TableTypeWidget::addData(const QString& data)
{
    QRegExp prefix("TABULAR\\[");
    QStringList with_right_bracket = data.split(prefix, QString::SkipEmptyParts);
    QRegExp suffix("\\]");
    QStringList raw = with_right_bracket[0].split(suffix, QString::SkipEmptyParts);
	if (raw.size() == 0) return;
    QRegExp delimiter("\\;");
    QStringList pairs = raw[0].split(delimiter, QString::SkipEmptyParts);

    foreach (QString p, pairs) {
		QStringList coordinates = p.trimmed().split(QRegExp("\\,"), QString::SkipEmptyParts);
        this->addRow(coordinates.toVector());
    }
}

QString TableTypeWidget::data()
{
    QLatin1String start("TABULAR[");
    QLatin1String end("]");

    QStringList pairs;
    for (int i = 0; i < mModel->rowCount() - 1; ++i)
    {
        QModelIndex indexX = mModel->index(i, 0);
        QString _x = mModel->data(indexX).toString();
        QModelIndex indexY = mModel->index(i, 1);
        QString _y = mModel->data(indexY).toString();
        QString line = _x + QLatin1String(",") + _y + QLatin1String(";");
        pairs << line;
    }

    return start + pairs.join(QLatin1String(" ")) + end;
}
