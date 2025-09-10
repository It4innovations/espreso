
#ifndef TEXTITEMWIDGETFACTORY_H
#define TEXTITEMWIDGETFACTORY_H

#include "textitemwidget.h"
#include "filepathwidget.h"
#include "optionwidget.h"
#include "textwidget.h"
#include "ivalidatableobject.h"

#include "gui/validators/validatorfactory.h"
#include "gui/declarations/datatypeeditwidget.h"

#include <QWidget>

namespace espreso
{

struct ECFParameter;

class TextItemWidgetFactory
{
public:
    TextItemWidgetFactory();
    virtual ~TextItemWidgetFactory() {}

    virtual TextItemWidget* create(QWidget* parent = 0) = 0;
};

class FilepathWidgetFactory : public TextItemWidgetFactory
{
public:
    TextItemWidget* create(QWidget* parent = 0) override;
};

class OptionWidgetFactory : public TextItemWidgetFactory
{
public:
    OptionWidgetFactory(const QStringList& options);

    TextItemWidget* create(QWidget* parent = 0) override;

private:
    QStringList m_options;
};

class TextWidgetFactory : public TextItemWidgetFactory
{
public:
    TextWidgetFactory(ValidatorFactory* validatorFactory = nullptr);
    virtual ~TextWidgetFactory();

    TextItemWidget* create(QWidget* parent = 0) override;

private:
    ValidatorFactory* m_factory;
};

class DataTypeEditWidgetFactory : public TextItemWidgetFactory, public IValidatableObject
{
public:
    DataTypeEditWidgetFactory(ECFParameter* expression = nullptr);

    TextItemWidget* create(QWidget* parent = 0) override;

    bool isValid() override;
    QString errorMessage() override;

private:
    DataTypeEditWidgetFactoryData m_data;
    ECFParameter* m_expr;
};

}

#endif // TEXTITEMWIDGETFACTORY_H
