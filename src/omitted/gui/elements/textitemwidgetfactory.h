#ifndef TEXTITEMWIDGETFACTORY_H
#define TEXTITEMWIDGETFACTORY_H

#include <QWidget>

#include "textitemwidget.h"
#include "filepathwidget.h"
#include "optionwidget.h"
#include "textwidget.h"
#include "validators/validatorfactory.h"
#include "ivalidatableobject.h"
#include "declarations/datatypeeditwidget.h"

namespace espreso
{

class ECFParameter;

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
    virtual TextItemWidget* create(QWidget* parent = 0) override;
};

class OptionWidgetFactory : public TextItemWidgetFactory
{
public:
    OptionWidgetFactory(const QStringList& options);

    virtual TextItemWidget* create(QWidget* parent = 0) override;

private:
    QStringList m_options;
};

class TextWidgetFactory : public TextItemWidgetFactory
{
public:
    TextWidgetFactory(ValidatorFactory* validatorFactory = nullptr);
    virtual ~TextWidgetFactory();

    virtual TextItemWidget* create(QWidget* parent = 0) override;

private:
    ValidatorFactory* m_factory;
};

class DataTypeEditWidgetFactory : public TextItemWidgetFactory, public IValidatableObject
{
public:
    DataTypeEditWidgetFactory(ECFParameter* expression = nullptr);

    virtual TextItemWidget* create(QWidget* parent = 0) override;

    virtual bool isValid() override;
    virtual QString errorMessage() override;

private:
	DataTypeEditWidgetFactoryData m_data;
    ECFParameter* m_expr;
};

}

#endif // TEXTITEMWIDGETFACTORY_H
